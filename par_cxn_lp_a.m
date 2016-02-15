function [costs]=par_cxn_lp_a(par_path,cxn_path,a='>=1',lamb=0,mu=1)
    % PAR_CXN_LP
    %
    % Solve partial connections via linear programming
    % par_path - path to file containing partials
    % cxn_path - path to file where connections will be saved
    % a        - adaptation method:
    %            '>=1' : if unequal number of in-partials and out-partials,
    %                    allow partials with multiple connections
    %            '<=1' : if unequal number of in-partials and out-partials,
    %                    allow partials with no connections
    % lamb     - weight of last costs incorporated (should be between 0 and 1)
    %            dampen discount of having similar last cost. lamb = 1 means
    %            having cost equal to last cost, multiplies current cost by 0,
    %            lamb = 0 means there's no benefit to having a similar last cost
    % mu       - if cost is greater than this value, consider connection as
    %            erroneous
    
    % Weight vector for amp,freq,damp,dfreq costs resp.
    L_PARTIAL_RECORD=4;
    w=[0;1;0;0];
    lastData=[];
    f=fopen(par_path,'r');
    g=fopen(cxn_path,'w');
    costs={};
    k=1;
    lastCosts=[];
    while ~feof(f)
        dtime=fread(f,1,'double');
        if(length(dtime)==0)
            break;
        end
        ndata=fread(f,1,'uint32');
        if(length(ndata)==0)
            break;
        end
        data=fread(f,ndata,'double');
        data=reshape(data,[length(data)/L_PARTIAL_RECORD,L_PARTIAL_RECORD]);
        % use only frequency data
        data=data(:,2);
        if (length(lastData)==0)
            lastData=data;
            continue;
        end
        m=size(data,1);
        n=size(lastData,1);
%        C=zeros(m,n);
        C=pt_build_cost_matrix(lastData(:),data(:),@divergence);
        if (length(lastCosts)==0)
            lastCosts=ones(1,n);
        end
        disc=(1 - lamb*exp(-(log(C./(ones(m,1)*lastCosts)).^2)));
        C=C.*disc;
        %C.*=(lastCosts.^(lamb));
        costs{k}=C;
        k+=1;
        % Transpose because we want to stack rows into a column
        C=C';
        c=C(:);
        size(c);
        % x is column vector where x((j-1)*m+i) = weight of connection (i,j)
        % cost of all connections into j sum to K
        % K1 is lower bound on number of connections into an output partial
        % K2 is lower bound on number of connections out of an input partial
        % K3 is upper bound on number of connections into an ouput partial
        % K4 is upper bound on number of connections out of an input partial 
        if strcmp(a,'>=1') == 1
            % Allow partials with more than one connection
            K1=1;
            K2=1;
            if (m > n)
                % more output partials than input, some input partials must
                % have more than 1 connection
                K3=1;
                K4=ceil(m/n);
            elseif (n > m)
                % more input partials than output, some output partials must
                % have more than 1 connection
                K3=ceil(n/m);
                K4=1;
            else
                % each input partial is connected to one output partial
                K3=1;
                K4=1;
            end
        elseif strcmp(a,'<=1') == 1
            % Allow partials with no connections
            K3=1;
            K4=1;
            if (m > n)
                % more output partials than input, some output partials have no
                % connection
                K1=0;
                K2=1;
            elseif (n > m)
                % more input partials than output, some input partials 
                % have no connection
                K1=1;
                K2=0;
            else
                % each input partial is connected to one output partial
                K1=1;
                K2=1;
            end
        else
            error('Bad adaptation method.');
        end
        % bounds on weight variables
        LB=0;
        UB=Inf;
        A1=zeros(m,n*m);
        for i_=1:m
            for j_=1:(n*m)
                A1(i_,j_)=(floor((j_-1)/n)+1)==i_;
            end
        end
        b1=K1*ones(m,1);
        A2=zeros(n,n*m);
        for i_=1:n
            for j_=1:(n*m)
                A2(i_,j_)=(mod((j_-1),n)+1)==i_;
            end
        end
        A2;
        b2=K2*ones(n,1);
        A3=A1;
        A4=A2;
        b3=K3*ones(m,1);
        b4=K4*ones(n,1);
        A=[A1;A2;A3;A4];
        b=[b1;b2;b3;b4];
        ctype1=zeros(length(b1),1);
        ctype1(:,:)='L';
        ctype2=zeros(length(b2),1);
        ctype2(:,:)='L';
        ctype3=zeros(length(b3),1);
        ctype3(:,:)='U';
        ctype4=zeros(length(b4),1);
        ctype4(:,:)='U';
        ctype=[ctype1;ctype2;ctype3;ctype4];
        ctype=char(ctype);
        vartype=zeros(m*n,1);
        vartype(:,:)='C';
        vartype=char(vartype);
        [xopt,fmin,errnum,extra]=glpk(c,...
                                      A,...
                                      b,...
                                      LB*ones(n*m,1),...
                                      UB*ones(n*m,1),...
                                      ctype,...
                                      vartype);
        fmin;
        errnum
%        xopt_out=(C' < mu).*reshape(xopt,[n,m])';
        % omitted bad costs still influence next costs
        lastCosts=sum(C'.*reshape(xopt,[n,m])',2)';
        %small=1e-10;
        %lastCosts=lastCosts'+small;
        %lcsum=sum(lastCosts);
        %lastCosts/=lcsum;
        % Output file is <uint32_t> number of rows m
        %                <uint32_t> number of columns n
        %                <double * m * n> x vector as solution for optimal weights
%        xopt_out=xopt_out';
%        xopt_out=xopt_out(:);
        fwrite(g,m,'uint32');
        fwrite(g,n,'uint32');
        fwrite(g,xopt,'double');
        lastData=data;
    end
    fclose(f);
    fclose(g);
end
%
function [d]=divergence(x,y)
    %d=exp(abs(log(x) - log(y)));
    %d=(log(x)-log(y))^2;
    d=(x.-y).^2;
end
