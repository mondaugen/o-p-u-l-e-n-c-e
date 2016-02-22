function [costs,errnums]=par_cxn_lp_b(par_path,
                              cxn_path,
                              alpha=0.5,
                              beta=0.5,
                              p=2,
                              lamb=0.95,
                              mu=1)
    % PAR_CXN_LP_B
    %
    % Solve partial connections via linear programming
    % A linear prediction of the next node is incorporated into the cost function.
    %
    % par_path - path to file containing partials
    % cxn_path - path to file where connections will be saved
    % alpha    - influence of difference between current node and next node in
    %            cost function
    % beta     - influence of difference between current prediction of next node
    %            and next node in cost function
    % p        - order of predictor
    % lamb     - forgetting factor for RLS linear predictor
    % mu       - if cost is greater than this value, consider connection as
    %            erroneous
    %
    % NOTE: alpha and beta best add up to 1
    % the only mode available is
    %            a = '<=1' : if unequal number of in-partials and out-partials,
    %                    allow partials with no connections
    
    % Weight vector for amp,freq,damp,dfreq costs resp.
    a='<=1';
    L_PARTIAL_RECORD=4;
    w=[0;1;0;0];
    data_k0={};
    f=fopen(par_path,'r');
    g=fopen(cxn_path,'w');
    costs={};
    k=1;
    costs_k0=[];
    D_k0=[];
    % RLS initializer
    delta=0.0001;
    while ~feof(f)
        dtime=fread(f,1,'double');
        if(length(dtime)==0)
            break;
        end
        ndata_k1=fread(f,1,'uint32');
        if(length(ndata_k1)==0)
            break;
        end
        data_k1=fread(f,ndata_k1,'double');
        data_k1=reshape(data_k1,[length(data_k1)/L_PARTIAL_RECORD,L_PARTIAL_RECORD]);
        % use only frequency data_k1
        data_k1=data_k1(:,2);
        m=size(data_k1,1);
        if (length(data_k0)==0)
            % There are no data for paths to come from, initialize and continue
            % again from the top of the loop
            data_k0=cell(m,1);
            W_k0=cell(m,1);
            P_k0=cell(m,1);
            for m_=1:m    
                data_k0{m_}=[data_k1(m_);zeros(p-1,1)];
                W_k0{m_}=[1;zeros(p-1,1)];
                P_k0{m_}=eye(p)/delta;
            end
            continue;
        end
        n=size(data_k0,1);
        W_k1=cell(m,n);
        P_k1=cell(m,n);
        E=zeros(m,n);
        data_k0_=zeros(n,1);
        for n_=1:n
            data_k0_(n_)=data_k0{n_}(1);
            for m_=1:m
                [W_k1{m_,n_},P_k1{m_,n_},E(m_,n_)]=rls_1(data_k1(m_),
                                                         data_k0{n_},
                                                         W_k0{n_},
                                                         P_k0{n_},
                                                         lamb);
            end
        end
        C1=pt_build_cost_matrix(data_k0_(:),data_k1(:),@divergence);
        C=alpha*C1+beta*abs(E);
        D=pt_build_cost_matrix(data_k0_(:),data_k1(:),@difference);
        if (length(costs_k0)==0)
            costs_k0=zeros(1,n);
        end
        if (length(D_k0)==0)
            D_k0=zeros(1,n);
        end
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
        errnums{k}=errnum;
        X=(abs(D) < mu).*reshape(xopt,[n,m])';
        W_k0=cell(m,1);
        data_k0_new=cell(m,1);
        data_k0;
        for m_=1:m
            if any(X(m_,:) > 0)
                n_=find(X(m_,:) > 0)(1);
                W_k0{m_}=W_k1{m_,n_};
                P_k0{m_}=P_k1{m_,n_};
                data_k0_new{m_}=[data_k1(m_);data_k0{n_}(1:(p-1),:)];
            else
                W_k0{m_}=[1;zeros(p-1,1)];
                P_k0{m_}=eye(p)./delta;
                data_k0_new{m_}=[data_k1(m_);zeros(p-1,1)];
            end
        end
        data_k0=data_k0_new;
        costs_k0=sum(C'.*X,2)';
        D_k0=sum(D.*X,2)';
        %small=1e-10;
        %costs_k0=costs_k0'+small;
        %lcsum=sum(costs_k0);
        %costs_k0/=lcsum;
        % Output file is <uint32_t> number of rows m
        %                <uint32_t> number of columns n
        %                <double * m * n> x vector as solution for optimal weights
        xopt_out=X';
        xopt_out=xopt_out(:);
        fwrite(g,m,'uint32');
        fwrite(g,n,'uint32');
%        fwrite(g,xopt,'double');
        fwrite(g,xopt_out,'double');
    end
    fclose(f);
    fclose(g);
end
%
function [d]=divergence(x,y)
    d=(x.-y).^2;
end
%
function [d]=difference(x,y)
    d=(x.-y);
end
