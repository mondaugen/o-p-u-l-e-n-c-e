function []=par_cxn_lp_c(par_path,cxn_path)
    % PAR_CXN_LP
    %
    % Solve partial connections via linear programming
    % par_path - path to file containing partials
    % cxn_path - path to file where connections will be saved
    %
    
    % Weight vector for amp,freq,damp,dfreq costs resp.
    L_PARTIAL_RECORD=4;
    w=[0;1;0;0];
    lastData=[];
    f=fopen(par_path,'r');
    g=fopen(cxn_path,'w');
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
        if (length(lastData)==0)
            lastData=data;
            continue;
        end
        m=size(data,1);
        n=size(lastData,1);
        C=zeros(m,n);
        for i=1:m
            for j=1:n
                C(i,j)=[divergence(lastData(j,1),data(i,1)),...
                        divergence(lastData(j,2),data(i,2)),...
                        divergence(lastData(j,3),data(i,3)),...
                        divergence(lastData(j,4),data(i,4))]*w;
            end
        end
        C;
        c=C(:);
        size(c);
        % x is column vector where x((j-1)*m+i) = weight of connection (i,j)
        % cost of all connections into j sum to K
        % currently K is just a constant, but could be altered to represent the
        % "importance" of a particular partial
        K=1;
        LB=0;
        UB=Inf;
        A=zeros(m,n*m);
        for i_=1:m
            for j_=1:(n*m)
                A(i_,j_)=((floor((j_-1)/n)+1)==i_)*C(i_,mod(j_-1,n)+1);
            end
        end
        A;
        b=K*ones(m,1);
        %A=zeros(n,n*m);
        %for i_=1:n
        %    for j_=1:(n*m)
        %        A(i_,j_)=(mod((j_-1),n)+1)==i_;
        %    end
        %end
        %A;
        %b=K*ones(n,1);
        ctype=zeros(length(b),1);
        ctype(:,:)='S';
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
        fmin
        errnum
        % Output file is <uint32_t> number of rows m
        %                <uint32_t> number of columns n
        %                <double * m * n> x vector as solution for optimal weights
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
    d=abs(log(x) - log(y));
    %d=(x-y)^2;
end
