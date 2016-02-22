function [x,x_,x_ls,c,p,e]=l1_norm_approx_test()
    DEBUG=1;
    
    %% Build dictionary of cosine and sine wavelets 
    %sigma=0.05;
    %tau=0.002;
    %omega0=5;
    %K=51;
    %L=10;
    %M=51;
    %m=(0:(M-1))';
    %k=(0:(K-1));
    %m=m*ones(1,K);
    %ts=m./M;
    %k=ones(M,1)*k;
    %Ad=exp(-(ts.-k*tau).^2/(sigma^2));
    %A=Ad;
    %for l=(1:L)
    %    Ac=cos(l*omega0*(ts.-k*tau));
    %    As=sin(l*omega0*(ts.-k*tau));
    %    Ac.*=Ad;
    %    As.*=Ad;
    %    A=[A Ac As];
    %end

    M=128;
%    N=101;
%    A=eye(M);
    m=(0:(M-1));
    K=128;
    A=[sin(-2*pi*m'*(1:((K/2)-1))/M) cos(-2*pi*m'*(0:K/2)/M)];
    A=[eye(M) A];
    c=cond(A);
    N=size(A,2);
    % Make random signal
    x=rand(N,1);
    x=(x<0.05).*randn(N,1);
    b=A*x;

    % least-squares solution
    x_ls=ols(b,A);

    [z_,p]=l1_norm_aprx(A,b,1000);
    x_=z_(1:N);
    e=norm(A*x_-b);
end;
