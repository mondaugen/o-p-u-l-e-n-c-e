clear;
M=512;
%    N=101;
%    A=eye(M);
m=(0:(M-1));
K=512;
A=[sin(-2*pi*m'*(1:((K/2)-1))/M) cos(-2*pi*m'*(0:K/2)/M)];
A=[eye(M) A];
c=cond(A);
N=size(A,2);
% Make random signal
x=rand(N,1);
x=(x<0.05).*randn(N,1);
b_orig=A*x;
b=b_orig+randn(size(A,1),1);
eig(A*A');

% least-squares solution
x_ls=ols(b,A);

[z_]=bpdn(A,b,50,100);
x_=z_(1:N);
e=norm(A*x_-b);
newplot(figure(1));
figure(1);
plot(1:length(x),x,1:length(x_ls),x_ls);
legend();
newplot(figure(2));
figure(2);
plot(1:length(b_orig),b_orig,'b',1:length(b),b,'r',1:length(b),A*x_,'g');
legend();
newplot(figure(3));
figure(3);
plot(1:length(x_),x_);
