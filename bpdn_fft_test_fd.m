% Test with fourier dictionary.
clear;
M=512;
% Resolution of dictionary 1 is normal fourier dictionary, 2 is twice the
% resolution, etc.
res=4;
m=(0:(M-1));
K=512;
inc=1/res;
A=[sin(-2*pi*m'*(1:inc:((K/2)-1))/M) cos(-2*pi*m'*(0:inc:K/2)/M)];
c=cond(A);
N=size(A,2);
% Make random signal
x=rand(N,1);
x=(x<0.02).*randn(N,1);
b_orig=A*x;
b=b_orig+sqrt(0.5)*randn(size(A,1),1);
eig(A*A');
nu=-70; % signals lower than this are zeroed 
sig=100^(0.8); % sparsity / reconstruction error trade off
% least-squares solution
x_ls=ols(b,A);
la=sig*log(2*size(A,2));
[z_,ep_]=bpdn_fft(res,b,la,1000,1e-2);
display(sprintf('Suboptimality: %f\n',ep_));
x_=z_(1:N);
e=norm(A*x_-b);
newplot(figure(1));
figure(1);
plot(1:length(x),log(10)^(-1)*20*log(abs(x+1e-6)),1:length(x_ls),log(10)^(-1)*20*log(abs(x_ls)),'r');
legend();
newplot(figure(2));
figure(2);
plot(1:length(b_orig),b_orig,'b',1:length(b),b,'r',1:length(b),A*x_,'g');
legend();
newplot(figure(3));
figure(3);
x_.*=(log(10)^(-1)*20*log(abs(x_))) > nu;
x_db=log(10)^(-1)*20*log(abs(x_+1e-6));
plot(1:length(x_),x_db);
