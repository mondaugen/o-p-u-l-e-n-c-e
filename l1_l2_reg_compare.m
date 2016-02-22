clear;
% L1_L2_REG_COMPARE
% Compare line fitting minimizing L2 and L1 error
m=10.2;
b=1;
x=(0:0.001:1)';
M=length(x);
sigma1=0.2;
sigma2=10;
y=m*x+b+randn(M,1)*sigma1;
% corrupt random indices
n_=randi([1 M/2],floor(M*0.1),1);
y(n_)+=rand(floor(M*0.1),1)*sigma2;
A=[ones(M,1) x];
r=ols(y,A);
z=l1_norm_aprx(A,y,1000);
r_l1=z(1:2);
plot(x,y,x,r(1)+x*r(2),x,r_l1(1)+x*r_l1(2));
