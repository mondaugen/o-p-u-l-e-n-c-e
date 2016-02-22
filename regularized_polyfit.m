N=50;
N1=10;
N2=40;
x=1:N;
p=3;
gam=(exp(-1:0.01:1)-1)/(exp(1)-1);
%gam=[0 1]
z=zeros(length(gam),N);
y=[ zeros(1,length(1:N1)),...
    randn(1,length((N1+1):(N2))),...
    zeros(1,length((N2+1):N))    ];
x_=x(N1:N2);
A=x_'.^(0:p);
for k=1:length(gam)
    R=diag([6*gam(k) zeros(1,p)]);
    om=linsolve((A'*A + R'*R),A'*y(N1:N2)');
    z(k,:)=polyval(flipud(om),x);
end
plot(x,y,'b',x,z,'g;a;');
