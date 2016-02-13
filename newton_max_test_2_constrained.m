clear;
a=[1;1];
A=[3 -1];
[xx,yy]=meshgrid(0.1:0.1:0.9,0.1:0.1:0.9);
fx=-log(1-a(1)*xx)-log(1-a(2)*yy)-log(xx)-log(yy);
%df=@(z) [1/(1-z(1))-1/z(1);1/(1-z(2))-1/z(2)];
df=@(x_,y_) [1/(1-x_)-1/x_;1/(1-y_)-1/y_];
ddf=@(x_,y_) diag([1/((1-x_)^2)+1/(x_^2),1/((1-y_)^2)+1/(y_^2)]);
contour(xx,yy,fx);
z=[.1;.3];
ep=1e-3;
lamb2=df(z(1),z(2))'*ddf(z(1),z(2))*df(z(1),z(2));
pa=[];
while lamb2/2 > ep
    H=ddf(z(1),z(2));
    tmp=[-df(z(1),z(2));
         0];
    KKT=[H  A'
         A 0];
    tmp=linsolve(KKT,tmp);
    u=tmp(1:2);
    w=tmp(3)
    u=u(:)
    pa=[pa z];
    lamb2=u'*ddf(z(1),z(2))*u;
    z=z+u
end
hold on
plot([0 1/3],[0 1]);
plot(pa(1,:),pa(2,:),'r');
hold off
