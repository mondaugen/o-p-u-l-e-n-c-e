a=[1;1];
[xx,yy]=meshgrid(0.01:0.01:0.99,0.01:0.01:0.99);
fx=-log(1-a(1)*xx)-log(1-a(2)*yy)-log(xx)-log(yy);
%df=@(z) [1/(1-z(1))-1/z(1);1/(1-z(2))-1/z(2)];
df=@(x_,y_) [1/(1-x_)-1/x_;1/(1-y_)-1/y_];
ddf=@(x_,y_) diag([1/((1-x_)^2)+1/(x_^2),1/((1-y_)^2)+1/(y_^2)]);
contour(xx,yy,fx);
z=[0.9;0.7];
ep=1e-3;
lamb2=df(z(1),z(2))'*ddf(z(1),z(2))*df(z(1),z(2));
pa=[];
while lamb2/2 > ep
    pa=[pa z];
    u=-(ddf(z(1),z(2))^(-1))*df(z(1),z(2))
    lamb2=df(z(1),z(2))'*ddf(z(1),z(2))*df(z(1),z(2));
    z=z+u
end
hold on
plot(pa(1,:),pa(2,:));
hold off
