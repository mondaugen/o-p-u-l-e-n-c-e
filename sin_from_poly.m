w=3*pi;
A=1;
phi=2*pi*0.25;
x=(-2:0.01:2);
y=A*cos(w*x+phi);
phis=zeros(length(x),1);
newplot(figure(1));
figure(1);
last_phi=0;
for n=2:(length(x)-1)
    x0=x(n);
    x_=x((n-1):(n+1))-x0;
    y_=y((n-1):(n+1));
    p_=polyfit(x_,y_,length(x_)-1);
%    w_=(-2*p_(1)/p_(3))^(1/2)
%    phi_=atan2(-p_(2),(p_(1)*w_))-w_*x0
%    A_=p_(1)/cos(w_*x0+phi_)
    %plot(x,y,'b',x,A_*cos(w_*x+phi_),'g');
    p(1)=-A*(w^2)*cos(w*x0+phi)/2;
    p(2)=-A*w*sin(w*x0+phi);
    p(3)=A*cos(w*x0+phi);
    w_=(-2*p(1)/p(3))^(1/2)
    phi_=atan2(-p(2),(p(1)*w_))-w_*x0;
    A_=p(1)/cos(w_*x0+phi_)
    phis(n)=phi_-last_phi;
    last_phi=phi_;
%    plot(x,y,'b',x,polyval(p,x-x0),'g',x,polyval(p_,x-x0),'r');
    axis([x(1) x(end) -1 1]);
%    sleep(0.5);
end
newplot(figure(2));
figure(2);
plot(x,phis);
