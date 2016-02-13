phi=2*pi*rand(1)
w=2*pi*rand(1)
t=(0:0.01:1);
N=length(t);
x=sin(w*t+phi);
dx=w*cos(w*t+phi);
plot(t,x,t,dx);
n=randi(N)
w_=((dx(n)^2)/(1 - x(n)^2))^(1/2)
phi_1=asin(x(n))-w_*t(n);
phi_2=fmod(pi-asin(x(n)),pi)-w_*t(n);
while (phi_1 < 0)
    phi_1 += pi;
end
while (phi_2 < 0)
    phi_2 += pi;
end
phi_1
phi_2
