x=100;
dx=@(x_) exp(x)-exp(1-x);
ddx=@(x_) exp(x)+exp(1-x);
L=10;
for l=1:L
    u=-dx(x)/ddx(x)
    t=(1-2*x)/(2*u)
    x=x+t*u
end
