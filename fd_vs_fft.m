clear;
M=512;
% Resolution of dictionary 1 is normal fourier dictionary, 2 is twice the
% resolution, etc.
res=4;
L=M*res;
m=(0:(M-1));
k_real=(0:(L/2));
k_imag=(1:(L/2-1));
A=[cos(-2*pi*m'/L*k_real) sin(-2*pi*m'/L*k_imag)];
x=randn(M,1);
X_fd=A'*x;
X_fft_=fft(x,L);
X_fft=[real(X_fft_(k_real+1));
       imag(X_fft_(k_imag+1))];
plot(1:length(X_fd),X_fd,1:length(X_fft),X_fft);
norm(X_fd-X_fft)

