function [C] = pt_build_cost_matrix(frame0,frame1,cost=@(x,y) (x.-y).^2)
% Cost of connection from frame0(n) to frame1(m) is in C(m,n)
M=length(frame1);
N=length(frame0);
C=cost(frame1(:)*ones(1,N),ones(M,1)*frame0(:).');
