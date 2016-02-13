% Use linear programming to solve partial connection paths (?)
% A(i,j) is cost of going from partial peak i to partial peak j
N=4; % number of variables
C = [ 1, 3;
      2, 2 ];
c=C(:);
A=[ 1 1 0 0;
    0 0 1 1 ];
[xopt,fmin,errnum,extra]=glpk(c,A,ones(2,1),zeros(N,1),ones(N,1),'SS','CCCC');
xopt
