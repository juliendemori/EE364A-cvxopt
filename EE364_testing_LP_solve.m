%Testing
m = 10;
n = 50;

A = randn(m-1,n);
A = [A; rand(1,n)];
xi = 10^(-5) * rand(n,1);
b = A*xi;
c = randn(n,1);

[xstar, vstar, count] = equal_constraint_newt_elim(A,b,c,xi);