%This program implements a general LP solver with the constraint x >= 0

%Dimension parameters
n = 50;
m = 10;

%Initializing data
%Have to add additional row and column to x1 and A respectively to solve
%phase 1 in part (b) formulation
A = randn(m,n);
A1 = [A, -sum(A,2)/2];
c1 = zeros(n+1,1);
c1(n+1) = 1;
x1 = randn(n+1,1);
x1(n+1) = 2 - min(x1);
z = x1 + (x1(n+1)+1)*ones(n+1,1);
b1 = A1*z;

%We now implement phase1 to determine starting point for phase2
[history, xf] = solve_simple_LP(A1, b1, c1, z);

%Now implement phase2
c2 = randn(n,1);
[history, xstar] = solve_simple_LP(A, b1, c2, xf(1:n));