function [ history, xstar ] = solve_simple_LP( A, b, c, xi )
%This function computes the optimal solution of the simple constrained
%problem when we start on a feasible point

%Some relevant parameters
t = 0.1;
mu = 25;

[m, n] = size(A);

tol = 10^(-3);
xpred = xi;

%Matrix of interesting schtuff
history = [];

while true
   [xstar, vstar, count] = equal_constraint_newt_elim(A, b, t*c, xpred);
   temp = zeros(2,1);
   temp(1,1) = count;
   temp(2,1) = n/t;
   history = [history temp];
   xpred = xstar;
   if (n/t) < tol;
       break;
   end
   t = mu*t;
end

%xstar = xpred;
end

