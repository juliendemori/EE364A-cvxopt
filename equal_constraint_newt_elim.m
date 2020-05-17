function [ xstar, vstar, count ] = equal_constraint_newt_elim( A, b, c, xi)
%Function minimizes a a linear function -(sum of logs) subject to equality
%constraints, using KKT block elimination


[m,n] = size(A);

%Algorithm params
alpha = 0.1;
beta = 0.3;


%Defining the function
f = @(x) c'*x - sum(log(x));

%Defining the tolerance
tol = 10^(-6);

%Keeping track
count = 0;
vals = [];

%Initializing x
x = xi;

while true
   g = 1./(x)
   grad = c - g;
   hess = diag(g.^2)
   %Block elimination solve
   boo = hess\(A');
   baa = hess\grad;
   S = -A*boo;
   w = S\(A*baa);
   dx = hess\(-A'*w - grad);
   lamsq = - dx'*grad;
   if lamsq/2 <= tol;
      break; 
   end
   t = 1;
   while min((x + t*dx)) < 0
      t = beta*t; 
   end
   while f(x+t*dx) > (f(x) + alpha*t*grad'*dx)
      t = beta*t;
   end
   x = x + t*dx;
   count = count + 1;
   vals = [vals; lamsq/2];
end



%Have to compute the value of xstar
xstar = x;
gradxs = c'*xstar - 1./(xstar);
vstar = -inv(A*A')*A*gradxs;

%Check that KKT conditions are satisfied
[hess, A'; A, zeros(m,m)] * [xstar; vstar] + [gradxs; zeros(length(b),1)]
end
