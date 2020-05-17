function [ xstar, vstar, count ] = equal_constraint_newton( A, b, c, xi)
%Function minimizes a a linear function -(sum of logs) subject to equality
%constraints


[m,n] = size(A);

%Algorithm params
alpha = 0.05;
beta = 0.3;

F = null(A);
cf = c'*F;


%Defining the function
f = @(z) cf*z + c'*xi - sum(log(F*z + xi));

%Defining the tolerance
tol = 10^(-6);

%Initializing vector
z = zeros(n-m,1);

%Keeping track
count = 0;
vals = [];


while true
   g = 1./(F*z + xi);
   grad = cf' - F'*g; 
   hess = F'*diag(g.^2)*F;
   dz = - hess\grad;
   lamsq = - dz'*grad;
   if lamsq/2 <= tol;
      break; 
   end
   t = 1;
   while min((F*z + t*F*dz + xi)) < 0
      t = beta*t; 
   end
   while f(z+t*dz) > (f(z) + alpha*t*grad'*dz)
      t = beta*t;
   end
   z = z + t*dz;
   count = count + 1;
   vals = [vals; lamsq/2];
end

%Have to compute the value of x and nu from z
xstar = F*z + xi;
gradxs = c'*xstar - 1./(xstar);
vstar = -inv(A*A')*A*gradxs;
end

