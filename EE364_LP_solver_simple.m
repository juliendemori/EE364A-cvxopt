%LP solver with simple inequality constraint using my previous solver

%Parameters
tee = 1;
mu = 15;

%Generate problem data.
m = 100;
n = 200;
A = randn(m,n);
A = [A; rand(1,n)];
xi = rand(n,1);
b = A*xi;
c = randn(n,1);

%Algorithm params
alpha = 0.01;
beta = 0.4;

F = null(A);
%T = F*F';
cf = c'*F;
xc = c'*xi;

%Defining the function
f = @(z) (cf*z + c'*xi) - sum(log(F*z + xi));

%New function
%bla = @(x) (tee)*

%Defining the tolerance
tol_lam = 10^(-6);
tol_tee = 10^(-3);



xnew = xi;
while true
   %Centering step
   z = zeros(n-(m+1),1);
   while true
       g = 1./(F*z + xnew);
       grad = (tee)*cf' - F'*g; 
       hess = F'*diag(g.^2)*F;
       dz = - hess\grad;
       lamsq = - dz'*grad;
       if lamsq/2 <= tol_lam;
          break; 
       end
       t = 1;
       while (F*z + t*F*dz + xnew) <= 0
          t = beta*t; 
       end

       while f(z + t*dz) > f(z) + alpha*t*grad'*dz;
       z = z + t*dz;
       end
   end
   xnew = F*z + xnew
   if n/tee < tol_tee;
       break;
   end
   tee = tee*mu;
end