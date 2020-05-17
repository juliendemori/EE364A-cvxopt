%Trying to efficiently solve the regularized least squares problem.

k = 100;
n = 4000;
A = randn(k, n);
bstart = randn(k,1);
b = A'*bstart;

%Creating sparse tridiagonal matrix
bla = [-ones(n,1), 3*ones(n,1), -ones(n,1)];
S = spdiags(bla, -1:1, n, n);

B = A';
C = A;


t = cputime; 
x = (S+A'*A)\b;
e = cputime - t

[h, w] = size(B);
B_new = [];
for i = 1:w;
   B_new = [B_new B(:,i)]; 
end


%Faster method, exploiting matrix structure
tnew = cputime;
bob = S\b;
boob = S\B_new;
y = (speye(k,k) + C*boob)\(C*bob);
xnew = S\(b - B*y);
eff = cputime - tnew