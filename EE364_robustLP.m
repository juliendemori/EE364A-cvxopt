%We solve the single linear program that we have reformulated the robust LP
%into.

%Initializing A,b, and cnom
rand('seed', 0);
A = rand(30,10);
b = rand(30,1);
c_nom = 1+rand(10,1); % nominal c values

n = 10;

%Constructing F and g from cnom
F = [eye(n,n); -eye(n,n); ones(n,1)'/n; -ones(n,1)'/n];
g = [1.25*c_nom; -0.75*c_nom; ones(n,1)'*c_nom/n*1.1; -ones(n,1)'*c_nom/n*0.9];


cvx_begin
    variables x(n) lam(2*n + 2);
    minimize(lam'*g);
    subject to
        A*x >= b;
        F'*lam == x;
        lam >= 0;
cvx_end

display(x)
display(lam)


%Simplest case: C = {c_nom}

cvx_begin
    variable x_new(n);
    minimize(c_nom'*x_new);
    subject to
        A*x_new >= b;
cvx_end

display(x_new)