%Implementing Newton's descent


%Generate Data
n = 100;
m = 400;
randn('state', 1);
A = randn(m,n);
f = @(x) - sum(log(1-A*x)) - sum(log(1-x.^2));

%Backtracking parameters
alp = 0.05;
bet = 0.1;

eta = 10^(-2);

%Initializing x
x = zeros(n,1);


%Keeping track of number of iterations and such
count = 0;
obj = [];
st = [];

while true;
    hess = A'*(diag((1./(1-A*x)).^2))*A + diag(1./(1+x).^2);
    grad = A'*(1./(1-A*x)) + 2*x./(1-x.^2);
    dx = - hess\grad;
    lam = -grad'*(dx);
    if (lam^2/2 <= eta)
        break;
    end
    t = 1;
    while max(A*(x + t*dx)) >= 1 | max(x + t*dx).^2 >=1
        t = bet*t;
    end
    while f(x + t*dx) > f(x) + t*alp*grad'*dx;
        t = t*bet;
    end
    x = x + t*dx;
    obj = [obj f(x)];
    st = [st t];
    count = count + 1;
end

%Determining optimum value using cvx
cvx_begin
    cvx_quiet(true);
    variable x(n)
    minimize sum(-log(1-A*x)) + sum(-log(1-x.^2))
cvx_end

p_star = cvx_optval*ones(1, count);

figure(1)
plot(1:count, obj, 'b', 1:count, p_star, 'r');
legend('objective', 'P star');

figure(2)
plot(1:count, st)
xlabel('Iteration number');
ylabel('Step Length');