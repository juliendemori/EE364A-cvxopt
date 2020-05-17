%We implement both gradient descent and Newton's method on a minimization
%problem


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

count = 0;
obj = [];
sl = [];

%////////////////////////////////////////
%(a) Gradient Descent

while true;
    grad = A'*(1./(1-A*x)) + 2*x./(1-x.^2);
    if norm(grad) <= eta;
        break
    end
    dx = -grad;
    t = 1;
    while max(A*(x + t*dx)) >= 1 | max(x + t*dx).^2 >=1
        t = bet*t;
    end
    while f(x + t*dx) > f(x) + t*alp*grad'*dx;
        t = t*bet;
    end
    x = x + t*dx;
    obj = [obj f(x)];
    sl = [sl t];
    count = count + 1;
end

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
plot(1:count, sl)
xlabel('Iteration number');
ylabel('Step Length');