%We solve the sparse index tracking problem

% data generation for sparse index tracking problem
% generates n, c, rbar, Sigma

rand('state',1);
randn('state',1);

n = 500; % number of stocks
k = 20;  % number of factors (not described or needed in problem statement)

% generate rbar, Sigma
F = 1.5*rand(n,k)-0.5; % factor matrix, entries uniform on [-0.5,1]
Sigma = 0.5*F*diag(rand(k,1))*F' + diag(0.1*rand(n,1));
mu_rf = 0.2;                          % risk-free return (weekly, 5% annual)
SR = 0.4;                             % Sharpe ratio
rbar = mu_rf + SR*sqrt(diag(Sigma));  % expected return
c = 5 + exp(2*randn(n,1));            % market capitalization (index weights)

cvx_begin
    variable w(n);
    minimize(sum(abs(w)));
    subject to 
        (c'*(Sigma + rbar*rbar')*c + quad_form(w, Sigma + rbar*rbar') - 2*w'*(Sigma + rbar*rbar')*c) <= 0.01*c'*(Sigma + rbar*rbar')*c;
cvx_end

%Count the number of nonzero entries
cardw = sum(abs(w) > 0.01);

%cardw = 52;