%Now we want to compute the minimum-risk portfolios with the same expected
%return as the uniform portfolio, for a few different constraint cases.
%All cases have the constraint 1'x = 1

%/////////////////////////
%Part (a)


%data
rand('state', 5);
randn('state', 5);
n=20;
pbar = ones(n,1)*.03+[rand(n-1,1); 0]*.12;
S = randn(n,n);
S = S'*S;
S = S/max(abs(diag(S)))*.2;
x_unif = ones(n,1)/n;

one = ones(n,1);


%First: no additional contraint
cvx_begin
    variable x_o(n);
    minimize(x_o'*S*x_o);
    subject to
        one' * pbar == 1;
        pbar'*(x_o) == pbar'*x_unif;
cvx_end
disp(x_o)


%Second: x >= 0
cvx_begin
    variable x_t(n);
    minimize(x_t'*S*x_t);
    subject to
        one' * pbar == 1;
        pbar'*x_t == pbar'*x_unif;
        x_t >= 0;
cvx_end
disp(x_t)


%Third: limit on total short position
cvx_begin
    variables x_th(n);
    minimize(x_th'*S*x_th);
    subject to
        one' * pbar == 1;
        pbar'*x_th == pbar'*x_unif;
        one'* max(-x_th, zeros(n,1)) <= 0.5;
cvx_end
disp(x_th)



