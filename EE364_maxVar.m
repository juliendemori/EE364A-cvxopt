%Problem 13.4
%We seek to obtain the covariance matrix that gives us the maximum possible
%variance for the portfolio

n = 4;
x = [0.1; 0.2; - 0.05; 0.1]
cvx_begin
    variable S(n,n);
    S == semidefinite(n);
    maximize(x' * S * x);
    subject to
        S(1,1) == 0.2;
        S(1,2) >= 0;
        S(1,3) >= 0;
        S(2,1) >= 0;
        S(2,2) == 0.1;3
        S(2,3) <= 0;
        S(2,4) <= 0;
        S(3,1) >= 0;
        S(3,2) <= 0;
        S(3,3) == 0.3;
        S(3,4) >= 0;
        S(4,2) <= 0;
        S(4,3) >= 0;
        S(4,4) == 0.1;
cvx_end
  
disp(S)
  
diags = [0.2, 0.1, 0.3, 0.1];
S_diag = diag(diags)
p_diag = x'*S_diag*x