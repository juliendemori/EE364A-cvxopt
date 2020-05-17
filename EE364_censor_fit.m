%We implement the method illustrated in fit-censored-data to find the
%optimal parameters for the linear fit to some known, and some unknown
%data.

% data for censored fitting problem.
randn('state',0);

n = 20;  % dimension of x's
M = 25;  % number of non-censored data points
K = 100; % total number of points
c_true = randn(n,1);
X = randn(n,K);
y = X'*c_true + 0.1*(sqrt(n))*randn(K,1);

% Reorder measurements, then censor
[y, sort_ind] = sort(y);
X = X(:,sort_ind);
D = (y(M)+y(M+1))/2;
y = y(1:M);

cens = K-M;

%Our optimization
cvx_begin
    variables ytot(K) c(n);
    minimize(sum((ytot - (c'*X)').^2));
    subject to
        ytot(M+1:K) >= D;
        ytot(1:M) == y;
cvx_end


%We also obtain a least squares fit on c by simply ignoring the censored
%data

%Make new X matrix containing only the x values up to M
Xnew = X(:,1:M);

cvx_begin
    variables c_unc(n);
    minimize(sum((y - (c_unc'*Xnew)').^2));
cvx_end

RE_cens = sum((c_true - c).^2)/sum(c_true.^2)
RE_uncens = sum((c_true - c_unc).^2)/sum(c_true.^2)