
%//////////////////////////////////////////////
%Part (b)

%data
rand('state', 5);
randn('state', 5);
n=20;
pbar = ones(n,1)*.03+[rand(n-1,1); 0]*.12;
S = randn(n,n);
S = S'*S;
S = S/max(abs(diag(S)))*.2;
x_unif = ones(n,1)/n;
S(:,n) = zeros(n,1);
S(n,:) = zeros(n,1)';
one = ones(n,1);
%Vectors to store the data that we will collect for various values of mu
mu = 0.1:0.1:1;
nmu = length(mu);
data1 = zeros(nmu,2);
data2 = zeros(nmu,2);



%Second: x >= 0
for i = 1:nmu;
cvx_begin
    variable x_t(n);
    minimize(-pbar'*x_t + mu(i) * x_t'*S*x_t);
    subject to
        one' * pbar == 1;
        x_t >= 0;
cvx_end
data1(i,1) = x_t'*S*x_t;
data1(i,2) = pbar'*x_t;

%Third: limit on total short position
cvx_begin
    variables x_th(n);
    minimize(-(1-mu(i))*pbar'*x_th + mu(i)*x_th'*S*x_th);
    subject to
        one' * pbar == 1;
        one'* max(-x_th, zeros(n,1)) <= 0.5;
cvx_end
data2(i,1) = x_th'*S*x_th;
data2(i,2) = pbar'*x_th;
end


plot(data1(:, 1), data1(:,2),'red', data2(:, 1), data2(:, 2), 'blue');
xlabel('Standard Deviation');
ylabel('Mean return');
title('Optimal trade-off curves for two cases');
legend('x >= 0', 'limit on total short position');