%We upload the relevant problem data, and proceed to find the most sparse
%descriptions of A and B that are plausible for the general LDS with
%Gaussian noise. 

% Generate a random, sparse system.
% Initialize both random generators, because sprandn uses both.
randn('state', 118); rand('state', 118);
n = 8;
A = full(sprandn(n, n, 0.2));
A = 0.95*A/max(abs(eig(A))); % make A stable.
m = 4;
B = full(sprandn(n, m, 0.3));
T = 100;
W = eye(n);
Whalf = sqrtm(W);

us = 10*randn(m, T-1); % input.
ws = Whalf*randn(n, T); % noise process.

xs = zeros(n, T);
xs(:, 1) = 50*randn(n, 1); % initial x.

% Simulate the system.
for t = 1:T-1
    xs(:,t + 1) = A*xs(:,t) + ws(:,t) + B*us(:,t);
end

%Generate a matrix of x(t) values so that all the constraints can be
%written in one cvx_line
Xt = zeros(n, T-2);
X = zeros(n, T-2);
U = zeros(m, T-2);
for i = 1:T-2;
   X(:,i) = xs(:,i);
   Xt(:,i) = xs(:,i+1);
   U(:,i) = us(:,i);
end


%Optimize
cvx_begin
    variables At(n,n) Bt(n,m) t v;
    minimize(sum(abs(At(:))) + sum(abs(Bt(:))));
    subject to
        sum(sum((Whalf^(-1)*(Xt - At*X - Bt*U)).^2)) <= n*(T-1) + 2*sqrt(2*n*(T-1));
cvx_end


%We want to count false positives and false negatives in At and Bt
afp = 0;
afn = 0;
bfp = 0;
bfn = 0;
for i = 1:n;
    for j = 1:n;
        if abs(At(i,j)) >= 0.01 && abs(A(i,j)) <= 0.01;
            afp = afp + 1;
        end
        if abs(At(i,j)) <= 0.01 && abs(A(i,j)) >= 0.01;
            afn = afn + 1;
        end
    end
end


for i = 1:n;
    for j = 1:m;
        if abs(Bt(i,j)) >= 0.01 && abs(B(i,j)) <= 0.01;
            bfp = bfp + 1;
        end
        if abs(Bt(i,j)) <= 0.01 && abs(B(i,j)) >= 0.01;
            bfn = bfn + 1;
        end
    end
end