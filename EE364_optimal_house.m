%We solve part (c) of the optimal house betting with the data provided. 

%Initializing
n = 5;
m = 5;

p = [0.5; 0.6; 0.6; 0.6; 0.2];
q = [10; 5; 5; 20; 10];
A = [1, 1, 0, 0 , 0; 0, 0 , 0, 1, 0; 1, 0, 0, 1, 1; 0, 1, 0, 0, 5; 0, 0, 1, 0, 0];

cvx_begin
    variables x(n) t;
    minimize(-p'*x + t);
    subject to
        x <= q;
        x >= 0;
        A*x <= t*ones(n,1);
cvx_end


%Simple case x = q
cvx_begin
    variables t_new;
    minimize(-p'*q + t_new);
    subject to
        A*q <= t_new*ones(n,1);
cvx_end

