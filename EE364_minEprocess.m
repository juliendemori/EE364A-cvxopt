%Minimum energy processor speed scheduling

% Data for minimum energy processor speed scheduling.

n = 12;  % number of jobs.
T = 16;  % number of time periods.

Smin = 1;  % min processor speed.
Smax = 4;  % max processor speed.
R = 1;  % max slew rate.

% Parameters in power/speed model.
alpha = 1;
beta = 1;
gamma = 1;


% Job arrival times and deadlines.
A = [1; 3;  4; 6; 9;  9; 11; 12; 13; 13; 14; 15];
D = [3; 6; 11; 7; 9; 12; 13; 15; 15; 16; 14; 16];
% Total work for each job.
W = [2; 4; 10; 2; 3; 2; 3; 2; 3; 4; 1; 4];

%Matrix to help with work constraints
B = zeros(n,T);
for i = 1:n;
   a = A(i);
   d = D(i);
   B(i, a:d) = 1;
end


%The optimization
cvx_begin
    variables theta(n,T) s(T)
    minimize(T + sum(s) + sum(square(s)))
    subject to
       ones(n,1)'*theta == ones(1,T);
        theta(:,1:T) >= 0; 
        s(2:T) - s(1:T-1) <= ones(T-1,1)*R;
        s(2:T) - s(1:T-1) >= -ones(T-1,1)*R;
        s <= Smax;
        s >= Smin;
        B(:,1:T).*theta(:,1:T) == theta(:,1:T);
        for i = 1:n;
            B(i, :)*s >= W(i);
        end
cvx_end

% Plot showing job availability times & deadlines.
figure;
hold on;
scatter(A,[1:n],'k*');
scatter(D+1,[1:n],'ko');
for i=1:n
    plot([A(i),D(i)+1],[i,i],'k-');
end
hold off;
xlabel('time t');
ylabel('job i');


%Minimal energy = 74

%Plots
bar((s*ones(1,n)).*theta',1,'stacked');
