%We solve a specific case of the quasi-convex optimization problem
%discussed in problem 6.9, for the case that f(t) is a rational of
%quadratic functions of t in numerator and denominator

%Define tolerance
eps = 10^(-3);
k = 201;
%Define bounds
u = 1;
l = -1;

%Create vector of t values
teas = zeros(k,1);
for i = 1:k;
    teas(i) = -3 + 6*(i-1)/(k-1);
end


while (u-l) > eps;
    av = (l+u)/2;
    cvx_begin
        variable a(3);
        variable b(2);
        minimize(0);
        subject to
            max(abs(a(1)+a(2)*teas+a(3)*teas.^2 - exp(teas).*(1 + b(1)*teas + b(2)*teas.^2)) - av*(1+b(1)*teas + b(2)*teas.^2)) <= 0;
    cvx_end
    if (strcmp(cvx_status,'Solved')) 
        u = av;
        a_star = a;
        b_star = b;
    else
        l = av;
    end
end


%Plot the result of this quasi-convex optimization problem
f = zeros(k,1);
for i = 1:k;
    f(i) = (a_star(1) + a_star(2)*teas(i) + a_star(3)*teas(i)^2)/(1+b_star(1)*teas(i) + b_star(2)*teas(i)^2);
end

figure(1)
plot(teas, f, 'black x', teas, exp(teas), 'red');
xlabel('t');
ylabel('Value');
title('Plots of f(t) and y');
legend('f(t)', 'y');


%Also plot the fitting error
figure(2)
plot(teas, f - exp(teas), 'red');
xlabel('t');
ylabel('fitting error');
title('Plot of f - y');