%We proceed to solve a simple quadratic program in cvx, for a particular
%choice of parameters
%Part A)

n = 2;
Q = [1, -1; 0, 2];
v = [-1; 0];
w = 0;

cvx_begin
    variable x(n);
    dual variables y1 y2 y3;
    minimize(quad_form(x, Q, v, w));
    subject to
        y1 : x(1) + 2*x(2) <= -2;
        y2 : x(1) - 4*x(2) <= -3;
        y3 : 5*x(1) + 76*x(2) <= 1;
cvx_end

disp(y1)
disp(y2)
disp(y3)
disp(x)
