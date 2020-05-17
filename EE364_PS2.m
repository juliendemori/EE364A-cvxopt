%This is the code I used to minimize the various objectives subject to the
%same constraints. I simply modified the argument of minimize() for each
%case a-e.
cvx_begin
    variables x1 x2;
    minimize(x1 + x2);
    subject to
        1 <= 2*x1 + x2;
        1 <= x2 + 3*x2;
        0 <= x1;
        0 <= x2; 
cvx_end

disp(x1);
disp(x2);