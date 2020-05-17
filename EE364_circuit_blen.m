% data file for circuit blending problem

n = 10; % number of variables
k = 6;  % number of designs

% component widths from known designs
% each column of W is a different design
W = [1.8381    1.5803   12.4483    4.4542    6.5637    5.8225;
     1.0196    3.0467   18.4965    3.6186    7.6979    2.3292;
     1.6813    1.9083   17.3244    4.6770    4.6581   27.0291;
     1.3795    2.6250   14.6737    4.1361    7.1610    7.5759;
     1.8318    1.4526   17.2696    3.7408    2.2107   10.3642;
     1.5028    3.0937   14.9034    4.4055    7.8582   20.5204;
     1.7095    2.1351   10.1296    4.0931    2.9001    9.9634;
     1.4289    3.5800    9.3459    3.8898    2.7663   15.1383;
     1.3046    3.5610   10.1179    4.3891    7.1302    3.8139;
     1.1897    2.7807   13.0112    4.2426    6.1611   29.6734];

W_min = 1;
W_max = 30;

% objective values for the different designs
% entry j gives the objective for design j
P = [29.0148   46.3369  282.1749   78.5183  104.8087  253.5439];
D = [15.9522   11.5012    4.8148    8.5697    8.0870    6.0273];
A = [22.3796   38.7908  204.1574   62.5563   81.2272  200.5119];

% specifications
P_spec = 60;
D_spec = 10;
A_spec = 50;


%We proceed to blend the circuits
%i.e. to determine a vector of theta parameters that add to one, such that
%we obtain a theta in the feasible set

cvx_begin
    variable theta(k);
    subject to
        theta'*log(P') <= log(P_spec);
        theta'*log(D') <= log(D_spec);
        theta'*log(A') <= log(A_spec);
        theta'*ones(k,1) == 1;
        theta >= 0;
cvx_end

disp(theta)

%Given the obtained theta, we compute the optimal set of widths w_opt that is
%the convex combination of circuit designs that lead to the above
%objectives
lw_opt = 0;
for i = 1:k;
    lw_opt = lw_opt + theta(i)*log(W(:,i));
end

%Finally computing the optimal circuit design
w_opt = exp(lw_opt);