%We proceed to smoothly interpolate the unknown parts of an image by
%minimizing first the L1 variation squared, and then the total variation


%Some data first
% tv_img_interp.m
% Total variation image interpolation.
% EE364a
% Defines m, n, Uorig, Known.

% Load original image.
Uorig = double(imread('tv_img_interp.png'));

[m, n] = size(Uorig);

% Create 50% mask of known pixels.
rand('state', 1029);
Known = rand(m,n) > 0.5;


% Calculate and define Ul2 and Utv.
%L1 Variation

[m,n] = size(Uorig);
s1 = (m-1)*n;
s2 = (n-1)*m;


%Constrain every value in the matrices of differences to be the difference
%between adjacent (either vertical or horizontal) values.
%In addition, constrain each value in U_l to be the value of Uorig, if it
%is known.
cvx_begin
    variables U_l(m,n) U1(s1) U2(s2);
    minimize(sum_square(U1) + sum_square(U2));
    subject to   
        for i = 1:m;
            for j = 1:n;
                if Known(i,j) == 1;
                    U_l(i,j) == Uorig(i,j);
                end
            end
        end
        for i = 2:m;
            for j = 1:n;
                U1((i-2)*(n-1) + j) == U_l(i, j) - U_l(i-1, j);
            end
        end

        for j = 2:n;
            for i = 1:m;
                U2((j-2)*(m-1) + i) == U_l(i,j) - U_l(i, j-1);
            end
        end
cvx_end

%total variation
cvx_begin
    variables U_t(m,n) U1(s1) U2(s2);
    minimize(sum(abs(U1)) + sum(abs(U2)));
    subject to
        for i = 1:m;
            for j = 1:n;
                if Known(i,j) == 1;
                    U_t(i,j) == Uorig(i,j);
                end
            end
        end
        for i = 2:m;
            for j = 1:n;
                U1((i-2)*(n-1) + j) == U_t(i, j) - U_t(i-1, j);
            end
        end

        for j = 2:n;
            for i = 1:m;
                U2((j-2)*(m-1) + i) == U_t(i,j) - U_t(i, j-1);
            end
        end
cvx_end
    

%%%%%

% Graph everything.
figure(1); cla;
colormap gray;

subplot(221);
imagesc(Uorig)
title('Original image');
axis image;

subplot(222);
imagesc(Known.*Uorig + 256-150*Known);
title('Obscured image');
axis image;

subplot(223);
imagesc(U_l);
title('l_2 reconstructed image');
axis image;

subplot(224);
imagesc(U_t);
title('Total variation reconstructed image');
axis image;