%Spacecraft landing problem

h = 1;
g = 0.1;
m = 10;
Fmax = 10;
p0 = [50,50,100]';
v0 = [-10,0,-10]';
alpha = 0.5;
gamma = 1;
K = 35;

%Making necessary matrices
I = eye(3,3);
e3 = I(3,:);



%PART a)
cvx_begin
    variables f(3,K) p(3,K) v(3,K);
    minimize(sum(norm(f, 'fro')));
    subject to
        for i = 1:K;
           p(3,i) >= alpha *norm([p(1,i), p(2,i)]); 
        end
       (sum(f.^2)) <= Fmax^2;
        p(:,1) == p0;
        v(:,1) == v0;
        p(:,K) == 0;
        v(:,K) == 0;
        p(:, 2:K) == p(:,1:K-1) + (1/2)*(v(:,1:K-1) + v(:,2:K));
        v(:,2:K) == v(:,1:K-1) + (1/m)*f(:,1:K-1) - g*e3'*ones(1,K-1);
cvx_end

% 
% T = 0;
% %PART b)
% for k = 1:K
%    cvx_begin
%     variables f(3,K) p(3,K) v(3,K);
%     minimize(0);
%     subject to
%         for i = 1:K;
%            p(3,i) >= alpha *norm([p(1,i), p(2,i)]); 
%         end
%         (sum(f.^2)) <= Fmax^2;
%         p(:,1) == p0;
%         v(:,1) == v0;
%         p(:,K) == 0;
%         v(:,K) == 0;
%         p(:, 2:K) == p(:,1:K-1) + (1/2)*(v(:,1:K-1) + v(:,2:K));
%         v(:,2:K) == v(:,1:K-1) + (1/m)*f(:,1:K-1) - g*e3'*ones(1,K-1);
%    cvx_end
%    if cvx_optval == 0;
%        T = k;
%        break;
%    end
% end

% use the following code to plot your trajectories
% plot the glide cone (don't modify)
% -------------------------------------------------------
x = linspace(-40,55,30); y = linspace(0,55,30);
[X,Y] = meshgrid(x,y);
Z = alpha*sqrt(X.^2+Y.^2);
figure; colormap autumn; surf(X,Y,Z);
axis([-40,55,0,55,0,105]);
grid on; hold on;

% INSERT YOUR VARIABLES HERE:
% -------------------------------------------------------
plot3(p(1,:),p(2,:),p(3,:),'b','linewidth',1.5);
quiver3(p(1,1:K),p(2,1:K),p(3,1:K),...
        f(1,:),f(2,:),f(3,:),0.3,'k','linewidth',1.5);
