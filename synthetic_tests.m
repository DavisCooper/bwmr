%% Define constants
addpath('funs')
d = 10;
r = 5;
n = 50;
eps = 10; %noise level
iter = 1e5;

%% Generate Problem
U_sol = randn(d,r);
X = randn(n,d);
Y = vecnorm(X*U_sol,2,2).^2;
Y_hat = Y + eps*(rand(n,1) - 0.5);

[Z,D,~] = svd((Y_hat.*X)'*X);
s = diag(D);
% U0FR = Z*sqrt(D);
% U0 = Z(:,1:r)*sqrt(diag(s(1:r)- s(r+1)));
U0FR = Z;
U0 = Z(:,1:r);

SigmaZ = eye(d);
%% Create polytope constraints
Y_u = Y_hat + eps/2;
Y_l = max(Y_hat - eps/2,0.1);
Y_ = [Y_u;Y_l];
sgn = [ones(n,1);-ones(n,1)];
X_ = [X;X];

%% Run methods
tic
[B_bwgd,~,~] = bwgd(X, Y_hat, r, iter, U0, SigmaZ);
toc

[B_bwgdfr,~,~] = bwgd(X, Y_hat, d, iter, U0FR, SigmaZ);

% [B_bwsgd,~,~] = bwsgd(X, Y, r, 1, U0, U_sol);
tic
[U_proj,~,U_i,l_i] = bw_proj(X_,Y_,U0,sgn,iter,Z(:,1:r));
B_proj = U_proj*U_proj'*U0*U0'*U_proj*U_proj';
toc

% [T_proj,~,U_i,l_i] = bw_proj(X_,Y_,U0*U0',sgn,iter,U0*U0');
% B_Tproj = T_proj*U0*U0'*T_proj;

[B_egdsqrt,~] = egd_sqrt_fr(X, Y_hat, 1, iter, U0FR, SigmaZ);

tic
[B_egd,~, ~] = egd(X, Y_hat.^2, r, 0.01, iter, U0, SigmaZ, 1);
toc

[B_egdfr,~, ~] = egd(X, Y_hat.^2, d, 0.01, iter, U0FR, SigmaZ, 1);

%% Compare Convergence





