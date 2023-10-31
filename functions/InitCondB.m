function [ A, C, Q, R, initZ, initV] = InitCond4B(X,K,p)

X = X(:,1:end);

% Standardise
[T,N] = size(X);
Mx = nanmean(X);
Wx = (nanstd(X));
X_st = (X-repmat(Mx,T,1))./repmat(Wx,T,1);

[T,N]=size(X_st);

OPTS.disp = 0;
[v,d] = eigs(cov(X_st(:,1:end)),K,'lm',OPTS);
L = v;
f = X_st(:,1:end)*v;
F0 = f(:,1:K);
R = eye(N);

% Run a VAR in F, obtain initial B and Q
[B,Bc,v,Q,invFF]=estvar(F0,p,[]);

% adjust for lags in state equation, Q is KxK
Q = [Q zeros(K,K*(p-1));zeros(K*(p-1),K*p)];
B = [B(:,:);eye(K*(p-1)) zeros(K*(p-1),K)];

% start with
Sp = zeros(K*p,1);
Pp = eye(K*p);

% Final set of initial conditions
A = B;
C = L;
Q = Q;
R = R;
initZ = Sp;
initV = Pp;
