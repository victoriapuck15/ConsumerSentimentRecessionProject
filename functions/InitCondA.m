function [ A, C, Q, R, initZ, initV] = InitCond5B(X,K,p,Transformation)
 
[T,N] = size(X);

% Zero step - filter gaps using HP filter (only if indicated)
X_trends = repmat(NaN, size(X));
X_gaps = repmat(NaN, size(X));
[~, ind] = ismember(Transformation,'gap');
%X_trends(:,find(ind)) = hpf(X_st(:,find(ind)), 'lambda=', 1600);
X_trends = hpf(X, 'lambda=', 1600);
X_gaps(:,find(ind)) = X(:,find(ind)) - X_trends(:,find(ind));
X_gaps(:,find(abs(ind-repmat(1,size(ind))))) = X(:,find(abs(ind-repmat(1,size(ind)))));

% Standardise
[T,N] = size(X_gaps);
Mx = nanmean(X_gaps);
Wx = (nanstd(X_gaps));
X_gaps_st = (X_gaps-repmat(Mx,T,1))./repmat(Wx,T,1);

% ...
OPTS.disp=0;
[v,d] = eigs(cov(X_gaps_st),K,'lm',OPTS);
L = v;
f = X_gaps_st*v;
F0 = f(:,1:K);
R0 = eye(N);

% Run a VAR in F, obtain initial B and Q
[B,~,~,Q,~] = estvar(F0,p,[]);

% adjust for lags in state equation, Q is KxK
Q0 = [Q zeros(K,K*(p-1));zeros(K*(p-1),K*p)];
B = [B(:,:);eye(K*(p-1)) zeros(K*(p-1),K)];

% 1. Observation equation
% y = C*Z + e
I0 = eye(N,N);
I0(:,find(abs(ind-repmat(1,size(ind))))) = 0;
C = [eye(N,N) I0 zeros(N,K)];
R = zeros(N,N);

% 2. Transition equation
% Z = A*Z(-1) + v
A = [zeros(N,N) zeros(N,N) L;
     zeros(N,N) I0 zeros(N,K);
     zeros(K,N) zeros(K,N) B];

Q = [R0 zeros(N,N) zeros(N,K);
     zeros(N,N) I0 zeros(N,K);
     zeros(K,N) zeros(K,N) Q0];

% 3. Starting values for KF
initZ = zeros(K*p+N*2,1); % starting value for vector Z
%initV = eye(K*p);
initV = eye(K*p+N*2); % starting matrix P

