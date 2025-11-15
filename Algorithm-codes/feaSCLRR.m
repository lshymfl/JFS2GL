function [Z,p,W,W1] = feaSCLRR(X,lambda1,lambda2,lambda3,rho,K,T1,T2)
% min_{p,Z} 1/2||diag(p)X-diag(p)XZ||_F^2+lambda1||Z||_*+lambda2||W.*Z||_1 + lambda3||ZRW1||_{2,1}
%                 s.t.   ||p||_0 = M,  z_i'*1=1, z_i >= 0

tol =  1e-4;
[m, n] = size(X);
maxIter = 500;
max_mu = 1e8;
mu = 1e-4;
I = eye(n);
R = (triu(ones(n,n-1),1) - triu(ones(n, n-1))) + (triu(ones(n, n-1),-1)-triu(ones(n, n-1)));
%% Initializing optimization variables
p = ones(m,1);
Z = eye(n);
% Q = zeros(n);
% M = zeros(n); 
% N = zeros(n,n-1);
Y1 = zeros(n);Y3 = zeros(n,n-1);Y2 = zeros(n);

s=8;
[W] = denoteW(X,s,T1);
%% Start main loop
iter = 0;

while iter<maxIter
    iter = iter + 1;
%% update Q
Xp = X'*diag(p)*X;
Q1 = Xp + mu*I;
Q2 = Xp + mu*Z - Y1;
Q = Q1\Q2;

%% update M
temp = Z - Y2/mu;
tau1 = lambda1/mu;
    [U,sigma,V] = svd(temp,'econ');
    Sigma = diag(sigma);
    svp = length(find(Sigma>tau1));
    if svp>=1
        Sigma = Sigma(1:svp)-tau1;
    else
        svp = 1;
        Sigma = 0;
    end
    M = U(:,1:svp)*diag(Sigma)*V(:,1:svp)';
 %% update p
H = X-X*Q;
for j =1:m
h1(j) = norm(H(j,:));
end
[idx,~] = sort(h1);
for i = 1:m
 if norm(H(i,:)) <= idx(K)
     p(i) = 1;
 else
     p(i) = 0; 
 end
end
 W1 = makeupW1(X,p,T2);
 
%% update N 
B = Z*R*W1 - Y3/mu;
gamma1 = lambda3/mu;
N = solve_l21(B,gamma1);
%% update Z
% Zk = Z;
R1 = R*W1;
Rw = (mu*Z*R1-mu*N-Y3)*R1';
Deth = 2*mu*Z + Rw - mu*(Q+M) - Y1 - Y2;
eta = 2*mu + mu*norm(R1)^2 +1;
sigma = mu*eta;   
C = Z - Deth/sigma;
gamma = lambda2/sigma;
for k = 1:n
V(:,k) = C(:,k) - gamma*W(:,k);
Z(:,k) = EProjSimplex(V(:,k));
end

%% update multiplier
dY1 = Q - Z; dY2 = M - Z; dY3 = N - Z*R*W1; 
Y1 = Y1 + mu*dY1;
Y2 = Y2 + mu*dY2;
Y3 = Y3 + mu*dY3;
mu = min(max_mu, mu*rho);

%% Check convergence
     stopC = max([max(max(abs(dY1))),max(max(abs(dY2))),max(max(abs(dY3)))]);
%      stopC = norm(Z-Zk,'fro')/norm(Zk,'fro');
    convergenced = stopC;
 if convergenced<tol
        break;
 end
end

function [W] = denoteW(X,s,T1)
[~,n] = size(X);
for i = 1:n
    for j = 1:n
       if abs(i-j) <= s/2
       W(i,j) = 0;
      else
       W(i,j) = 1 - exp(-abs(i-j)/T1);
      end
    end
end

function [W1] = makeupW1(X,p,T2)
[~,n] = size(X);
for k = 1:n-1
    Xx(:,k) = diag(p)*(X(:,k+1)-X(:,k));
end
for j=1 : n-1
    w(j) = exp(-norm(Xx(:,j))^2/T2);
end
W1 = diag(w);