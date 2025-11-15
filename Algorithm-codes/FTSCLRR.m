function [Z,p,W1] = FTSCLRR(X,lambda1,lambda2,lambda3,rho,K,W,T2)
% min_{p,Z} 1/2||diag(p)X-diag(p)XZ||_F^2+lambda1||Z||_*+lambda2||W.*Z||_1 + lambda3||ZRW1||_{2,1}
%                 s.t.   ||p||_0 = M,Z>=0

tol = 1e-5;
[m, n] = size(X);
% maxIter = 500;
max_mu = 1e8;
mu = 0.1;
I = eye(n);
R = (triu(ones(n,n-1),1) - triu(ones(n, n-1))) + (triu(ones(n, n-1),-1)-triu(ones(n, n-1)));
%% Initializing optimization variables
p = ones(m,1);
Z = zeros(n);
% L = zeros(n);
% N = zeros(n,n-1);
Y1 = zeros(n);
Y2 = zeros(n);
Y3 = zeros(n,n-1);
%% Start main loop
max_iterations = 400;
% func_vals = [];
for k = 1 : max_iterations
%% update Q
Xp = (diag(p)*X)'*diag(p)*X;
Q1 = Xp + mu*I;
Q2 = Xp + mu*Z + Y1;
Q = Q1\Q2;

%% update weighted W1 
X1 = diag(p)*X;
w = makeupW1(X1,T2);
W1 = diag(w);
% W1 = eye(n-1);
%% update L
tempe = Z + Y2/mu;
tau1 = lambda1/mu;
L = solve_nn( tempe, tau1);
%% update N 
B = Z*R*W1 + Y3/mu;
gamma1 = lambda3/mu;
N = solve_l21(B,gamma1);
%% update Z
R1 = R*W1;
Rw = mu*Z*(R1*R1') - mu*N*R1' + Y3*R1';
Deth = 2*mu*Z + Rw - mu*(Q+L) + Y1 + Y2;
eta = 2*mu + mu*norm(R1)^2 +1; sigma = mu*eta;   
temp = Z - Deth/sigma;
tau = lambda2/sigma;
Z = max(0,temp - W*tau)+min(0,temp + W*tau);
Z = max(Z,0);
%% update p
% H = X-X*Q;
% for j =1:m
% q1(j) = norm(H(j,:))^2;
% end
% q = sqrt(q1);
% eta1 = 1./q;
% eta = 1/sum(eta1);
% p = max(eta./q,0);
% % p = p';
% [idx,~] = sort(p,'descend');
% for i = 1:m
%  if p(i) >= idx(K)
%      p(i) = p(i);
%  else
%      p(i) = 0; 
%  end
% end

H = X-X*Q;
for j =1:m
h1(j) = norm(H(j,:));
end
[idx,~] = sort(h1);
for i = 1:m
 if h1(i) <= idx(K)
     p(i) = 1;
 else
     p(i) = 0; 
 end
end
%% update multiplier
dY1 = Z - Q; dY2 = Z - L; dY3 = Z*R*W1 - N; 
Y1 = Y1 + mu*dY1;
Y2 = Y2 + mu*dY2;
Y3 = Y3 + mu*dY3;
mu = min(max_mu, mu*rho);

[~,Sigma,~] = svd(Z);
% func_vals(k) = 0.5 * norm(diag(p)*X - diag(p)*X*Z,'fro')^2 + lambda1*sum(sum(Sigma)) + lambda2*sum(sum(W.*Z)) + lambda3*norm_l1l2(Z*R*W1);
%% Check convergence
     stopC = max([max(max(abs(dY1))),max(max(abs(dY2))),max(max(abs(dY3)))]);
%      stopC = norm(Z-Zk,'fro')/norm(Zk,'fro');
    convergenced = stopC;
 if convergenced<tol
        break;
 end
end

% fprintf(' func_vals\n');
% fprintf(' %7.4f \n',func_vals);
% end
