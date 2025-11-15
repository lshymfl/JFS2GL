function [Z,p] = LapFLRR_v2(X,alpha,beta,options,M,rho,kappa,maxiter,toler)
%--------------------------------------------------------------------------
% This is the function to call the Feature Selection Embedded Subspace 
% Clustering with Nuclear Norm. The objective function is as follows:
%     min_{Z,p} 1/2||diag(p)X - diag(p)XZ||_F^2 + alpha ||Z||_*
% X is the data set of size d x n. 
% alpha is the regularization parameter.
% M controls the number of features to select.
%--------------------------------------------------------------------------
% For use of this code, please cite the following reference:
% 
% C. Peng; Z. Kang; M. Yang; Q. Cheng, "Feature Selection Embedded 
% Subspace Clustering," in IEEE Signal Processing Letters , vol.PP, 
% no.99, pp.1-1 doi: 10.1109/LSP.2016.2573159

% Copyright @ Chong Peng, 2016
%--------------------------------------------------------------------------

[d,n] = size(X);
I = eye(n);
Z0 = I; W0 = I; U0 = I; p0 = ones(d,1); Theta_W = zeros(n);
% L = eye(n);
for t = 1:maxiter
 
     A = constructW((diag(p0)*X),options);
    DCol = full(sum(A,2));
    B = spdiags(DCol,0,n,n);
    L = A-B;
    
    D = diag(p0)*X;
    DD = D'*D;
    DD = full(DD);
    L = full(L);
    W = lyap( 2*DD+rho*I,beta*(L+L'),-(2*DD+rho*Z0-Theta_W) );
%     W = sylvester(2*DD+rho*I,beta*(L+L'),2*DD+rho*Z0-Theta_W);

    Z = cal_thresh_nuclear( W+Theta_W/rho,alpha/rho );
    
    p = solve_p( X-X*W,M,d);
    
     
    
    Theta_W = Theta_W + rho*(W-Z);
    rho = min(rho*kappa,1e8);

    err = max([norm(Z-Z0,'fro'), norm(W-W0,'fro'), norm(Z-W,'fro'), norm(p-p0,'fro')]);
    
    if err <= toler
        break
    end
    
    W0 = W; Z0 = Z; p0 = p; 
    
end

end

function p = solve_p(G,M,d)
p = zeros(d,1);
G = sum(G.^2,2);

[a,~] = sort(G);
th = a(M);

p(G<=th) = 1;

p = sparse(p);

end

function Y = cal_thresh_nuclear(M,lambda)
[U,S,V] = svd(M);
S = diag(S);
S = S - lambda;
S = max(S,0);
Y = U*diag(S)*V';

end

function [E] = solve_L21_column(W, lambda)
n = size(W,2);
E = W;
for i = 1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
end

function [E] = solve_L21_row(W, lambda)
n = size(W,1);
E = W;
for i = 1:n
    E(i,:) = solve_l2(W(i,:),lambda);
end
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end
end
