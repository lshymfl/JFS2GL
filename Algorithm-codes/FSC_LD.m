function [Z,p] = FSC_LD(X,alpha,M,rho,kappa,maxiter,toler)
%--------------------------------------------------------------------------
% This is the function to call the Feature Selection Embedded Subspace 
% Clustering with Lod-Determinant rank approximation. The objective 
% function is as follows:
%     min_{Z,p} 1/2||diag(p)X - diag(p)XZ||_F^2 + alpha ||Z||_{ld},
% where ||Z||_{ld} = \sum_i log(1+\sigma_i(Z)) = logdet(I + Z'Z)
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
Z0 = I; J0 = I; p0 = ones(d,1); Theta = zeros(n); 

for t = 1:maxiter
    
    D = diag(p0)*X;
    DD = D'*D;
    J = (DD+rho*I)\(DD+rho*Z0-Theta);
    
    Z = XLogDet_solver(J+Theta/rho,rho/alpha/2);
	% Z = cal_thresh_logdet(J+Theta/rho,alpha/rho);
    
    p = solve_p( X - X*J,M,d);
    
    Theta = Theta + rho*(J-Z);
    rho = min(rho*kappa,1e8);

    err = max([norm(J-J0,'fro'), norm(Z-Z0,'fro'), norm(p-p0,'fro')]);
    
    if err <= toler
        break
    end
    
    J0 = J; Z0 = Z; p0 = p;
    
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

function [ X ] = XLogDet_solver(D,rho)

[U,S,V] = svd(D);
S0 = diag(S);
r = length(S0);

P = [ones(r,1), 1-S0, 1/2/rho-S0];

rt = zeros(r,1);

for t = 1:r
    
    p = P(t,:);
    Delta = p(2)^2-4*p(1)*p(3);
    if Delta <= 0
        rt(t) = 0;
    else
        
        rts = roots(p);
        rts = sort(rts);
        if rts(1)*rts(2)<=0
            rt(t) = rts(2);
        elseif rts(2)<0
            rt(t) = 0;
        else
            funval = log(1+rts(2))+rho*(rts(2)-S0(t)).^2;
            if funval>log(1+0)+rho*(0-S0(t)).^2;
                rt(t) = 0;
				else
				rt(t) = rts(t);
            end
        end
    end
    
end

sig = diag(rt);

X = U*sig*V';

end


function Y = cal_thresh_logdet(M,lambda)
[U,S,V] = svd(M,'econ');
S = diag(S);

ks = (S-1)/2+sqrt((1+S).^2/4-lambda);
fval = funval(ks,S,lambda);
f0 = (S.*S)/2/lambda;

th = (1+S).*(1+S)/lambda;

c = (( (fval <= f0).*(th>4) ) == 1);

S(:) = 0;
S(c) = ks(c);

Y = U*sparse(diag(S))*V';

end



