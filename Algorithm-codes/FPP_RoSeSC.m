function [ Y ] =FPP_RoSeSC(X,H,lambda1,lambda2,rho)
%
% argmin 0.5||X'-WX'||_{F}^2 +\lambda_1||W||_{1} + \lambda_2||R'W||_{2,1}
% Outputs:
% X: observed data X. 
% Inputs:
%% Initializing
[~,n] = size(X);
U=zeros(n-1,n);  % rand(n-1,n)
% Y=normalize(rand(n));
 Y = zeros(n);  
% Y=0.01*rand(n,n);      %0.001*rand(n)
R = (triu(ones(n,n-1),1) - triu(ones(n, n-1))) + (triu(ones(n, n-1),-1)-triu(ones(n, n-1)));
mu = 2;
mu_max = 1e6;
% % rho_0 = 1.1;
B = H'*R';
C = X'*X;
M = norm(B)^2;
N = norm(X')^2;
tol_2 = 1e-4;
max_iterations = 200;
%  func_vals = [];
% isConverged=false;
for k = 1 : max_iterations

    U_prev = U;
    Y_prev = Y;
   gamma=1.1*(M/mu+N);  
    %% Solve for U
 tempU = B*Y_prev/mu+U_prev;
     U = tempU - solve_l2l1(tempU*mu,mu*lambda2)/mu;
    
    %% Solve for Y
    Y_hat = C-Y_prev*C-B'*(2*U - U_prev);
    tempY = Y_prev+Y_hat/gamma;
     Y = solve_l1( tempY, lambda1/gamma );
%         Y_hat = max(tempY - lambda1/gamma, 0);
%     Y = Y_hat + min(tempY + lambda1/gamma, 0);
    
    % Update mu
    
 
    
%     if ( norm(Y - Y_prev,'fro')/norm(Y ,'fro') < tol_2)
%         rho = rho_0;
%     else
%         rho = 1;
%     end
    mu_bar = rho * mu;
    mu = min(mu_max, mu_bar);
%    func_vals(k) = 0.5 * norm(X' - Y*X','fro')^2 + lambda1*norm_l1(Y) +lambda2*norm_l1l2(H'*R'*Y);  
    %% Check convergence
 
    
    stopCriterion=norm(Y-Y_prev,'fro')/norm(Y,'fro');
     if stopCriterion<tol_2
        break;
    end
    
end
%  fprintf(' func_vals\n');
% fprintf(' %7.4f \n',func_vals);
% end