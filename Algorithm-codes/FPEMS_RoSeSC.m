function [ Y ] =FPEMS_RoSeSC(X,H,lambda1,lambda2,mu)
%
% argmin 0.5||X'-WX'||_{F}^2 +\lambda_1||W||_{1} + \lambda_2||R'W||_{2,1}
% Outputs:
% X: observed data X. 
% Inputs:

%% Initializing
[~,n] = size(X);
U=randn(n-1,n);
% U_prev=U;
Y=randn(n);
% Y_prev=Y;
R = (triu(ones(n,n-1),1) - triu(ones(n, n-1))) + (triu(ones(n, n-1),-1)-triu(ones(n, n-1)));

% norm_two = lansvd(X, 1, 'L');
% mu = 1.1;
% mu_bar = 1e6;
B = H'*R';

% mu = 0.1;
% mu_max = 1;
% 
% rho = (norm(X,2)^2) * 1.2;
% 
% gamma_0 = 1.2;

tol_2 = 1e-4;
isConverged=false;
%% Start Looping
while ~isConverged
    U_prev=U;
    Y_prev=Y;

% mu = norm(B)^2/(gamma-norm(X')^2);
gamma=1.1*(norm(B)^2/mu+norm(X')^2); 
 %% update U
%     tempU = B*Y_prev/mu+U_prev;
%      U1 = proxl21(tempU*mu,mu*lambda2);
%      U = tempU - U1/mu;
     tempU = B*Y_prev/mu+U_prev;
     U1 = proxl21(tempU*mu,mu*lambda2);
     U = tempU - U1/mu;
 %% update W   
    U2 = 2*U - U_prev;
    Y_hat = (X'-Y_prev*X')*X-B'*U2;
    tempY = Y_prev+Y_hat/gamma;
     Y = solve_l1( tempY, lambda1/gamma );
%      Y_bar=max(tempY - lambda1/gamma,0);
%      Y = Y_bar + min(tempY + lambda1/gamma,0);
     
  %% update mu
%  if (mu * max(sqrt(rho) * norm(U - U_prev,'fro'), norm(Y - Y_prev))/norm(X,'fro') < tol_2)
%         gamma1 = gamma_0;
%     else
%         gamma1 = 1;
%  end
%     
%     mu = min(mu_max, gamma1 * mu);
%     gamma=1.1*(norm(B)^2/mu+norm(X')^2); 
%% check strop criterion
    stopCriterion=norm(Y-Y_prev,'fro')/norm(Y,'fro')
     if stopCriterion<tol_2
         isConverged=true;
    end
    
end



