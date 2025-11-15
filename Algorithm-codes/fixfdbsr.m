function [Z,E] = fixfdbsr(X,lambda,canshu,he)


tol = 1e-4;
maxIter = 1e3;
[d ,n] = size(X);
rho = 1.1;

max_mu = 1e30;
mu = 1e-6;
xtx = X'*X;
inv_x = inv(xtx+2*eye(n));
%% Initializing optimization variables
% intialize
J = zeros(n,n);
Z = zeros(n,n);
E = sparse(d,n);

Y1 = zeros(d,n);
Y2 = zeros(n,n);
Y3 = zeros(n,n);
%% Start main loop
iter = 0;
DEBUG = 1;
if DEBUG
disp(['initial,rank=' num2str(rank(Z))]);
end
while iter<maxIter
    iter = iter + 1;
    
    temp = Z + Y2/mu;
    J = max(0,temp - 1/mu)+min(0,temp + 1/mu);   
    
    if trace(J'*J)==0
        J=diag(sqrt(he/n)*ones(n,1));
    else
        temptrace=sqrt(he/(trace(J'*J)));
        J=J*temptrace;
    end 
    
    temp = Z + Y3/mu;
    [UH,sigmaH,VH] = svd(temp,'econ');
    for time=1:n
        if sum(diag(sigmaH(1:time,1:time)))>canshu/mu
            if time<n
                if ((sum(diag(sigmaH(1:time,1:time)))-canshu/mu)/time)>=sigmaH(time+1,time+1)
                    tsigma=(sum(diag(sigmaH(1:time,1:time)))-canshu/mu)/time;
                    break
                end
            else
                tsigma=(sum(diag(sigmaH(1:time,1:time)))-canshu/mu)/time;
            end
        else
            if time==n
                tsigma=0;
            end
        end
    end
    
    for cba=1:time
        sigmaH(cba,cba)=tsigma;
    end
    H = UH*sigmaH*VH';
    

    Z = inv_x*(xtx-X'*E+J+H+(X'*Y1-Y2-Y3)/mu);
    
    xmaz = X-X*Z;
    temp = X-X*Z+Y1/mu;
    E = solve_l1l2(temp,lambda/mu);
    
    leq1 = xmaz-E;
    leq2 = Z-J;
    leq3 = Z-H;
    stopC = max([max(max(abs(leq1))),max(max(abs(leq2))),max(max(abs(leq3)))]);
    if DEBUG
        if iter==1 || mod(iter,50)==0 || stopC<tol
            disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e') ...
                ',rank=' num2str(rank(Z,1e-3*norm(Z,2))) ',stopALM=' num2str(stopC,'%2.3e')]);
        end
    end
    if stopC<tol
        break;
    else
        Y1 = Y1 + mu*leq1;
        Y2 = Y2 + mu*leq2;
        Y3 = Y3 + mu*leq3;
        mu = min(max_mu,mu*rho);
    end
end

function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end


function [x] = solve_l2(w,lambda)
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end