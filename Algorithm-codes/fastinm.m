function [Z,E] = fastinm(X,lambda,kkk)

if nargin<2
    lambda = 1;
end
tol = 1e-4;
maxIter = 1e3;
[d n] = size(X);

if nargin<3
    kkk = floor(n*n/2);
end

rho = 1.1;
max_mu = 1e30;
mu = 1e-6;


xtx = X'*X;

[tempXU,tempXS,tempXV]=svd(xtx);
weidu=length(find(diag(tempXS)>=1e-6));
XU=tempXU(:,1:weidu);
XS=tempXS(1:weidu,1:weidu);
newXS=(1./((XS<1e-6)+XS)).*(XS>=1e-6);
tempU=XU*inv(newXS+XU'*XU);



%% Initializing optimization variables
% intialize
J = zeros(n,n);
Z = zeros(n,n);
E = sparse(d,n);

Y1 = zeros(n,n);
Y2 = zeros(d,n);
%% Start main loop
iter = 0;
DEBUG = 1;
if DEBUG
disp(['initial,rank=' num2str(rank(Z))]);
end
while iter<maxIter
    iter = iter + 1;
        
    tempj = Z + Y1/mu;
    abstempj=reshape(abs(tempj),[1,n*n]);
    [allvalue,allindex]=sort(abstempj,'descend');

    value=allvalue(1:kkk);
    index=allindex(1:kkk);
    cumww=cumsum(ones(1,kkk));
    cumori=cumsum(value)-1/mu;
    valtemp(1:kkk-1)=value(2:kkk);
    valtemp(kkk)=0;
    cumtemp=(cumori./cumww)>valtemp;
    
    J=tempj;
    
    [mm,nn]=find(cumtemp==1);
    if numel(nn)>0
        tkw=cumori(nn(1))/cumww(nn(1));
        J(index(1:nn(1)))=(sign(J(index(1:nn(1))))*tkw);
        J(allindex(kkk+1:end))=0;
    else
        J=zeros(n);
    end
    
    tempZ=xtx-X'*E+J+(X'*Y2-Y1)/mu;
    Z=tempZ-tempU*(XU'*tempZ);
    
    xmaz = X-X*Z;
    tempe = xmaz+Y2/mu;
    
    E = solve_l1l2(tempe,lambda/mu);
    
    leq1 = Z-J;
    leq2 = xmaz-E;
    
    stopC = max([max(max(abs(leq1))),max(max(abs(leq2)))]);
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