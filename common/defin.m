function W = defin(X)
[~,n] = size(X);
for i = 1:n
    for j = 1:n
        NX(i,j) = norm(X(:,i)-X(:,j));
    end
end
sigma = mean(mean(NX));
K = exp(-NX/sigma);
W = 1 - K;
