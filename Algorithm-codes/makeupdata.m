function D = makeupdata (X)
[~,n] = size(X);
for i=1 : n-1
    D(i) = norm(X(:,i+1)-X(:,i))^2;
end