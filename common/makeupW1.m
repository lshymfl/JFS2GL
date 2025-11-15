function [w] = makeupW1(X1,T2)
[~,n] = size(X1);
for k = 1:n-1
    Xx(:,k) = X1(:,k+1)-X1(:,k);
end
for j=1 : n-1
    w(j) = exp(-(norm(Xx(:,j))^2)/T2);
end