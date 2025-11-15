function A = construcA(Z,X)
[~,n] = size(Z);
for i=1 : n
    for j = 1:n
        G1(i,j) = Z(:,i)'*Z(:,j);
        G2(i,j) = norm(X(:,i))*norm(X(:,j));
        A(i,j) = abs(G1(i,j)/G2(i,j));
    end
end