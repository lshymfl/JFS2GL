function G = construcG(Z)
[~,n] = size(Z);
for i=1 : n
    for j = 1:n
        G1(i,j) = Z(:,i)'*Z(:,j);
        G2(i,j) = norm(Z(:,i))*norm(Z(:,j));
        G(i,j) = G1(i,j)/G2(i,j);
    end
end