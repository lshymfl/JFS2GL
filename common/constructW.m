function G = constructW(T,options)
[~,n] = size(T);
for i=1 : n
    for j = 1:n
        G1(i,j) = T(:,i)'*T(:,j);
        G2(i,j) = norm(T(:,i))*norm(T(:,j));
        G(i,j) = G1(i,j)/G2(i,j);
    end
end