function A = construcAff(Z,fai)
[m,n] = size(Z);
U = Z*Z';
for i=1 : m
    for j = 1:n
         A(i,j) = abs(U(i,j))^fai;
    end
end