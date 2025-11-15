function [W] = MACW(X1,sig)
[~,n] = size(X1);
t = sig^2;
for i = 1:n
    for j = 1:n
        if  i == j
            W(i,j) = 0;
        else
        W(i, j) = exp(-(norm(X1(:, i) - X1(:, j))^2)/t); 
        end
    end
end