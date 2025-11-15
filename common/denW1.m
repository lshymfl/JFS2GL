function [W] = denW1(X,T1,s)
[~,n] = size(X);
for i = 1:n
    for j = 1:n
        if abs(i-j)<=s
            W(i,j) = 0;
        else
       W(i,j) = 1 - exp(-abs(i-j)/T1);
        end
    end
end