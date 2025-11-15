function [W] = denoteW(Z)
 [m,n] = size(Z);
 for i = 1:m
     for j = 1:n
         W(i,j) = exp(-Z(i,j));
     end
 end


