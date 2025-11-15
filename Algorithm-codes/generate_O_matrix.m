function O = generate_O_matrix(nSmp)
%Generate the smooth regularization matrix.
% 
% tm = -2*eye(nSmp);
% te = eye(nSmp-1);
% tem = [zeros(1,nSmp);[te, zeros(nSmp-1,1)]];
% 
% O = tm+ tem +tem';
% O(1) = -1;
% O(end) = -1;

tmp=diag(ones(nSmp-1,1),1);
O=tmp+tmp'-2*eye(nSmp);
O(1) = -1;
O(end) = -1;

end