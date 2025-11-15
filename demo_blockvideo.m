clear 
close all
addpath Ncut_9
addpath common
addpath Algorithm-codes
addpath data

% H=10;
% acc=[];
% NMI=[];
% for i=1:H
%% real data human active 
n_space = 8;
load video8
A1 = fea;
corruption=0;
% 0 0.00001 0.001  0.1
N = randn(size(A1)) * corruption;
X1 = A1 + N;
% X = X1;
X = normalize(X1);
label = gnd; 
%% solve PSNR
s=max(max(A1));
[u,v]=size(A1);
C=A1-X1;
c=norm(C,'fro');
PSNR = 10*log10(s*s/(c^2/u/v));

%% running FeaSCLRR
[~,n]=size(X);
T1 = n;
s = 12;
[W] = denW1(X,T1,s);
rho =1.8;
lambda1 = 0.05;
lambda2 = 0;
lambda3 = 0.5;
% T1 = n;
T2 = 0.1;
K = 400;
% K = size(A1,1); 
tic; 
[Z9,p,W1] = FTSCLRR(X,lambda1,lambda2,lambda3,rho,K,W,T2);
t9=toc;
% plot(diag(W1)');
A9 = (abs(Z9)+abs(Z9'))/2;
clusters9 = ncutW(A9, n_space);
final_clusters9 = condense_clusters(clusters9, 1);
NcutDiscrete9 = clusters9;
accuracy9 = compute_accuracy(NcutDiscrete9, label);
n9 = nmi(label, final_clusters9);

figure,imagesc(A9),xlabel('feaSCLRR');

plotClusters(final_clusters9);xlabel('feaSCLRR');
