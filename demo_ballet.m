clear 
close all
addpath Ncut_9
addpath common
addpath Algorithm-codes
addpath data

% H=8;
% acc=[];
% NMI=[];
% for i=1:H
%% real data human active 
n_space = 8;
load ba8         % ba3  ba6  ba8
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
%% running SSC
% tic;
% Z1 = ssc_exact_fro(X, 0.1);
% t1=toc;
% clusters1 = ncutW((abs(Z1) + abs(Z1'))/2, n_space);
% final_clusters1 = condense_clusters(clusters1, 1);
% NcutDiscrete1 = clusters1;
% accuracy1 = compute_accuracy(NcutDiscrete1, label);
% n1 = nmi(label, final_clusters1);
% %% running LRR
% tic;
% Z2 = lrr_relaxed(X,0.01);
% t2=toc;
% clusters2 = ncutW((abs(Z2) + abs(Z2'))/2, n_space);
% final_clusters2 = condense_clusters(clusters2, 1);
% NcutDiscrete2 = clusters2;
% accuracy2 = compute_accuracy(NcutDiscrete2, label);
% n2 = nmi(label, final_clusters2);
% %% running v-LADMAP-OSC
% tic;
% lambda_1 = 0.099;
% lambda_2 = 0.001;
% Z3 = osc_relaxed(X, lambda_1, lambda_2);
% t3=toc;
% clusters3 = ncutW((abs(Z3) + abs(Z3'))/2, n_space);
% final_clusters3 = condense_clusters(clusters3, 1);
% NcutDiscrete3 = clusters3;
% accuracy3 = compute_accuracy(NcutDiscrete3, label);
% n3 = nmi(label, final_clusters3);
% %% running ADMM-TSC 
% paras = [];
% paras.lambda1 = 0.01;
% paras.lambda2 = 15;
% paras.n_d = 60;
% paras.ksize = 7;
% paras.tol = 1e-4;
% paras.maxIter = 120;
% paras.stepsize = 0.1;
% tic;
% %%%---Learn representations Z---%%%
% [D, Z4, err] = TSC_ADMM(X,paras);
% t4 = toc;
% % disp('Segmentation...');
% %%%---Graph construction and segmentation---%%%
% 
% nbCluster = length(unique(label));
% vecNorm = sum(Z4.^2);
% W2 = (Z4'*Z4) ./ (vecNorm'*vecNorm + 1e-6);
% 
% [oscclusters,~,~] = ncutW(W2,nbCluster);
% final_clusters4 = denseSeg(oscclusters, 1);
% 
% accuracy4 = compacc(final_clusters4, label);
% n4 = nmi(label, final_clusters4);
% %% running BD_QOSC   0.01  0.61
% k = n_space;
% p = 1.1;
% lambda_1 = 0.5;
% lambda_2 =0.1;
% gamma_1 = 0.01;
% maxIterations = 500;
% tic;
% Z5 = BD_QOSC( X,k,lambda_1, lambda_2, gamma_1,p, maxIterations);
% t5=toc;
% clusters5 = ncutW((abs(Z5)+abs(Z5'))/2 , n_space);
% final_clusters5 = condense_clusters(clusters5, 1);
% NcutDiscrete5 = clusters5;
% accuracy5 = compute_accuracy(NcutDiscrete5, label);
% n5 = nmi(label, final_clusters5);
% %% running ADM-SCLRR  
% tic;
% NX=X;
% for r=1:size(NX,2)
%         NX(:,r)=NX(:,r)/norm(NX(:,r));
% end
%  MX = abs(NX'*NX);
%  thta = mean(mean(1-MX));
%  K = exp(-(1-MX)/thta);
%  W=1-K;
%  lambda=0.6;
%  beta=1.5;
%  [Z6,E] = admsclrr(X,lambda,W,beta);
% t6=toc;
% clusters6 = ncutW((abs(Z6)+abs(Z6'))/2 , n_space);
% final_clusters6 = condense_clusters(clusters6, 1);
% NcutDiscrete6 = clusters6;
% accuracy6 = compute_accuracy(NcutDiscrete6, label);
% n6 = nmi(label, final_clusters6);
% %% running FPEMS-WOSC  
% rho = 0.5;
% t = 10;
% [~,n] = size(A1);
% V = A1 + N;
% H1=makeupW(V,t);
% H = diag(H1);
% R = (triu(ones(n,n-1),1) - triu(ones(n, n-1))) + (triu(ones(n, n-1),-1)-triu(ones(n, n-1)));
% lambda1 = 0.5;
% lambda2 = 0.5;
% tic;
% Y = FPP_RoSeSC(X,H,lambda1,lambda2,rho);
% Z7=Y';
% t7=toc;
% clusters7 = ncutW(abs(Z7) + abs(Z7'), n_space);
% final_clusters7 = condense_clusters(clusters7, 1);
% NcutDiscrete7 = clusters7;
% accuracy7 = compute_accuracy(NcutDiscrete7, label);
% n7 = nmi(label, final_clusters7);
%% running FSCNN
% rho = 1.2;  kappa = 1.5;    alpha = 0.4;   
% M = 2600;    toler = 1e-3;   maxiter = 500;
% tic; 
% [Z8,~] = FSC_NN(X,alpha,M,rho,kappa,maxiter,toler);
% t8=toc;
% % plot(diag(W1)');
% A8 = (abs(Z8)+abs(Z8'))/2;
% clusters8 = ncutW(A8, n_space);
% final_clusters8 = condense_clusters(clusters8, 1);
% NcutDiscrete8 = clusters8;
% accuracy8 = compute_accuracy(NcutDiscrete8, label);
% n8 = nmi(label, final_clusters8);
% %% running FeaMAC
% rho = 0.5;  kappa = 1.1;    alpha = 0.3;   beta = 10;
% M = 20;    toler = 1e-3;   maxiter = 500;
% options = 0; 
% tic; 
% [Z81,p] = LapFLRR_v2(X,alpha,beta,options,M,rho,kappa,maxiter,toler);
% t81=toc;
% % plot(diag(W1)');
% A81 = (abs(Z81)+abs(Z81'))/2;
% clusters81 = ncutW(A81, n_space);
% final_clusters81 = condense_clusters(clusters81, 1);
% NcutDiscrete81 = clusters81;
% accuracy81 = compute_accuracy(NcutDiscrete81, label);
% n81 = nmi(label, final_clusters81);
%% running FeaSCLRR
[~,n]=size(X);
T1 = n;
s = 12;
[W] = denW1(X,T1,s);

rho =1.8;
lambda1 = 0.05;
lambda2 = 0.005;
lambda3 = 0.1;
% T1 = n;
T2 = 0.01;
K = 200; %% 200 220  200
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
% 
% acc(i,:)=[PSNR,accuracy1,accuracy2,accuracy3,accuracy4,accuracy5, accuracy6, accuracy7, accuracy8, accuracy9];
% NMI(i,:)=[n1,n2,n3,n4,n5,n6,n7,n8,n9];
% end

% figure,imagesc((abs(Z1)+abs(Z1'))/2),xlabel('SSC');
% figure,imagesc((abs(Z2)+abs(Z2'))/2),xlabel('LRR');
% figure,imagesc((abs(Z3)+abs(Z3'))/2),xlabel('OSC');
% figure,imagesc(W2),xlabel('TSC');
% figure,imagesc((abs(Z5)+abs(Z5'))/2),xlabel('BD-QOSC');
% figure,imagesc((abs(Z6)+abs(Z6'))/2),xlabel('SC-LRR');
% figure,imagesc((abs(Z7)+abs(Z7'))/2),xlabel('SSC-L1TG');
% figure,imagesc((abs(Z8)+abs(Z8'))/2),xlabel('FSC-NN');
% figure,imagesc((abs(Z81)+abs(Z81'))/2),xlabel('FeaMAC');
% figure,imagesc(A9),xlabel('feaSCLRR');
% 
% plotClusters(final_clusters1);xlabel('SSC');
% plotClusters(final_clusters2);xlabel('LRR');
% plotClusters(final_clusters3);xlabel('OSC');
% plotClusters(final_clusters4);xlabel('TSC');
% plotClusters(final_clusters5);xlabel('BD-QOSC');
% plotClusters(final_clusters6);xlabel('SC-LRR');
% plotClusters(final_clusters7);xlabel('SSC-L1TG');
% plotClusters(final_clusters8);xlabel('FSC-NN');
% plotClusters(final_clusters81);xlabel('FeaMAC');
% plotClusters(final_clusters9);xlabel('feaSCLRR');
