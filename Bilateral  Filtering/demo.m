%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%主函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 读取
clear
load W1
X = W1;
A = diag(X);
f = A;
% f = imread('pflower.jpg');
 
% 设置参数
r = 5;% 滤波半径
a = 3;% 全局方差
b = 0.1;% 局部方差
 
% 判断二维图还是三维图
if ismatrix(f)
    g = bfilt_gray(f,r,a,b);
else
    g = bfilt_rgb(f,r,a,b);
end
 
% 显示
% subplot(121)
% imshow(f)
% subplot(122)
% imshow(g)
 
 v=diag(g)';
 

