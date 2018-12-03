close all;
clear;
clc;

n=2;
H=Hilbert(n);
x=ones(n,1);
b=H*x;

%Gauss消去法
fprintf(1,'Gauss消去法相对误差:')
x1=Gauss(H,b);
error1=norm(x1-x)/norm(x);
disp(error1);

%Jacobi迭代法
x0 = zeros(n,1);
[x2,k] = Jacobi(H,b,x0,1e-6);
error2=norm(x2-x)/norm(x);
disp(error2);disp(k);