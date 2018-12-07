close all;
clear;
clc;

%% 绘制原始图像
xh = (-1:0.01:1)';
runge = 1./(1+25*xh.^2);
figure
plot(xh,runge,'r','Linewidth',2)
xlabel('x')
ylabel('y')
title('Runge函数图像')
saveas(gcf, 'result/runge.png')

%% 等距Lagrange插值
x = (-1:0.1:1)';
y = 1./(1+25*x.^2);

xh2 = (-0.5:0.01:0.5)';
runge2 = 1./(1+25*xh2.^2);

lag = lagrange(x,y,xh);
figure
plot(xh,lag, 'b', 'Linewidth',2)
xlabel('x')
ylabel('y')
title('20次Lagranged等距节点插值Runge函数图像')
saveas(gcf, 'result/lagrunge.png')

lag2 = lagrange(x,y,xh2);
figure
plot(xh2,lag2, 'b', 'Linewidth',2)
hold on
plot(xh2,runge2, 'g', 'Linewidth',2)
xlabel('x')
ylabel('y')
legend('Lagrange插值图像','原图像')
title('20次Lagrange等距节点插值Runge函数局部图像与原图像对比')
saveas(gcf, 'result/lagrunge2.png')

%% 等距Newton插值
newton = newtonInterpol(x,y,xh);
figure
plot(xh,newton, 'g', 'Linewidth',2)
xlabel('x')
ylabel('y')
title('20次Newton等距节点插值Runge函数图像')
saveas(gcf, 'result/newtonrunge.png')

newton2 = newtonInterpol(x,y,xh2);
figure
plot(xh2,newton2, 'b', 'Linewidth',2)
hold on
plot(xh2,runge2, 'r', 'Linewidth',2)
xlabel('x')
ylabel('y')
legend('Newton插值图像','原图像')
title('20次Newton等距节点插值Runge函数局部图像与原图像对比')
saveas(gcf, 'result/newtonrunge2.png')

%% Chebyshev多项式零点Lagrange插值
k = (0:20)';
x = cos((2*k+1)/42*pi);
y = 1./(1+25*x.^2);

lagche = lagrange(x,y,xh);
figure
hold on
plot(xh,runge,'c','Linewidth',2)
plot(xh,lagche, 'r', 'Linewidth',2)
xlabel('x')
ylabel('y')
title('Chebyshev多项式零点插值Runge函数图像')
legend('原图像','插值图像')
saveas(gcf, 'result/lagche.png')

%% 分段线性插值
x = (-1:0.1:1)';
y = 1./(1+25*x.^2);

ylin = interp1(x,y,xh);
figure
plot(xh,ylin, 'b', 'Linewidth',2)
xlabel('x')
ylabel('y')
title('分段线性插值Runge函数图像')
saveas(gcf, 'result/linrunge.png')

%% 三次自然样条插值
x = (-1:0.1:1)';
y = 1./(1+25*x.^2);

yspline = interp1(x,y,xh,'spline');
figure
plot(xh,yspline, 'm', 'Linewidth',2)
xlabel('x')
ylabel('y')
title('三次自然样条插值Runge函数图像')
saveas(gcf, 'result/splinerunge.png')
