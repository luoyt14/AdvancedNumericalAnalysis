close all;
clear;
clc;

%% 绘制原始图像
xh = (-1:0.01:1)';
yt = f(xh);
figure
plot(xh,yt,'r','Linewidth',1)
xlabel('x')
ylabel('y')
title('f(x)函数图像')
saveas(gcf, 'result/f.png')

%% Newton插值
x = (-1:0.1:1)';
y = f(x);

xh2 = (-0.5:0.01:0.5)';
runge2 = f(xh2);

lag = lagrange(x,y,xh);
newton = newtonInterpol(x,y,xh);
figure
plot(xh,newton, 'g', 'Linewidth',1)
xlabel('x')
ylabel('y')
title('20次Newton等距节点插值f(x)图像')
saveas(gcf, 'result/newtonf.png')

newton2 = newtonInterpol(x,y,xh2);
figure
plot(xh2,newton2, 'b', 'Linewidth',1)
hold on
plot(xh2,runge2, 'r', 'Linewidth',1)
xlabel('x')
ylabel('y')
legend('Newton插值图像','原图像')
title('20次Newton等距节点插值f(x)局部图像与原图像对比')
saveas(gcf, 'result/newtonf2.png')

% Chebyshev多项式零点Lagrange插值
k = (0:20)';
x = cos((2*k+1)/42*pi);
y = f(x);

lagche = lagrange(x,y,xh);
errorlag = abs(lagche-yt)./abs(yt+1e-6);
figure
hold on
plot(xh,yt,'r','Linewidth',1)
plot(xh,lagche, 'c', 'Linewidth',1)
xlabel('x')
ylabel('y')
title('Chebyshev多项式零点插值f(x)图像')
legend('原图像','插值图像')
saveas(gcf, 'result/lagchef.png')

%% 分段Chebyshev多项式零点插值
t = (0:9)';
xl = (cos((2*t+1)/20*pi)-1)/2;
yl = f(xl);
k = (0:10)';
xr = (cos((2*k+1)/22*pi)+1)/2;
yr = f(xr);
x = [xl;xr];
y = [yl;yr];

lagchel = lagrange(xl,yl,xh(1:100));
lagcher = lagrange(xr,yr,xh(101:201));
lagche=[lagchel;lagcher];
errorlag1 = abs(lagche-yt)./abs(yt+1e-6);
figure
hold on
plot(xh,yt,'r','Linewidth',1)
plot(xh,lagche, 'c', 'Linewidth',1)
xlabel('x')
ylabel('y')
title('分段Chebyshev多项式零点插值f(x)图像')
legend('原图像','插值图像')
saveas(gcf, 'result/lagchef1.png')

%% 分段线性插值
x = (-1:0.1:1)';
y = f(x);

ylin = interp1(x,y,xh);
figure
hold on
plot(xh,yt,'r','Linewidth',1)
plot(xh,ylin, 'b', 'Linewidth',1)
xlabel('x')
ylabel('y')
title('分段线性插值f(x)图像')
legend('原图像','插值图像')
saveas(gcf, 'result/linf.png')

%% 三次自然样条插值
x = (-1:0.1:1)';
y = f(x);

yspline = interp1(x,y,xh,'spline');
figure
hold on
plot(xh,yt,'b','Linewidth',1)
plot(xh,yspline, 'm', 'Linewidth',1)
xlabel('x')
ylabel('y')
title('三次自然样条插值f(x)图像')
legend('原图像','插值图像')
saveas(gcf, 'result/splinef.png')

%% 以间断点为界，左右两部分分别进行三次自然样条插值
x1 = (-1:0.1:-0.1)';
y1 = f(x1);
xh1 = (-1:0.01:-0.01)';
yspline1 = interp1(x1,y1,xh1,'spline');
x2 = (0:0.1:1)';
y2 = f(x2);
xh2 = (0:0.01:1)';
yspline2 = interp1(x2,y2,xh2,'spline');
yspline=[yspline1;yspline2];
xh = [xh1;xh2];
yt = f(xh);

figure
hold on
plot(xh,yt,'b','Linewidth',1)
plot(xh,yspline, 'm', 'Linewidth',1)
xlabel('x')
ylabel('y')
title('不跨过间断点的三次自然样条插值f(x)图像')
legend('原图像','插值图像')
saveas(gcf, 'result/splinef1.png')