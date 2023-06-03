%% Van der Pol方程 
clear all;		
close all;
clc;
clf;

%% example4-1 ddx + 0.2*(x^2-1)*dx -x + x^3  = 0.53*cos(t)
%% 选择状态变量x1 = x, x2 =dx
Van = @(t,x)[x(2); -0.2*(x(1)^2-1)*x(2)+x(1)-x(1)^3+0.53*cos(t)]; 
x0 = [0.1; -0.2];
dt = 0.05;
t_end = 50;
[t, x] = ode45(Van, [0:dt:t_end], x0);
x1 = x(:,1);
x2 = x(:,2);
figure(1);
plot(t,x1);
figure(2);
plot(x1,x2);


%% example4-2 ddx + 0.3^2*(x + 0.2^2*x^3)  = 0.2*sin(2t)
%% 选择状态变量x1 = x, x2 =dx
% Van = @(t,x)[x(2); -0.3^2*(x(1) + 0.2^2*x(1)^3) + 0.2*sin(2*t)]; 
% x0 = [0.15; 0];
% dt = 0.05;
% t_end = 45;
% [t, x] = ode45(Van, [0:dt:t_end], x0);
% x1 = x(:,1);
% x2 = x(:,2);
% figure(1);
% plot(t,x1);
% figure(2);
% plot(x1,x2);

%% example4-3 ddx + 0.1*(x^2-1)*dx + 0.5x + 0.5x^3  = 0.5*cos(0.79t)
%% 选择状态变量x1 = x, x2 =dx
% Van = @(t,x)[x(2); -0.1*(x(1)^2-1)*x(2)-0.5*x(1)-0.5*x(1)^3+0.5*cos(0.79*t)]; 
% x0 = [0; 0];
% dt = 0.05;
% t_end = 50;
% [t, x] = ode45(Van, [0:dt:t_end], x0);
% x1 = x(:,1);
% x2 = x(:,2);
% figure(1);
% plot(t,x1);
% figure(2);
% plot(x1,x2);


%% example4-4 ddx + 0.5*dx -x + x^3  = 0.8725*cos(t)+0.0005*cos(t)
%% 选择状态变量x1 = x, x2 =dx
% Van = @(t,x)[x(2); -0.5*x(2)+ x(1)- x(1)^3+0.8725*cos(t)+0.0005*cos(t)]; 
% x0 = [1; 1];
% dt = 0.05;
% t_end = 500;
% [t, x] = ode45(Van, [0:dt:t_end], x0);
% x1 = x(:,1);
% x2 = x(:,2);
% figure(1);
% plot(t,x1);
% figure(2);
% plot(x1,x2);

%% example4-5 ddx + 5*(x^2-1)*dx +x + 0.01*x^3  = 4.9*cos(2.463*t)  
%% 选择状态变量x1 = x, x2 =dx
Van = @(t,x)[x(2);  -5*(x(1)^2-1)*x(2)- x(1)- 0.01*x(1)^3+4.9*cos(2.463*t)]; 
x0 = [0.1; 0.1];
dt = 0.01;
t_end = 300;
[t, x] = ode45(Van, [0:dt:t_end], x0);
x1 = x(:,1);
x2 = x(:,2);
figure(1);
plot(t,x1);
figure(2);
plot(x1,x2);











%% 参数值
% alpha = 16;
% beta = 15;
% gamma = 0.5;
% xi = 1.4;
% a = 0.2;
% b = 0.4;
% 
% %% 微分方程
% Chua = @(t,x)[alpha*(x(2)-x(1)+xi*x(1)-(a+3*b*x(4)^2)*x(1)); x(1)-x(2)+x(3); -beta*x(2)-gamma*x(3); x(1)];
% 
% %% 解微分方程
% x0 = [0;10^(-10);0;0];
% dt = 0.001;
% t_end = 10;
% 
% [t, x] = ode45(Chua, [0:dt:t_end], x0);
% % t = [0:dt:t_end]';
% % x = ode4(Chua, [0:dt:t_end], x0);
% 
% x1 = x(:,1);
% ii = x1 ./ -1.4;
% x2 = x(:,2);
% x3 = x(:,3);
% x4 = x(:,4);
% 
% figure(1);
% plot(t,x3); % plot x3
% axis([0 10 -3 3]);
%figure(2);
% plot(t,-beta*x2-gamma*x3); % plot  x3dot = -beta*x2-gamma*x3
% %axis([0 10 -3 3]);
% figure(3);
% i = (a+3*b*x4.^2).*x1;
% v = x1;
% plot(v,i);

% figure(1)
% plot(x1,x2);
% xlabel('\itx\rm_1','FontSize',12,'FontName','Times New Roman');
% ylabel('\itx\rm_2','FontSize',12,'FontName','Times New Roman');
% 
% figure(2)
% plot(x1,x4);
% xlabel('\itx\rm_1','FontSize',12,'FontName','Times New Roman');
% ylabel('\itx\rm_4','FontSize',12,'FontName','Times New Roman');

% figure(3)
% plot(x1,ii);
% xlabel('\itv\rm','FontSize',12,'FontName','Times New Roman');
% ylabel('\iti\rm','FontSize',12,'FontName','Times New Roman');

% save chua.mat x

% %% %% Van der Pol方程 ddy + mu*(y^2-1)*dy + y = 0
%    %% 选择状态变量x1 = y x2 =dy
% Van = @(t,x,mu)[x(2); -mu*(x(1)^2-1)*x(2)-x(1)];
% x0 = [-0.1; -0.5];
% dt = 0.01;
% t_end = 50;
% mu = 1;
% [t, x] = ode45(Van, [0:dt:t_end], x0, [], mu);
% x1 = x(:,1);
% x2 = x(:,2);
% figure;
% plot(x1,x2);


