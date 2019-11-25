% Matthew Kobilas
% Professor Horacio Rotstein
% MATH 430-001
% 18 September 2019

% Solution of the logistic growth with a threshold equation
% Numerical method: modified Euler (Runge-Kutta, order 2)

clearvars;
close all;

% ODE
% Equation (2) from the homework
% X' = -r (1-X/T) ( 1 - X/K) X + I;
% r, K, T, and I are constants
% T and K represent the threshold (V_th) and saturation (V_sat) V levels
%   for the unbiased case (I = 0).
% I represents the current applied to the system
%   For I = 0, V_th = T and V_sat = K

% Plot 1 appears to be the Voltage (X) vs. Time (t)
% Plot 2 appears to be the Voltage (X) vs. F

Tmax = 100;
dt = 0.01;
t = 0:dt:Tmax;
r = 1;
T = 0.25;
K = 1;
logthr=@(x) -r*x.*(1-x/T).*(1-x/K);

I = 0;
x = zeros(1,length(t));
x(1) = 0.01;
for j=1:length(t)-1
    k1x = logthr(x(j))+I;
    ax = x(j)+k1x*dt;
    k2x = logthr(ax)+I;
    x(j+1)=x(j)+(k1x+k2x)*dt/2;
end
xx = -1:0.01:2;
figure(2)
hFig = figure(2);
set(hFig, 'Position', [40 400 1000 500]); 
subplot(1,2,1)
hold on
plot(t,x,'-b','linewidth',2);
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',20);
title('V(t): I = 0');
xlabel('t');
ylabel('X');
subplot(1,2,2)
hold on
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]);
plot(xx,logthr(xx)+I,'linewidth',2);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',20);
xlabel('X');
ylabel('F');
title('F(V): I = 0');

I = 0.05;
x = zeros(1,length(t));
x(1) = 0.01;
for j=1:length(t)-1
    k1x = logthr(x(j))+I;
    ax = x(j)+k1x*dt;
    k2x = logthr(ax)+I;
    x(j+1)=x(j)+(k1x+k2x)*dt/2;
end
xx = -1:0.01:2;
figure(3)
hFig = figure(3);
set(hFig, 'Position', [40 400 1000 500]); 
subplot(1,2,1)
hold on
plot(t,x,'-b','linewidth',2);
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',20);
title('V(t): I = 0.05');
xlabel('t');
ylabel('X');
subplot(1,2,2)
hold on
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]);
plot(xx,logthr(xx)+I,'linewidth',2);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',20);
xlabel('X');
ylabel('F');
title('F(V): I = 0.05');

I = 0.1;
x = zeros(1,length(t));
x(1) = 0.01;
for j=1:length(t)-1
    k1x = logthr(x(j))+I;
    ax = x(j)+k1x*dt;
    k2x = logthr(ax)+I;
    x(j+1)=x(j)+(k1x+k2x)*dt/2;
end
xx = -1:0.01:2;
figure(4)
hFig = figure(4);
set(hFig, 'Position', [40 400 1000 500]); 
subplot(1,2,1)
hold on
plot(t,x,'-b','linewidth',2);
plot([0 Tmax],[T T],'--');
plot([0 Tmax],[0 0]','--');
plot([0 Tmax],[K K],'--');
axis([0 Tmax -0.1 1.1]);
set(gca,'fontsize',20);
title('V(t): I = 0.1');
xlabel('t');
ylabel('X');
subplot(1,2,2)
hold on
plot([xx(1)-1 xx(end)+1],[0 0],'--','Color',[.6 .6 .6]);
plot(xx,logthr(xx)+I,'linewidth',2);
axis([-0.1 1.1 -0.5 0.5]);
set(gca,'fontsize',20);
xlabel('X');
ylabel('F');
title('F(V): I = 0.1');