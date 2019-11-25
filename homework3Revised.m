% Matthew Kobilas
% Professor Horacio Rotstein
% MATH 430-001
% 6 October 2019

close all;
currFigure = 1;

C = 1;
El = -60;
Gl = 0.1;
Vzero = -60;

% Problem 1-a
R = 1/Gl;
tau = R*C;
fprintf("Tau=RC is equal to %1.4f milliseconds\n", tau);

% Problem 1-b
params = [-0.5 0.5];
Tmax = 1000;
dt = 0.1;
t = 0:dt:Tmax;
ti = dt;
tf = 10000;
for k=1:length(params)
    Iapp = params(k);
    Vinf = R*Iapp;

    H = zeros(1,length(t));
    H(floor(ti/dt):floor(tf/dt))=1;

    V = zeros(1,length(t));
    V(1) = El;

    for j=1:length(t)-1
        kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
        av = V(j)+kv1*dt;
        kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
        V(j+1) = V(j) + (kv1+kv2)*dt/2;
    end

    absErr = zeros(1, length(t));
    analyt = zeros(1, length(t));
    for i=1:length(t)
        answr = analyticalPme(i, El, Vzero, tau, Vinf);
        absErr(i) = abs(V(i) - answr);
        analyt(i) = answr;
    end

    figure(currFigure)
    currFigure = currFigure + 1;
    subplot(1, 2, 1)
    hold on
    plot(t,V,'b','linewidth',2);
    plot(t, analyt, 'r', 'linewidth', 2);
    axis([0 Tmax -66 -54]);
    set(gca,'fontsize',20);
    xlabel('t (ms)');
    ylabel('V (mV)');
    title(strcat('PME Iapp = ', num2str(Iapp, 2)));
    subplot(1, 2, 2)
    hold on
    plot(t, absErr, 'g', 'linewidth', 2);
    axis([0 Tmax -1 5]);
    set(gca, 'fontsize', 20);
    xlabel('t (ms)');
    ylabel('Absolute Error (x 10^-^4)');
    title(strcat('AbsErr Iapp = ', num2str(Iapp, 2)));
    hold off
end

% Problem 1-c

Tmax = 400;
dt = 0.1;
t = 0:dt:Tmax;

ti = 100;
tf = 200;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;
Iapp = 0.5;

V = zeros(1,length(t));
V(1) = El;

for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;
end

heaviAnalyt = El + Vinf * ((heaviside(t-ti).*(1-exp(-1*(t-ti/tau)))) - ...
    (heaviside(t-tf).*(1-exp(-1*(t-tf/tau)))));

figure(currFigure);
currFigure = currFigure + 1;
subplot(1, 2, 1);
hold on;
plot(t,V,'b','linewidth',2);
plot(t, heaviAnalyt, 'r', 'linewidth', 2);
axis([0 Tmax -80 0]);
set(gca,'fontsize',20);
xlabel('t');
ylabel('V');
title("Square Pulse of current");
subplot(1, 2, 2);
hold on;
plot(t, abs(V-heaviAnalyt), 'g', 'linewidth', 2);
axis([0 Tmax -20 20])
set(gca, 'fontsize', 20);
xlabel('t (ms)');
ylabel('Absolute Error (mV)');
title('AbsErr Iapp = Heaviside');
hold off

% Problem 1-d
params = [1 5 10 20 100];
Iapp0 = 5;
Tmax = 1000;
dt = 0.1;
t = 0:dt:Tmax;
ti = dt;
tf = 10000;
ampCount = zeros(1, length(params));
maxAmp = [-200 -200 -200 -200 -200];
for k=1:length(params)
    w = params(k);

    H = zeros(1, length(t));
    H(floor(ti/dt):floor(tf/dt))=1;

    Iapp = Iapp0*sin(2*pi*w*t/1000);
    ampSwitch = true;

    for i=1:length(t)-1
        kv1 = (-Gl*(V(i)-El)+Iapp(i)*H(i))/C;
        av = V(i)+kv1*dt;
        kv2 = (-Gl*(av-El)+Iapp(i+1)*H(i+1))/C;
        V(i+1) = V(i) + (kv1+kv2)*dt/2;
        if (V(i) >= -55) && ampSwitch
            ampCount(k) = ampCount(k) + 1;
            ampSwitch = false;
        elseif V(i) < -60
            ampSwitch = true;
        end
        if (V(i) > maxAmp(k))
            maxAmp(k) = V(i);
        end
    end

    figure(currFigure);
    currFigure = currFigure + 1;
    hold on;
    plot(t, V, 'b', 'linewidth', 2);
    plot(t, Iapp, 'r', 'linewidth', 2);
    axis([0 Tmax -120 120]);
    set(gca, 'fontsize', 20);
    xlabel('t');
    ylabel('V');
    title(strcat('Sine wave resp w=', num2str(w)));
    hold off
end

figure(currFigure);
currFigure = currFigure + 1;
hold on;
plot(params, ampCount, 'b', 'linewidth', 2);
set(gca, 'fontsize', 20);
xlabel('Input Frequency');
ylabel('Output Frequency');
title("Input vs. Output Frequency");
hold off;

figure(currFigure);
currFigure = currFigure + 1;
hold on;
plot(params, maxAmp, 'b', 'linewidth', 2);
set(gca, 'fontsize', 20);
xlabel('Input Frequency');
ylabel('Output Amplitude');
title('Input Frequency vs. Output Amplitude');
hold off;

% Problem 2-a
params = [0.5 1 1.01 2];
C = 1;
El = -60;
Gl = 0.1;
Vth = -50;
Vrst = -65;
T = 1000;
dt = 0.1;
t = 0:dt:T;
ti = dt;
tf = 10000;
for k=1:length(params)
    Iapp = params(k);
    H = zeros(1, length(t));
    H(floor(ti/dt):floor(tf/dt))=1;
    V = zeros(1, length(t));
    V(1) = El;
    for i=1:length(t)-1
        kv1 = (-Gl*(V(i)-El)+Iapp*H(i))/C;
        av = V(i)+kv1*dt;
        kv2 = (-Gl*(av-El)+Iapp*H(i+1))/C;
        V(i+1) = V(i) + (kv1+kv2)*dt/2;
        if V(i+1) >= Vth
            V(i) = 60;
            V(i+1) = Vrst;
        end
    end
    figure(currFigure);
    currFigure = currFigure + 1;
    hold on;
    plot(t, V, 'r', 'linewidth', 2);
    axis([0 T -80 80]);
    set(gca, 'fontsize', 20);
    xlabel('t');
    ylabel('V');
    title(strcat('Integrate and Fire Model Iapp=', num2str(Iapp)));
    R = 1/Gl;
    tau = R*C;
    Vinf = R*Iapp;
    rIsi = (1/tau)*((log(((Vrst-Vinf-El)/(Vth-Vinf-El))))^-1);
    if (Vrst-Vinf-El)/(Vth-Vinf-El) <= 0
        fprintf("r_isi = 0 spikes/ms at %1.4f microAmps\n", Iapp);
    else
        fprintf("r_isi = %1.4f spikes/ms at %1.4f microAmps\n", rIsi*1000, Iapp);
    end
    hold off
end

% Problem 2-b
params = [10 100];
C = 1;
El = -60;
Ek = -85;
Gl = 0.1;
Iapp = 2;
Vth = -50;
Vrst = -65;
T = 1000;
dt = 0.1;
t = 0:dt:T;
ti = dt;
tf = 10000;
for k=1:length(params)
    H = zeros(1,length(t));
    H(floor(ti/dt):floor(tf/dt))=1;
    V = zeros(1,length(t));
    V(1) = El;
    gsra = zeros(1, length(t));
    gsra(1) = 0;
    deltaGsra = 0.1;
    tauSra = params(k);

    for i=1:length(t)-1
        kg1 = (-gsra(i)/tauSra);
        kv1 = (((-V(i)+El) - R .* gsra(i) .* (V(i)-Ek) + R .* Iapp .* H(i)) / tauSra);
        
        ag = gsra(i) + kg1 * dt;
        av = V(i) + kv1 * dt;
        
        kg2 = (-ag / tauSra);
        kv2 = (((-av + El) - R .* ag .* (av-Ek) + R .* Iapp .* H(i+1))/tauSra);
        
        V(i+1) = V(i) + (kv1+kv2)*dt/2;
        gsra(i+1) = gsra(i) + (kg1+kg2)*dt/2;
        
        if V(i+1) >= Vth
            V(i) = 60;
            V(i+1) = Vrst;
            gsra(i+1) = gsra(i) + deltaGsra;
        end
    end

    figure(currFigure)
    currFigure = currFigure + 1;
    subplot(2,1,1);
    hold on
    plot(t,V,'r','linewidth',2);
    axis([0 T -80 80]);
    set(gca,'fontsize',20);
    xlabel('t')
    ylabel('V');
    title("Spike Rate Adaptation");
    subplot(2,1,2);
    hold on
    plot(t, gsra, 'r', 'linewidth', 2);
    axis([0 T 0 0.15]);
    set(gca, 'fontsize', 20);
    xlabel('t');
    ylabel('G_S_R_A');
    title(strcat('Tau_S_R_A = ', num2str(tauSra), 'ms'));
end

function V = analyticalPme(t, El, Vzero, tau, Vinf)
    V = Vinf + El + (Vzero - El - Vinf)*exp(-t/tau);
end