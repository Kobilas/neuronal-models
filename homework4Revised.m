% Matthew Kobilas
% Professor Horacio Rotstein
% MATH 430-001
% 22 October 2019

close all;

% Question 1 Functions
% Steady-state (in)activation functions
m_infin = @(V) 1 ./ (1 + exp(-((V + 40) ./ 9)));
h_infin = @(V) 1 ./ (1 + exp((V + 62) ./ 10));
n_infin = @(V) 1 ./ (1 + exp(-((V + 53) ./ 16)));
% Voltage-dependent time constants
tau_m = @(V) 0.3;
tau_h = @(V) 1 + (11 ./ (1 + exp((V + 62) ./ 10)));
tau_n = @(V) 1 + (6 ./ (1 + exp((V + 53) ./ 16)));

figureNum = 1;

% Question 1
% Part (a)
figure(figureNum);
figureNum = figureNum + 1;
hold on;
fplot(m_infin, 'b', 'linewidth', 2);
fplot(h_infin, 'r', 'linewidth', 2);
fplot(n_infin, 'g', 'linewidth', 2);
xlabel("V (mV)");
title("Steady-state (in)activation functions");
legend("m_\infty", "h_\infty", "n_\infty");
axis([-50 100 0 1]);
hold off;

% Part (b)
warning('off', 'MATLAB:fplot:NotVectorized');
figure(figureNum);
figureNum = figureNum + 1;
hold on;
fplot(tau_m, 'b', 'linewidth', 2);
fplot(tau_h, 'r', 'linewidth', 2);
fplot(tau_n, 'g', 'linewidth', 2);
xlabel("V (mV)");
title("Voltage-dependent time constants");
legend("\tau_m", "\tau_h", "\tau_n");
axis([-50 100 0 5]);
hold off;
warning('on', 'MATLAB:fplot:NotVectorized');

% Part (c)
params = 0:0.25:4;
numSpikes = zeros(1, length(params));
for k=1:length(params)
    C = 1;
    E_L = -52;
    E_Na = 55;
    E_K = -75;
    G_L = 0.3;
    G_Na = 120;
    G_K = 36;
    m_Na = 3;
    h_Na = 1;
    n_K = 4;

    Tmax = 4000;
    dt = 0.1;
    t = 0:dt:Tmax;

    ti = 1000;
    tf = 3000;
    H = zeros(1,length(t));
    H(floor(ti/dt):floor(tf/dt)) = params(k);

    m_dot = zeros(1, length(t));
    m_dot(1) = 0;
    h_dot = zeros(1, length(t));
    h_dot(1) = 0;
    n_dot = zeros(1, length(t));
    n_dot(1) = 0;

    V = zeros(1,length(t));
    V(1) = E_L;
    spikeSwitch = true;

    for i=1:length(t)-1
        kv1 = (H(i) - (G_K * n_dot(i)^n_K * (V(i) - E_K)) ...
            - (G_Na * m_dot(i)^m_Na * h_dot(i)^h_Na * (V(i) - E_Na)) ...
            - (G_L * (V(i) - E_L))) / C;
        av = V(i) + kv1 * dt;
        kv2 = (H(i+1) - (G_K * n_dot(i+1)^n_K * (av - E_K)) ...
            - (G_Na * m_dot(i+1)^m_Na * h_dot(i+1)^h_Na * (av - E_Na)) ...
            - (G_L * (av - E_L))) / C;
        V(i+1) = V(i) + (kv1+kv2)*dt/2;
        
        if V(i) > 0 && spikeSwitch
            numSpikes(k) = numSpikes(k) + 1;
            spikeSwitch = false;
        elseif V(i) < 0
            spikeSwitch = true;
        end
        numSpikes(k) = 1/(numSpikes(k)/2000);
        
        k1 = (m_infin(V(i)) - m_dot(i)) / tau_m(V(i));
        ax = m_dot(i) + k1 * dt;
        k2 = (m_infin(V(i+1)) - ax) / tau_m(V(i+1));
        m_dot(i+1) = m_dot(i) + (k1+k2)*dt/2;

        k1 = (h_infin(V(i)) - h_dot(i)) / tau_h(V(i));
        ax = h_dot(i) + k1 * dt;
        k2 = (h_infin(V(i+1)) - ax) / tau_h(V(i+1));
        h_dot(i+1) = h_dot(i) + (k1+k2)*dt/2;

        k1 = (n_infin(V(i)) - n_dot(i)) / tau_n(V(i));
        ax = n_dot(i) + k1 * dt;
        k2 = (n_infin(V(i+1)) - ax) / tau_n(V(i+1));
        n_dot(i+1) = n_dot(i) + (k1+k2)*dt/2;
    end
    if mod(params(k), 1) == 0
        figure(figureNum);
        figureNum = figureNum + 1;
        hold on;
        plot(t, V, 'b', 'linewidth', 2);
        set(figure(figureNum-1), 'Position', [40 400 1000 500]);
        axis([0 Tmax E_K E_Na]);
        set(gca,'fontsize',24);
        xlabel('t (ms)');
        ylabel('V (mV)');
        title(strcat("Hodgkin-Huxley; I_a_p_p = ", num2str(params(k))));
        hold off;
    end
end

% Part (d)
figure(figureNum);
figureNum = figureNum + 1;
hold on;
plot(params, numSpikes, 'b', 'linewidth', 2);
set(gca, 'fontsize', 24);
xlabel('I_a_p_p Values');
ylabel('Frequency');
title("I_a_p_p vs. Spike Freq.");
hold off;

% Question 2
% Part (a)
V = linspace(-150, 150);
warning('off', 'MATLAB:fplot:NotVectorized');
figure(figureNum);
figureNum = figureNum + 1;
hold on;
plot(V, m_infin_rsrch(V), 'b', 'linewidth', 2);
plot(V, h_infin_rsrch(V), 'r', 'linewidth', 2);
plot(V, n_infin_rsrch(V), 'g', 'linewidth', 2);
xlabel("V (mV)");
title("Paper's (in)activation functions");
legend("m_\infty", "h_\infty", "n_\infty");
axis([-150 150 -1 1]);
hold off;

% Part (b)
figure(figureNum);
figureNum = figureNum + 1;
hold on;
plot(V, tau_m_rsrch(V), 'b', 'linewidth', 2);
plot(V, tau_h_rsrch(V), 'r', 'linewidth', 2);
plot(V, tau_n_rsrch(V), 'g', 'linewidth', 2);
xlabel("V (mV)");
title("Voltage-dependent time constants");
legend("\tau_m", "\tau_h", "\tau_n");
axis([-150 150 -5 100]);
hold off;
warning('on', 'MATLAB:fplot:NotVectorized');

% Part (c)
%numSpikes = zeros(1, length(params));
an = @(V) (0.032.*(V+52))./(1-(exp(-(V+52)./5)));
bn = @(V) (0.5.*(exp(-(57+V)./(40))));
ninf = @(V) an(V)./(an(V)+bn(V));
am = @(V) (0.32.*(54+V))./(1-exp(-(V+54)./4));
bm = @(V) (0.28.*(V+27))./(exp((V+27)./5)-1);
minf = @(V) am(V)./(am(V)+bm(V));
ah = @(V) 0.128.*exp(-(50+V)./18);
bh = @(V) 4./(1+exp(-(V+27)./5));
hinf = @(V) ah(V)./(ah(V)+bh(V));
taun = @(V) 1./(an(V)+bn(V));
taum = @(V) 1./(am(V)+bm(V));
tauh = @(V) 1./(ah(V)+bh(V));
dt = 0.05;
x = 0:dt:10000;
C = 1;
E_L = -67;
E_Na = 50;
E_K = -100;
G_Na = 100;
G_K = 80;
G_L = 0.1;
step = (0.12 - 0.1190) / 11;
Iapp0(1:12) = 0.1190:step:0.12;
V = zeros(1,length(x));
n = zeros(1,length(x));
h = zeros(1,length(x));
m = zeros(1,length(x));
V(1) = -70;

figureCount = figureNum;
spikeCount = zeros(1, length(Iapp0));

for j=1:length(Iapp0)-1
    for i=1:length(x)-1
        dV1t1 = Iapp0(j)-(G_K*n(i)^4*(V(i)-E_K))-(G_Na*(m(i)^3)*(h(i))*(V(i)-E_Na))-(G_L*(V(i)-E_L));
        dn1t1 = (ninf(V(i))-n(i))/taun(V(i));
        dm1t1 = (minf(V(i))-m(i))/taum(V(i));
        dh1t1 = (hinf(V(i))-h(i))/tauh(V(i));

        V1t1 = V(i)+dV1t1*dt;
        n1t1 = n(i)+dn1t1*dt;
        m1t1 = m(i)+dm1t1*dt;
        h1t1 = h(i)+dh1t1*dt;

        dV1t2 = Iapp0(j)-(G_K*((n1t1)^4)*(V1t1-E_K))-(G_Na*((m1t1)^3)*h1t1*(V1t1-E_Na))-G_L*(V1t1-E_L);
        dn1t2 = ((ninf(V1t1)-n1t1)/taun(V1t1));
        dm1t2 = ((minf(V1t1)-m1t1)/taum(V1t1));
        dh1t2 = ((hinf(V1t1)-h1t1)/tauh(V1t1));

        V(i+1) = V(i)+(dV1t1+dV1t2)/2*dt;
        n(i+1) = n(i)+(dn1t1+dn1t2)/2*dt;
        m(i+1) = m(i)+(dm1t1+dm1t2)/2*dt;
        h(i+1) = h(i)+(dh1t1+dh1t2)/2*dt;
    end
    figure(figureCount);
    figureCount = figureCount + 1;
    plot(x, V, 'linewidth', 2);
    set(figure(figureCount-1), 'Position', [40 400 1000 500]);
    axis([0 10000 E_K E_Na]);
    title(strcat("I_a_p_p = ", num2str(Iapp0(j))));
    peaks = findpeaks(V);
    spikeCount(j) = length(peaks)/10;
end
figure(figureCount);
figureCount = figureCount + 1;
plot(Iapp0, spikeCount, 'linewidth', 2);
title("I_a_p_p vs. Frequency");
xlabel("I_a_p_p (mV)");
ylabel("Frequency (spikes/s)");

an = @(V) (0.032.*(V+52))./(1-(exp(-(V+52)./5)));
bn = @(V) (0.5.*(exp(-(57+V)./(40))));
ninf = @(V) an(V)./(an(V)+bn(V));
am = @(V) (0.32.*(54+V))./(1-exp(-(V+54)./4));
bm = @(V) (0.28.*(V+27))./(exp((V+27)./5)-1);
minf = @(V) am(V)./(am(V)+bm(V));
ah = @(V) 0.128.*exp(-(50+V)./18);
bh = @(V) 4./(1+exp(-(V+27)./5));
hinf = @(V) ah(V)./(ah(V)+bh(V));
taun = @(V) 1./(an(V)+bn(V));
taum = @(V) 1./(am(V)+bm(V));
tauh = @(V) 1./(ah(V)+bh(V));
dt = 0.05;
x = 0:dt:10000;
C = 1;
E_L = -67;
E_Na = 50;
E_K = -100;
G_Na = 100;
G_K = 80;
G_L = 0.1;
step = (0.15 - 0.12) / 11;
Iapp0(1:12) = 0.12:step:0.15;
V = zeros(1,length(x));
n = zeros(1,length(x));
h = zeros(1,length(x));
m = zeros(1,length(x));
V(1) = -70;

spikeCount = zeros(1, length(Iapp0));

for j=1:length(Iapp0)-1
    for i=1:length(x)-1
        dV1t1 = Iapp0(j)-(G_K*n(i)^4*(V(i)-E_K))-(G_Na*(m(i)^3)*(h(i))*(V(i)-E_Na))-(G_L*(V(i)-E_L));
        dn1t1 = (ninf(V(i))-n(i))/taun(V(i));
        dm1t1 = (minf(V(i))-m(i))/taum(V(i));
        dh1t1 = (hinf(V(i))-h(i))/tauh(V(i));

        V1t1 = V(i)+dV1t1*dt;
        n1t1 = n(i)+dn1t1*dt;
        m1t1 = m(i)+dm1t1*dt;
        h1t1 = h(i)+dh1t1*dt;

        dV1t2 = Iapp0(j)-(G_K*((n1t1)^4)*(V1t1-E_K))-(G_Na*((m1t1)^3)*h1t1*(V1t1-E_Na))-G_L*(V1t1-E_L);
        dn1t2 = ((ninf(V1t1)-n1t1)/taun(V1t1));
        dm1t2 = ((minf(V1t1)-m1t1)/taum(V1t1));
        dh1t2 = ((hinf(V1t1)-h1t1)/tauh(V1t1));

        V(i+1) = V(i)+(dV1t1+dV1t2)/2*dt;
        n(i+1) = n(i)+(dn1t1+dn1t2)/2*dt;
        m(i+1) = m(i)+(dm1t1+dm1t2)/2*dt;
        h(i+1) = h(i)+(dh1t1+dh1t2)/2*dt;
    end
    figure(figureCount);
    figureCount = figureCount + 1;
    plot(x, V, 'linewidth', 2);
    set(figure(figureCount-1), 'Position', [40 400 1000 500]);
    axis([0 10000 E_K E_Na]);
    title(strcat("I_a_p_p = ", num2str(Iapp0(j))));
    peaks = findpeaks(V);
    spikeCount(j) = length(peaks)/10;
end
figure(figureCount);
figureCount = figureCount + 1;
plot(Iapp0, spikeCount, 'linewidth', 2);
title("I_a_p_p vs. Frequency");
xlabel("I_a_p_p (mV)");
ylabel("Frequency (spikes/s)");

% Question 2 Functions
% Steady-state (in)activation functions
function research_m_infinity = m_infin_rsrch(V)
    alpha = (0.32 .* (V + 54)) ./ (1 - exp(-(V + 54)) ./ 4);
    beta = (0.28 .* (V + 27)) ./ ((exp((V + 27) / 5)) - 1);
    research_m_infinity = alpha ./ (alpha + beta);
end
function research_h_infinity = h_infin_rsrch(V)
    alpha = (0.128 .* exp((-(50 + V)) / 18));
    beta = 4 ./ (1 + exp((-(V + 27)) / 5));
    research_h_infinity = alpha ./ (alpha + beta);
end
function research_n_infinity = n_infin_rsrch(V)
    alpha = (0.032 .* (V + 52)) ./ (1 - exp(-(V + 52) ./ 5));
    beta = (0.5 .* exp(-(57 + V) ./ 40));
    research_n_infinity = alpha ./ (alpha + beta);
end
% Voltage-dependent time constants
function research_m_tau = tau_m_rsrch(V)
    alpha = (0.32 .* (V + 54)) ./ (1 - exp(-(V + 54)) ./ 4);
    beta = (0.28 .* (V + 27)) ./ ((exp((V + 27) / 5)) - 1);
    research_m_tau = 1 ./ (alpha + beta);
end
function research_h_tau = tau_h_rsrch(V)
    alpha = (0.128 .* exp((-(50 + V)) / 18));
    beta = 4 ./ (1 + exp((-(V + 27)) / 5));
    research_h_tau = 1 ./ (alpha + beta);
end
function research_n_tau = tau_n_rsrch(V)
    alpha = (0.032 .* (V + 52)) ./ (1 - exp(-(V + 52) ./ 5));
    beta = (0.5 .* exp(-(57 + V) ./ 40));
    research_n_tau = 1 ./ (alpha + beta);
end