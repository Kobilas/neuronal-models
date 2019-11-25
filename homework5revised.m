% Matthew Kobilas
% Professor Horacio Rotstein
% MATH 430-001
% 6 November 2019

close all;

% Question 1
% (a)
rho = 2.5;
R_a = 200;
R_m = 20000;
C_m = 1;
tau = C_m * R_m;
% lambda = sqrt((a*R_m)/(2*R_l))
lambda = sqrt((rho*R_m)/(2*R_a));
fprintf("Question 1, Part a:\n");
fprintf("tau: %d\nlambda: %3.2f\n", tau, lambda);
% (b)
d = 4;
R_a = 100;
R_m = 10000;
I = 0.25;
V_rest = -65;
lambda = sqrt(((d/2)*R_m)/(2*R_a));
figure(1);
hold on;
fplot(@(x) stdyState(lambda, R_a, d/2, I, V_rest, x), [0 60], 'b', 'linewidth', 2);
set(gca, 'fontsize', 24);
xlabel('Distance');
ylabel('Voltage');
title("Steady State Infinite Cable");
hold off;
fprintf("Question 1, Part b:\n");
fprintf("The plot of V_s(x) - V_rest will not reach 63 percent\n");
fprintf("of its maximum value since there is a horizontal\n");
fprintf("asymptote at y = 65 mV.\n");

% Question 2
% (a)
C_m = 1;
R_m = 3333;
R_a = 100;
E_L = -65;
tau = C_m * R_m;
fprintf("Question 2, Part a:\n");
fprintf("tau: %d\n", tau);
% (b)
rho = 1;
L = 0.159;
A = 2 * pi * rho * L;
AC_m = A * C_m;
fprintf("Question 2, Part b:\n");
fprintf("Area: %3.3f\nMembrane Capacitance: %3.3f\n", A, AC_m);
% (c)
rho_1 = 10^-4;
rho_2 = rho_1;
R_a = 100;
L_1 = 0.159;
L_2 = L_1;
C_m = 10^-6;
g_12 = (rho_1 * rho_2^2) / (R_a * L_1 * (((rho_2^2) * L_1) + ((rho_1^2) * L_2)));
comp = g_12 / C_m;
fprintf("Question 2, Part c:\n");
fprintf("Ratio of D_jk/C_m is %.4f 1/cm^2-ohm\n", comp);
% (d)
C_m = 1;
R_m = 3333;
E_L = -65;
rho = 10^-4;
L = 0.159;
A = pi * rho^2 * L;
g = (rho_1 * rho_2^2) / (R_a * L_1 * (((rho_2^2) * L_1) + ((rho_1^2) * L_2)));
Tmax = 1000;
dt = 0.1;
t = 0:dt:Tmax;
V0 = zeros(0, length(t));
V1 = zeros(0, length(t));
V2 = zeros(0, length(t));
V3 = zeros(0, length(t));
V4 = zeros(0, length(t));
Iapp = 0;
V0(1) = 0;
kv01 = (g * (V0(1) - E_L) + (Iapp / A) - (V0(1) / R_m)) / C_m;
av01 = V0(1) + kv01 * dt;
kv02 = (g * (av01 - E_L) + (Iapp / A) - (av01 / R_m)) / C_m;
V0(2) = V0(1) + (kv01 + kv02) * dt / 2;
V1(1) = V0(2);
kv11 = (g * (V1(1) - E_L) + (Iapp / A) - (V1(1) / R_m)) / C_m;
av11 = V1(1) + kv11 * dt;
kv12 = (g * (av11 - E_L) + (Iapp / A) - (av11 / R_m)) / C_m;
V1(2) = V1(1) + (kv11 + kv12) * dt / 2;
V2(1) = V1(2);
kv21 = (g * (V2(1) - E_L) + (Iapp / A) - (V2(1) / R_m)) / C_m;
av21 = V2(1) + kv21 * dt;
kv22 = (g * (av21 - E_L) + (Iapp / A) - (av21 / R_m)) / C_m;
V2(2) = V2(1) + (kv21 + kv22) * dt / 2;
V3(1) = V2(2);
kv31 = (g * (V3(1) - E_L) + (Iapp / A) - (V3(1) / R_m)) / C_m;
av31 = V3(1) + kv31 * dt;
kv32 = (g * (av31 - E_L) + (Iapp / A) - (av31 / R_m)) / C_m;
V3(2) = V3(1) + (kv31 + kv32) * dt / 2;
V4(1) = V3(2);
kv41 = (g * (V4(1) - E_L) + (Iapp / A) - (V4(1) / R_m)) / C_m;
av41 = V4(1) + kv41 * dt;
kv42 = (g * (av41 - E_L) + (Iapp / A) - (av41 / R_m)) / C_m;
V4(2) = V4(1) + (kv41 + kv42) * dt / 2;
for i=2:length(x)-1
    kv01 = (g * (V0(i) - E_L) + (Iapp / A) - (V0(i) / R_m)) / C_m;
    av01 = V0(i) + kv01 * dt;
    kv02 = (g * (av01 - E_L) + (Iapp / A) - (av01 / R_m)) / C_m;
    V0(i + 1) = V0(i) + (kv01 + kv02) * dt / 2;
    kv11 = (g * (V1(i) - E_L) + (Iapp / A) - (V1(i) / R_m)) / C_m;
    av11 = V1(i) + kv11 * dt;
    kv12 = (g * (av11 - E_L) + (Iapp / A) - (av11 / R_m)) / C_m;
    V1(i + 1) = V1(i) + (kv11 + kv12) * dt / 2;
    kv21 = (g * (V2(i) - E_L) + (Iapp / A) - (V2(i) / R_m)) / C_m;
    av21 = V2(i) + kv21 * dt;
    kv22 = (g * (av21 - E_L) + (Iapp / A) - (av21 / R_m)) / C_m;
    V2(i + 1) = V2(i) + (kv21 + kv22) * dt / 2;
    kv31 = (g * (V3(i) - E_L) + (Iapp / A) - (V3(i) / R_m)) / C_m;
    av31 = V3(i) + kv31 * dt;
    kv32 = (g * (av31 - E_L) + (Iapp / A) - (av31 / R_m)) / C_m;
    V3(i + 1) = V3(i) + (kv31 + kv32) * dt / 2;
    kv41 = (g * (V4(i) - E_L) + (Iapp / A) - (V4(i) / R_m)) / C_m;
    av41 = V4(i) + kv41 * dt;
    kv42 = (g * (av41 - E_L) + (Iapp / A) - (av41 / R_m)) / C_m;
    V4(i + 1) = V4(i) + (kv41 + kv42) * dt / 2;
end
fprintf("Question 2, Part d:\n");
fprintf("You can apply any amount of current to the first compartment to get a voltage reading in the last compartment.\n");

function SSInfinCableV = stdyState(lambda, r_L, a, I, rest, x)
    SSInfinCableV = (((lambda * r_L)/(pi * a^2)) * I * exp(-x/lambda)) - rest;
end