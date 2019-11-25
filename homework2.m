% Matthew Kobilas
% Professor Horacio Rotstein
% MATH 430-001
% 25 September 2019

% E = (RT/zF) * ln(ion_out/ion_in)

R = 8.314;
%R = 1.98;
K = 273.15;
F = 96480;

% log function is base e
% returns resting membrane potential in Volts
nernst=@(T, z, concOut, concIn) ((R*(T+K))/(z*F)) * log(concOut/concIn);

% K, K_out = 5 mmol/L, K_in = 150 mmol/L, valence = 1
fprintf("K equilibrium at 20 degrees Celsius: %1.4f mV\n", ...
    nernst(20, 1, 5, 150)*1000);
fprintf("K equilibrium at 25 degrees Celsius: %1.4f mV\n", ...
    nernst(25, 1, 5, 150)*1000);

% Na, Na_out = 150 mmol/L, Na_in = 15 mmol/L, valence = 1
fprintf("Na equilibrium at 20 degrees Celsius: %1.4f mV\n", ...
    nernst(20, 1, 150, 15)*1000);
fprintf("Na equilibrium at 25 degrees Celsius: %1.4f mV\n", ...
    nernst(25, 1, 150, 15)*1000);

% Cl, Cl_out = 125 mmol/L, Cl_in = 10 mmol/L, valence = -1
fprintf("Cl equilibrium at 20 degrees Celsius: %1.4f mV\n", ...
    nernst(20, -1, 125, 10)*1000);
fprintf("Cl equilibrium at 25 degrees Celsius: %1.4f mV\n", ...
    nernst(25, -1, 125, 10)*1000);

% Ca, Ca_out = 2 mmol/L, Ca_in = 0.0002 mmol/L, valence = 2
fprintf("Ca equilibrium at 20 degrees Celsius: %1.4f mV\n", ...
    nernst(20, 2, 2, 0.0002)*1000);
fprintf("Ca equilibrium at 25 degrees Celsius: %1.4f mV\n", ...
    nernst(25, 2, 2, 0.0002)*1000);

% returns concentration relation of an ion using nernst equation
oppNernst=@(T, z, E) exp((E*z*F)/(R*(T+K)));

% X, E_X = -60 mV, valence = 1
fprintf("X concentration relation at 20 degrees Celsius: %1.4f required to maintain potential at -60 mV\n", ...
    oppNernst(20, 1, -0.06));
fprintf("X concentration relation at 25 degrees Celsius: %1.4f required to maintain potential at -60 mV\n", ...
    oppNernst(25, 1, -0.06));

% tau*dV/dt = -(V - E_L) + R*I_app
% tau = RC
% V(0) = E_L, R = 100 MOhm, C = 100 pFarad, I_app = 0.25 nAmp, E_L = -60mV
% C: Capacitance (Farad or Coulomb/Volt)
% R: Resistance (Ohm or Volt/Ampere)
% I: Current (Ampere or Coulomb/second)
% Q: Charge (Coulomb)
% V: Potential difference (Volt or Joule/Coulomb)

close all;

% Biophysical parameters

C = 0.1;
El = -60;
Gl = 0.01;
Iapp = 0.25;

% tau=R*C=(V/A)*(C/V)=((J/C)/(C/s))*(C/(J/C))=(Js/CC)*(CC/J) = seconds
R = 1/Gl;
tau = R*C; % seconds
fprintf("Tau=RC is equal to %1.4f milliseconds\n", tau);
%RI_app=(V/A)*(C/s)=((J/C)/(C/s))*(C/s)=(Js/CC)(C/s)=J/C = Volts
RI_app = R*Iapp;
fprintf("R*I_app is equal to %1.4f milliVolts\n", RI_app);

timeTaken=@(V_target)-tau * log(1-(V_target-El)/RI_app);
if(1-(-50-El)/RI_app) > 0
    fprintf("Time until V = -50mV: %1.4f milliseconds\n", timeTaken(-50));
else
    fprintf("Time until V = -50mV: UNDEFINED\n");
end
if(1-(-30-El)/RI_app) > 0
    fprintf("Time until V = -30mV: %1.4f milliseconds\n", timeTaken(-30));
else
    fprintf("Time until V = -30mV: UNDEFINED\n");
end

% Time definitions
Tmax = 1000;
dt = 0.1;
t = 0:dt:Tmax;

% Square wave (Heaviside function)
ti = dt;
tf = 10000;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;

% Initial conditions (for D = 0)
V = zeros(1,length(t));
V(1) = El;

% Computation of the solution (for D = 0)
for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;
end

% Graph
figure
hold on
plot(t,V,'b','linewidth',2);
axis([0 Tmax -80 80]);
set(gca,'fontsize',24);
xlabel('t (ms)');
ylabel('V (mV)');
title("Passive Membrane Model");