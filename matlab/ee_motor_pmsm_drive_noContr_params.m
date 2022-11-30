clear all;
clc;

%% Modelo con LUTs inversas

% LUT from ANSYS Maxwell:
% id, iq, thetar --> fluxD, fluxQ, flux0, torque 
ee_ece_table;

PM = 0.1571;  % Permanent magnet flux
Ld = 0.0013;  % D-axis inductance
Lq = 0.0039;  % Q-axis inductance
L0 = 0.0013;  % Zero-sequence inductance
Rs = 0.07;    % Stator resistance
J = 0.02;


% LUTs inversas de corrientes de Juan
load('pmsm_phi2i_21p');

%% Controller parameters

tsim = 1;
Ts = 2e-6;                % Fundamental sample time
fsw = 2e3;                % Switching frequency (Hz)
fc = fsw*10;              % Control loop frequency (Hz)
Tsc = 1/fc;               % Control loop sample time
Imax = 1000;               % Assumed max stator current (peak value)
Tmax = 1.5*N*PM*Imax;     % Maximum electromagnetic torque
fnom = 480;               % Nominal frequency (Hz)
rpm_nom = 60*fnom/N;      % Nominal rotor speed in rpm
omegam_nom = 2*pi*fnom/N; % Nominal mechanical rotor speed (rad/s)
% Pmax = omegam_nom*Tmax;   % Maximum power
rpm0 = 7200;
f0 = rpm0*N/60;
torque0 = 0;

Vdc = 50;
% vpor = 265;
% 
% % Control de velocidad
% kpn       =    5000;
% kin       =    5000;
% 
% % Control de id
% kpf       =    5;
% kif       =    10;
% 
% % Control de iq
% kpt       =    2;
% kit       =    1;