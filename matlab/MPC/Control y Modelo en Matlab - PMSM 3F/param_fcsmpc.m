%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of machine's parameters and limits for PMSM_3ph_dq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global param_Ld1 param_flux_max_1 param_Lq1 param_Ld3 param_flux_max_3 param_Lq3 param_phi_3
% global param_R param_p Vmax_ph_ph Imax_ph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Machine's parameters (sacados de
% https://imperix.com/doc/implementation/field-oriented-control-of-pmsm)
param_R          = 0.07;               % Stator resistance
param_p          = 4;                 % Pole pairs
param_flux_m     = 0.1571;              % Permanent magnet flux
param_Ld         = 0.0013;          % Inductance Ld
param_Lq         = 0.0039;          % Inductance Lq
% Load's parameters
param_J   = 2.9;             % Inertia coefficient
param_f   = 0.0;              % Friction coefficient

In = 100;   % Corriente de fase nominal 
Vn = 460;   % Tensión de fase nominal
wm_n = 314; % [rad/s] velocidad mecánica nominal
Tn = 3.9;   % par nominal

% DC supply
Vbus = 500;

% Control parameters
Kp_w  = 1000;
Ki_w  = 15;

% % Machine's parameters (curso de Matlab)
% param_R          = 2*5.8264e-3;               % Stator resistance
% param_p          = 26;                 % Pole pairs
% param_flux_m     = 5.8264*sqrt(2);              % Permanent magnet flux
% param_Ld         = 2*1.5731e-3;          % Inductance Ld
% param_Lq         = 2*1.5731e-3;          % Inductance Lq
% % Load's parameters
% param_J   = 1e4;             % Inertia coefficient
% param_f   = 0.0;              % Friction coefficient
% 
% In = 1867.76;   % Corriente nominal
% 
% % DC supply
% Vbus = 1300;

% % Limits
% Imax_ph_VSI = 50;              % Maximum phase peak current (VSI)
% % Imax_ph_TH  = 64;               % Maximum phase RMS current (thermal)
% Imax_ph     = Imax_ph_VSI;      % Current limit
% Vmax_ph_VSI = Vbus/2;           % Maximum phase peak voltage in the VSI
% Vmax_ph_ph  = Vbus;             % Maximum phase-to-phase peak voltage

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concordia Transformation
ang = 2*pi/3;
param_M = 2/3*[     1        cos(ang)   cos(2*ang);
                    0        sin(ang)   sin(2*ang);
                  1/2            1/2          1/2];

% Star Coupling Transformation
param_T = (1/3)*[ 2 -1 -1;
                 -1  2 -1;
                 -1 -1  2]; 

% número de vectores de voltaje diferentes
param_NV = 8; 
% Posibles matrices de vectores de disparo considerados por el optimizador
XI8 = [ ...
     0 0 0;     % 0
     0 0 1;     % 1
     0 1 0;     % 2
     0 1 1;     % 3
     1 0 0;     % 4
     1 0 1;     % 5
     1 1 0;     % 6
     1 1 1];    % 7

% las ecs de la PMSM se escriben enel espacio de estados
% de forma que x=(id, iq)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% param_PMSM_A = [1/param_Ld, 0;
%                 0         , 1/param_Lq];
% param_PMSM_B = -param_R*param_PMSM_A;
% param_PMSM_C = -param_PMSM_A;
% param_PMSM_D = [0               , -param_p*param_Lq;
%                 param_p*param_Ld,  0];
% param_PMSM_E = [0;
%                 param_p*param_flux_m];

param_PMSM_A = [1/param_Ld, 0;
                0         , 1/param_Lq];
param_PMSM_B = -param_R*param_PMSM_A;
param_PMSM_C = -param_PMSM_A;
param_PMSM_D = [0               , +param_p*param_Lq;
                -param_p*param_Ld,  0];
param_PMSM_E = [0;
                -param_p*param_flux_m];

