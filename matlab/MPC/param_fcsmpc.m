%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of machine's parameters and limits for PMSM_3ph_dq

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LUT from ANSYS Maxwell:
% id, iq, thetar --> fluxD, fluxQ, flux0, torque 
ee_ece_table;

% LUTs inversas de corrientes de Juan:
load('pmsm_phi2i_21p');

% Cambio para adaptar a Park original
data.i_d = fliplr(data.i_d);
data.i_q = fliplr(-data.i_q);
torque = fliplr(torque);
fluxD = fliplr(fluxD);
fluxQ = fliplr(-fluxQ);

% Machine's parameters 
param_R          = 0.07;               % Stator resistance
param_p          = N;                 % Pole pairs
param_flux_m     = 0.1571;              % Permanent magnet flux
param_Ld         = 0.0013;          % Inductance Ld
param_Lq         = 0.0039;          % Inductance Lq

% Load's parameters
param_J   = 2.9;             % Inertia coefficient
param_f   = 0;              % Friction coefficient

In = 500;   % Corriente de fase nominal 
Vn = 460;   % Tensión de fase nominal
wm_n = 314; % [rad/s] velocidad mecánica nominal
Tn = 3.9;   % par nominal


% DC supply
Vbus = 500;

% Control parameters
Kp_w  = 700;
Ki_w  = 15;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concordia Transformation
ang = 2*pi/3;

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
% Modelo para Park original
param_PMSM_A = [1/param_Ld, 0;
                0         , 1/param_Lq];
param_PMSM_B = -param_R*param_PMSM_A;
param_PMSM_C = -param_PMSM_A;
param_PMSM_D = [0               , -param_p*param_Lq;
                param_p*param_Ld,  0];
param_PMSM_E = [0;
                param_p*param_flux_m];

% % Modelo para Park modificada
% param_PMSM_A = [1/param_Ld, 0;
%                 0         , 1/param_Lq];
% param_PMSM_B = -param_R*param_PMSM_A;
% param_PMSM_C = -param_PMSM_A;
% param_PMSM_D = [0               , +param_p*param_Lq;
%                -param_p*param_Ld,  0];
% param_PMSM_E = [0;
%                -param_p*param_flux_m];


