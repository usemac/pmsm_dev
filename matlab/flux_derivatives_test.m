%% Prueba de la función de cálculo de las derivadas parciales del flujo
clear all;
clc;

%% Creación de una matrix de flujos, en la fase a, para una PMSM ideal
% Parámetros de la máquinca
PM = 0.1;     % Permanent magnet flux
N = 6;        % Number of pole pairs
Ld = 0.0002;  % D-axis inductance
Lq = 0.0002;  % Q-axis inductance
L0 = 0.00018; % Zero-sequence inductance
Rs = 0.013;   % Stator resistance

% Definición de las variables de las que dependerá el flujo
iA = linspace(-250,250,5);
iB = iA;
iC = iA;
X = pi/180*linspace(0,360/N,180/N+1);

% Flujo de la fase A como función de ia, ib, ic y Theta_r
F1 = ee_generateIdealPMSMfluxData(PM,Ld,Lq,L0,iA,iB,iC,X);

%% Flujo y sus derivadas
[F,T,dFdA1,dFdB1,dFdC1,dFdX1] = ee_generateIdealPMSMfluxData(PM,Ld,Lq,L0,iA,iB,iC,X);

%% Derivadas del flujo
[dFdA,dFdB,dFdC,dFdX, iD, iQ] = ee_calculateFluxPartialDerivatives(iA,iB,iC,X,F);

%% Derivadas para un punto concreto
[dFdAp,dFdBp,dFdCp,dFdXp] = ee_calculateFluxPartialDerivatives(iA(1),iB(1),iC(1),X(:),F(1,1,1,:));

% NECESITO OBTENER LAS TABLAS DE FLUJOS COMO FUNCIÓN DE IA, IB, IC Y THETAR
