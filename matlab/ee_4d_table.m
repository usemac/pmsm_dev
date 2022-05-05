%*******************************************************************************************
% Parámetros del modelo 4D
%*******************************************************************************************

% Ideal PMSM parameters
PM = 0.1;     % Permanent magnet flux
N = 4;        % Number of pole pairs
Ld = 0.0002;  % D-axis inductance
Lq = 0.0005;  % Q-axis inductance
L0 = 0.00018; % Zero-sequence inductance
Rs = 0.013;   % Stator resistance

% Test condition

% Load resistance [Ohm]
Rload = 1000000;
% Speed
RPM = 100;


%B_Sweepings
iA = linspace(-250,250,5);
iB = iA;
iC = iA;
X = pi/180*linspace(0,360/N,180/N+1);
%E_Sweepings

%B_OutputMatrix F
[F,T,dFdA,dFdB,dFdC,dFdX] = ee_generateIdealPMSMfluxData(PM,Ld,Lq,L0,iA,iB,iC,X);
