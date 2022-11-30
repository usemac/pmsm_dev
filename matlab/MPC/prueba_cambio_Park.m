% INTENTO DE DEFINIR LUTs PARA LA PARK ORIGINAL
close all;
clear all;
clc;

%%
syms theta

P1 = 2/3*[ cos(theta)    cos(theta-2*pi/3)   cos(theta-4*pi/3);
          -sin(theta)   -sin(theta-2*pi/3)  -sin(theta-4*pi/3)];

invP1 = [ cos(theta)            -sin(theta);
          cos(theta-2*pi/3)     -sin(theta-2*pi/3);
          cos(theta-4*pi/3)     -sin(theta-4*pi/3)];

P2 = 2/3*[ cos(theta)    cos(theta-2*pi/3)   cos(theta-4*pi/3);
           sin(theta)    sin(theta-2*pi/3)   sin(theta-4*pi/3)];

invP2 = [ cos(theta)             sin(theta);
          cos(theta-2*pi/3)      sin(theta-2*pi/3);
          cos(theta-4*pi/3)      sin(theta-4*pi/3)];

% Si tenemos magnitudes en 2 y queremos para a 1:
% iabc = invP2*idq2 --> idq1 = P1*invP2*idq2

T = P1*invP2
simplify(T)