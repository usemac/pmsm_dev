clear all
close all
clc

%%
param_fcsmpc;

%% Prueba de interpolación trilinear
iD = data.i_d;
iQ = data.i_q;
fD = data.phi_d;
fQ = data.phi_d;
TH = data.theta;

% Elegimos un punto cualquiera
% x-->fd, y-->fq, z-->theta
x = 0.0356;
y = -0.09661;
z = 24.13;

value_D = trilinear_interp(x,y,z,fD,fQ,TH,iD)
value_Q = trilinear_interp(x,y,z,fD,fQ,TH,iQ)

% Valores por interporlación Matlab
value_D_interpn = interpn(data.phi_d,data.phi_q',data.theta,data.i_d,x,y,z,'linear')
value_Q_interpn = interpn(data.phi_d,data.phi_q',data.theta,data.i_q,x,y,z,'linear')

%% Segunda Prueba 
x = 52.3;
y = -10.15;
z = 5.134;

fD = trilinear_interp(x,y,z,idVec,iqVec,angleVec,fluxD);
fQ = trilinear_interp(x,y,z,idVec,iqVec,angleVec,fluxQ);

fD_matlab = interpn(idVec,iqVec',angleVec,fluxD,x,y,z,'linear');
fQ_matlab = interpn(idVec,iqVec',angleVec,fluxQ,x,y,z,'linear');

% for x=1:nx
%     % iD
%     fig1 = figure;
%     set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
%     hold on
%     surf(fqVec,fdVec,iD(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\iti_d}','FontSize',tamletra, 'FontName', font)
%     view([45 25])
% %     zlim([0.05 0.25])
% %     xlim([0 800])
%     grid
%     colormap('parula')
%     set(fig1,'Position',[20 50 360 270])
% 
%     % iQ
%     fig1 = figure;
%     set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
%     hold on
%     surf(fqVec,fdVec,iQ(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
%     view([215 25])
% %     zlim([0.05 0.25])
% %     xlim([0 800])
%     grid
%     colormap('parula')
%     set(fig1,'Position',[20 50 360 270])
% end
