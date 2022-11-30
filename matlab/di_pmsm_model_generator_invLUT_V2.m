clear all;
close all;
clc;

%% Parameters initialization

% LUT from ANSYS Maxwell:
% id, iq, thetar --> fluxD, fluxQ, flux0, torque 
ee_ece_table;

% Stator resistance [Ohm]
R_s = 0.07;

% Load resistance [Ohm]
% Rload = 1000000;
Rload = 10;

% Speed
RPM = 7200;
J = 0.02;

% %% PRUEBA: flujo dq para corrientes 0:
% % id = 0
% id = idVec(11);
% iq = iqVec(11);
% 
% for x=1:length(angleVec)
%     fd(x) = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(x));
%     fq(x) = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(x));
%     f0(x) = interpn(idVec,iqVec',angleVec,flux0,id,iq,angleVec(x));
% end
% 
% figure
% subplot(3,1,1)
% plot(angleVec,fd)
% legend('\phi_d')
% xlabel('\theta [º]')
% ylabel('\phi [Wb]')
% subplot(3,1,2)
% plot(angleVec,fq)
% legend('\phi_q')
% xlabel('\theta [º]')
% ylabel('\phi [Wb]')
% subplot(3,1,3)
% plot(angleVec,f0)
% legend('\phi_0')
% xlabel('\theta [º]')
% ylabel('\phi [Wb]')


% %% Inversión de las LUT de los flujos D y Q (MÉTODO 1)
% close all;
% 
% % Búsqueda de los rangos de valores de los flujos D y Q en las LUTs
% % Máximos de cada LUT
% fluxD_max = max(fluxD,[],'all');
% fluxQ_max = max(fluxQ,[],'all');
% % Mínimos de cada LUT
% fluxD_min = min(fluxD,[],'all');
% fluxQ_min = min(fluxQ,[],'all');
% % En vista de los resultados, definimos los siguientes valores
% fdVec = [-0.2:0.05:0.4];
% fqVec = [-0.4:0.05:0.4];
% % fdVec = [-0.1:0.05:0.3];
% % fqVec = [-0.3:0.05:0.3];
% % % Punto real de las LUTs, id=0 e iq=90 (FUNCIONA)
% % fdVec = fluxD(11,14,1)
% % fqVec = fluxQ(11,14,1)
% 
% % Reducimos el vector de ángulos
% % xVec = 0;
% % xVec = [0:5:90];
% 
% % Inversión mediante minimización
% nx  = length(angleVec);
% nfd = length(fdVec);
% nfq = length(fqVec);
% 
% % idVec_aux = [-350:5:350];
% % iqVec_aux = [-350:5:350];
% idVec_aux = idVec;
% iqVec_aux = iqVec;
% 
% nid = length(idVec_aux);
% niq = length(iqVec_aux);
% 
% iD = zeros(nfd,nfq,nx);
% iQ = zeros(nfd,nfq,nx);
% 
% metodo = 'linear';
% 
% % ESTA INVERSIÓN DA BASTANTES ERRORES
% for x=1:nx
%     for fd=1:nfd
%         fluxd0 = fdVec(fd);
%         for fq=1:nfq
%             fluxq0 = fqVec(fq);
%             % Buscamos los valores id e iq óptimos
%             flux_err_opt = 1000000;
%             for indexd=1:nid
%                 id = idVec_aux(indexd);
%                 for indexq=1:niq
%                     iq = iqVec_aux(indexq);
% %                     if id < idVec(1) || id > idVec(length(idVec)) || iq < iqVec(1) || iq > iqVec(length(iqVec))
% %                         metodo = 'makima';
% %                     else
% %                         metodo = 'linear';
% %                     end
%                     fluxd = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(x),metodo);
%                     fluxq = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(x),metodo);
%                     flux_err = abs(fluxd - fluxd0) + abs(fluxq - fluxq0);
%                     if flux_err < flux_err_opt
%                         flux_err_opt = flux_err;
%                         % Gruardamos el valor óptimo obtenido
%                         iD(fd,fq,x) = id;
%                         iQ(fd,fq,x) = iq;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % Representación de las corrientes para un valor de thetar
% % Parámetros de representación
% tamletra = 10;
% tamlinea = 1;
% font = 'Times';
% tammarker = 4;
% 
% nx  = length(angleVec);
% nfd = length(fdVec);
% nfq = length(fqVec);
% 
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
% %%
% % Comprobación de las nuevas LUTs
% metodo = 'linear';
% 
% % Comparando corrientes
% nid = length(idVec);
% niq = length(iqVec);
% id_err = zeros(nid,niq,nx);
% iq_err = zeros(nid,niq,nx);
% 
% for x=1:nx
%     for indexd=1:nid
%         id = idVec_aux(indexd);
%         for indexq=1:niq
%             iq = iqVec_aux(indexq);
%             fluxd = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(x),metodo);
%             fluxq = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(x),metodo);
%             id_LUT = interpn(fdVec,fqVec',angleVec,iD,fluxd,fluxq,angleVec(x),metodo);
%             iq_LUT = interpn(fdVec,fqVec',angleVec,iQ,fluxd,fluxq,angleVec(x),metodo);
%             id_err(indexd,indexq,x) = id-id_LUT;
%             iq_err(indexd,indexq,x) = iq-iq_LUT;
%             
%         end
%     end
%     % iD_err
%     fig1 = figure;
%     set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
%     hold on
%     surf(iqVec,idVec,id_err(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\iti_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\iti_d error}','FontSize',tamletra, 'FontName', font)
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
%     surf(iqVec,idVec,iq_err(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\iti_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\iti_q error}','FontSize',tamletra, 'FontName', font)
%     view([215 25])
% %     zlim([0.05 0.25])
% %     xlim([0 800])
%     grid
%     colormap('parula')
%     set(fig1,'Position',[20 50 360 270])
% end
% 
% % Comparando flujos
% fd_err = zeros(nfd,nfq,nx);
% fq_err = zeros(nfd,nfq,nx);
% 
% for x=1:nx
%     for indexd=1:nfd
%         fluxd = fdVec(indexd);
%         for indexq=1:nfq
%             fluxq = fqVec(indexq);
%             id = interpn(fdVec,fqVec',angleVec,iD,fluxd,fluxq,angleVec(x),metodo);
%             iq = interpn(fdVec,fqVec',angleVec,iQ,fluxd,fluxq,angleVec(x),metodo);
%             fluxd_LUT = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(x),metodo);
%             fluxq_LUT = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(x),metodo);
%             fd_err(indexd,indexq,x) = fluxd-fluxd_LUT;
%             fq_err(indexd,indexq,x) = fluxq-fluxq_LUT;
%         end
%     end
%     % fD_err
%     fig1 = figure;
%     set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
%     hold on
%     surf(fqVec,fdVec,fd_err(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\it\phi_d error}','FontSize',tamletra, 'FontName', font)
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
%     surf(fqVec,fdVec,fq_err(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\it\phi_q error}','FontSize',tamletra, 'FontName', font)
%     view([215 25])
% %     zlim([0.05 0.25])
% %     xlim([0 800])
%     grid
%     colormap('parula')
%     set(fig1,'Position',[20 50 360 270])
% end
% 
% %% Inversión de las LUT de los flujos D y Q (MÉTODO 2)
% close all;
% 
% % Creación de la matriz de datos de entrada a la función
% idVec_aux = [-350:5:350];
% iqVec_aux = [-350:5:350];
% 
% nx  = length(angleVec);
% nid = length(idVec_aux);
% niq = length(iqVec_aux);
% 
% % Tendremos una matriz de nid*niq filas, 4 columnas (id,iq,fluxd,fluxq), y
% % nx páginas (1 por ángulo)
% InputPoints = zeros(nid*niq,4,nx);
% 
% metodo = 'makima';
% 
% for x=1:nx
%     cont = 1;
%     thetar = angleVec(x);
%     for index_d=1:nid
%         id = idVec_aux(index_d);
%         for index_q=1:niq
%             iq = iqVec_aux(index_q);
%             % Buscamos los valores de flujos
%             fluxd = interpn(idVec,iqVec',angleVec,fluxD,id,iq,thetar,metodo);
%             fluxq = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,thetar,metodo);
%             % Rellenados la matriz de datos
%             InputPoints(cont,:,x) = [id,iq,fluxd,fluxq];
%             cont = cont+1;
%         end
%     end
% end
% 
% % Búsqueda de los rangos de valores de los flujos D y Q en las LUTs
% % Máximos de cada LUT
% fluxD_max = max(fluxD,[],'all');
% fluxQ_max = max(fluxQ,[],'all');
% % Mínimos de cada LUT
% fluxD_min = min(fluxD,[],'all');
% fluxQ_min = min(fluxQ,[],'all');
% % En vista de los resultados, definimos los siguientes valores
% fdVec = [-0.2:0.05:0.4];
% fqVec = [-0.4:0.05:0.4];
% % fdVec = [round(fluxD_min,3):0.01:round(fluxD_max,3)];
% % fqVec = [round(fluxQ_min,3):0.01:round(fluxQ_max,3)];
% % Inicialización de las LUT
% nfd = length(fdVec);
% nfq = length(fqVec);
% iD = zeros(nfd,nfq,nx);
% iQ = zeros(nfd,nfq,nx);
% 
% Smoothness = 0.005;
% 
% % Búsqueda de las LUT
% for i = 1:nx
%     z1 = RegularizeData3D(InputPoints(:,3,i), InputPoints(:,4,i), InputPoints(:,1,i), fdVec, fqVec, 'interp', 'bicubic', 'smoothness', Smoothness);
%     iD(:,:,i) = z1';
%     z2 = RegularizeData3D(InputPoints(:,3,i), InputPoints(:,4,i), InputPoints(:,2,i), fdVec, fqVec, 'interp', 'bicubic', 'smoothness', Smoothness);
%     iQ(:,:,i) = z2';
% end

% % Representación de las corrientes para un valor de thetar
% % Parámetros de representación
% tamletra = 10;
% tamlinea = 1;
% font = 'Times';
% tammarker = 4;
% 
% x = 1;
% % for x=1:7:nx
%     % iD
%     fig1 = figure;
%     set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
%     hold on
%     surf(fqVec,fdVec,iD(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',0.5)
%     scatter3(InputPoints(:, 4), InputPoints(:, 3), InputPoints(:, 1), 'fill');
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
%     surf(fqVec,fdVec,iQ(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',0.5)
%     scatter3(InputPoints(:, 4), InputPoints(:, 3), InputPoints(:, 2), 'fill');
%     ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
%     view([215 25])
% %     zlim([0.05 0.25])
% %     xlim([0 800])
%     grid
%     colormap('parula')
%     set(fig1,'Position',[20 50 360 270])
% % end
% 
% % Comprobación de las nuevas LUTs
% metodo = 'linear';
% 
% % Comparando corrientes
% id_err = zeros(nid,niq,nx);
% iq_err = zeros(nid,niq,nx);
% 
% for x=1:7:nx
%     for indexd=1:nid
%         id = idVec(indexd);
%         for indexq=1:niq
%             iq = iqVec(indexq);
%             fluxd = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(x),metodo);
%             fluxq = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(x),metodo);
%             id_LUT = interpn(fdVec,fqVec',angleVec,iD,fluxd,fluxq,angleVec(x),metodo);
%             iq_LUT = interpn(fdVec,fqVec',angleVec,iQ,fluxd,fluxq,angleVec(x),metodo);
%             id_err(indexd,indexq,x) = id-id_LUT;
%             iq_err(indexd,indexq,x) = iq-iq_LUT;
%         end
%     end
%     % iD_err
%     fig1 = figure;
%     set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
%     hold on
%     surf(iqVec,idVec,id_err(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\iti_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\iti_d error}','FontSize',tamletra, 'FontName', font)
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
%     surf(iqVec,idVec,iq_err(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\iti_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\iti_q error}','FontSize',tamletra, 'FontName', font)
%     view([215 25])
% %     zlim([0.05 0.25])
% %     xlim([0 800])
%     grid
%     colormap('parula')
%     set(fig1,'Position',[20 50 360 270])
% end
% 
% % Comparando flujos
% fd_err = zeros(nfd,nfq,nx);
% fq_err = zeros(nfd,nfq,nx);
% 
% for x=1:7:nx
%     for indexd=1:nfd
%         fluxd = fdVec(indexd);
%         for indexq=1:nfq
%             fluxq = fqVec(indexq);
%             id = interpn(fdVec,fqVec',angleVec,iD,fluxd,fluxq,angleVec(x),metodo);
%             iq = interpn(fdVec,fqVec',angleVec,iQ,fluxd,fluxq,angleVec(x),metodo);
%             fluxd_LUT = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(x),metodo);
%             fluxq_LUT = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(x),metodo);
%             fd_err(indexd,indexq,x) = fluxd-fluxd_LUT;
%             fq_err(indexd,indexq,x) = fluxq-fluxq_LUT;
%         end
%     end
%     % fD_err
%     fig1 = figure;
%     set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
%     hold on
%     surf(fqVec,fdVec,fd_err(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\it\phi_d error}','FontSize',tamletra, 'FontName', font)
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
%     surf(fqVec,fdVec,fq_err(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
%     ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\it\phi_q error}','FontSize',tamletra, 'FontName', font)
%     view([215 25])
% %     zlim([0.05 0.25])
% %     xlim([0 800])
%     grid
%     colormap('parula')
%     set(fig1,'Position',[20 50 360 270])
% end
% 
% %% Comprobación adicional de los flujos
% id = idVec(11);
% iq = iqVec(14);
% 
% metodo = 'makima';
% 
% fluxd = zeros(1,nx);
% fluxq = zeros(1,nx);
% fluxd_new = zeros(1,nx);
% fluxq_new = zeros(1,nx);
% 
% for i=1:nx
%     fluxd(i) = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(i),metodo);
%     fluxq(i) = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(i),metodo);
%     currentd(i) = interpn(fdVec,fqVec',angleVec,iD,fluxd(i),fluxq(i),angleVec(i),metodo);
%     currentq(i) = interpn(fdVec,fqVec',angleVec,iQ,fluxd(i),fluxq(i),angleVec(i),metodo);
%     fluxd_new(i) = interpn(idVec,iqVec',angleVec,fluxD,currentd(i),currentq(i),angleVec(i),metodo);
%     fluxq_new(i) = interpn(idVec,iqVec',angleVec,fluxQ,currentd(i),currentq(i),angleVec(i),metodo);
% end
% 
% figure
% hold on
% plot(angleVec,fluxd,angleVec,fluxq)
% plot(angleVec,fluxd_new,'--',angleVec,fluxq_new,'--')
% grid on

% %% Obtention of new LUT: fluxA, fluxB, and fluxC as functions of iA, iB, iC, and thetar
% 
% % Extension of the fluxD, fluxQ, flux0, and torque LUTs to a maximum thetar of 360/N
% % 30 (360/3/N) --> 90 (360/N)
% for i=1:30
%     fluxD(:,:,31+i) = fluxD(:,:,i+1);
%     fluxQ(:,:,31+i) = fluxQ(:,:,i+1);
%     flux0(:,:,31+i) = flux0(:,:,i+1);
%     torque(:,:,31+i) = torque(:,:,i+1);
% end
% for i=1:30
%     fluxD(:,:,61+i) = fluxD(:,:,i+1);
%     fluxQ(:,:,61+i) = fluxQ(:,:,i+1);
%     flux0(:,:,61+i) = flux0(:,:,i+1);
%     torque(:,:,61+i) = torque(:,:,i+1);
% end
% 
% angleVec = [0:1:90];
% niD = length(idVec);
% niQ = length(iqVec);
% nang = length(angleVec);
% 
% % Calulation of phase fluxes LUTs as funtions of id, iq, and thetar
% fluxA = zeros(niD,niQ,nang);
% fluxB = zeros(niD,niQ,nang);
% fluxC = zeros(niD,niQ,nang);
% for x=1:nang
%     for d=1:niD
%         for q=1:niQ
%             % Clarke-Park matrix
%             theta = N*angleVec(x)*pi/180;
%             T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
%                      sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
%                      1/2           1/2                    1/2];
%             invT = inv(T);
% 
%             % fluxD, fluxQ, flux0, thetar --> fluxA, fluxB, fluxC
%             fluxfase = invT*[fluxD(d,q,x);fluxQ(d,q,x);flux0(d,q,x)];
%             fluxA(d,q,x) = fluxfase(1);
%             fluxB(d,q,x) = fluxfase(2);
%             fluxC(d,q,x) = fluxfase(3);
% 
% %             % Phase currents
% %             ifase = invT(:,1:2)*[idVec(d);iqVec(q)];
% %             iA(d,q,x) = ifase(1);
% %             iB(d,q,x) = ifase(2);
% %             iC(d,q,x) = ifase(3);
%         end
%     end
% end
% 
% % ifa(1,1:91)=iA(11,15,:);
% % ifb(1,1:91)=iB(11,15,:);
% % ifc(1,1:91)=iC(11,15,:);
% % figure
% % plot(angleVec,ifa,angleVec,ifb,angleVec,ifc)
% % grid on
% % 
% % fd(1,1:91) = fluxD(11,15,:);
% % fq(1,1:91) = fluxQ(11,15,:);
% % figure
% % plot(angleVec,fd,angleVec,fq)
% % legend('\phi_{d}','\phi_{q}')
% % grid on
% % 
% % fa(1,1:91) = fluxA(11,15,:);
% % fb(1,1:91) = fluxB(11,15,:);
% % fc(1,1:91) = fluxC(11,15,:);
% % figure
% % plot(angleVec,fa,angleVec,fb,angleVec,fc)
% % legend('\phi_{a}','\phi_{b}','\phi_{c}')
% % grid on
% 
% % Selection of the phase currents values for the new LUTs
% % iaVec = [-300 -270 -240 -210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210 240 270 300];
% % ibVec = [-300 -270 -240 -210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210 240 270 300];
% % icVec = [-300 -270 -240 -210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210 240 270 300];
% iaVec = [-210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210];
% ibVec = [-210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210];
% icVec = [-210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210];
% 
% niA = length(iaVec);
% niB = length(ibVec);
% niC = length(icVec);
% 
% fluxA2 = zeros(niA,niB,niC,nang);
% fluxB2 = zeros(niA,niB,niC,nang);
% fluxC2 = zeros(niA,niB,niC,nang);
% 
% method = 'linear';
% 
% % Calculation of LTUs for the phase fluxes in terms of iA, iB, iC and
% % thetar
% for x=1:nang
%     for a=1:niA
%         for b=1:niB
%             for c=1:niC
%                 % Clarke-Park matrix
%                 theta = N*angleVec(x)*pi/180;
%                 T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
%                          sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
%                          1/2           1/2                    1/2];
% %                 T = 2/3*[ cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
% %                          -sin(theta)   -sin(theta - 2*pi/3)   -sin(theta + 2*pi/3)
% %                           1/2           1/2                    1/2];
%                 % d-q currents
%                 idq = T*[iaVec(a);ibVec(b);icVec(c)];
%                 % Interpolation form the phase fluxes LTUs obtained before
%                 fluxA2(a,b,c,x) = interpn(idVec,iqVec',angleVec,fluxA,idq(1),idq(2),angleVec(x),method);
%                 fluxB2(a,b,c,x) = interpn(idVec,iqVec',angleVec,fluxB,idq(1),idq(2),angleVec(x),method);
%                 fluxC2(a,b,c,x) = interpn(idVec,iqVec',angleVec,fluxC,idq(1),idq(2),angleVec(x),method);
%             end
%         end
%     end
% end
% 
% % idq = [idVec(10);iqVec(15)];
% % for i = 1:length(angleVec)
% %     theta = N*angleVec(i)*pi/180;
% %     T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
% %              sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
% %              1/2           1/2                    1/2];
% %     invT = inv(T);
% %     ifase(:,i) = invT(:,1:2)*idq;
% % 
% %     fa2(1,i) = interpn(iaVec,ibVec,icVec,angleVec,fluxA2,ifase(1),ifase(2),ifase(3),angleVec(i));
% %     fb2(1,i) = interpn(iaVec,ibVec,icVec,angleVec,fluxB2,ifase(1),ifase(2),ifase(3),angleVec(i));
% %     fc2(1,i) = interpn(iaVec,ibVec,icVec,angleVec,fluxC2,ifase(1),ifase(2),ifase(3),angleVec(i));
% % end
% 
% 
% % figure
% % plot(angleVec,fa,angleVec,fb,angleVec,fc)
% % hold on
% % plot(angleVec,fa2,':',angleVec,fb2,':',angleVec,fc2,':')
% % legend('\phi_{a2}','\phi_{b2}','\phi_{c2}')
% % grid on
% 
% %% Comprobamos que los flujos son iguales
% id = idVec(11);
% iq = iqVec(15);
% 
% for i=1:length(angleVec)
% %     % Flujos extraídos directamente de la LUT sin interpolar
% %     fd(i) = fluxD(11,15,i);
% %     fq(i) = fluxQ(11,15,i);
% %     f0(i) = flux0(11,15,i);
%     % Flujos de fase sacados de fluxD y fluxQ
%     flujod(i) = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(i));
%     flujoq(i) = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(i));
%     flujo0(i) = interpn(idVec,iqVec',angleVec,flux0,id,iq,angleVec(i));
%     % Pasamos a flujo abc
%     theta = N*angleVec(i)*pi/180;
%     T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
%              sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
%              1/2           1/2                    1/2];
% %     T = 2/3*[ cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
% %              -sin(theta)   -sin(theta - 2*pi/3)   -sin(theta + 2*pi/3)
% %               1/2           1/2                    1/2];
%     invT = inv(T);
%     
% %     fabc = invT*[fd(i);fq(i);f0(i)];
% %     fa(i) = fabc(1); 
% %     fb(i) = fabc(2); 
% %     fc(i) = fabc(3);
% 
%     flujosabc = invT*[flujod(i);flujoq(i);flujo0(i)];
%     flujoa(i) = flujosabc(1); 
%     flujob(i) = flujosabc(2); 
%     flujoc(i) = flujosabc(3);
%     
%     % Flujos de fase sacados de fluxA y fluxB y fluxC
% %     flujoa1_1(i) = fluxA(11,15,i);
% %     flujob1_1(i) = fluxB(11,15,i);
% %     flujoc1_1(i) = fluxC(11,15,i);
%     flujoa1(i) = interpn(idVec,iqVec',angleVec,fluxA,id,iq,angleVec(i));
%     flujob1(i) = interpn(idVec,iqVec',angleVec,fluxB,id,iq,angleVec(i));
%     flujoc1(i) = interpn(idVec,iqVec',angleVec,fluxC,id,iq,angleVec(i));
% 
%     % Flujos de fase sacados de fluxA2 y fluxB2 y fluxC2
%     iabc = invT(:,1:2)*[id;iq];
%     flujoa2(i) = interpn(iaVec,ibVec,icVec,angleVec,fluxA2,iabc(1),iabc(2),iabc(3),angleVec(i));
%     flujob2(i) = interpn(iaVec,ibVec,icVec,angleVec,fluxB2,iabc(1),iabc(2),iabc(3),angleVec(i));
%     flujoc2(i) = interpn(iaVec,ibVec,icVec,angleVec,fluxC2,iabc(1),iabc(2),iabc(3),angleVec(i));
% end
% 
% % Flujos de fase 
% figure
% % plot(angleVec,fa,'b',angleVec,fb,'r',angleVec,fc,'g')
% hold on
% plot(angleVec,flujoa,angleVec,flujob,angleVec,flujoc)
% % figure
% % plot(angleVec,flujoa1_1,'b',angleVec,flujob1_1,'r',angleVec,flujoc1_1,'g')
% % hold on
% plot(angleVec,flujoa1,'--',angleVec,flujob1,'--',angleVec,flujoc1,'--')
% plot(angleVec,flujoa2,':',angleVec,flujob2,':',angleVec,flujoc2,':')
% grid on
% 
% %% Flux partial derivatives calculation
% % A phase
% [dFAdA,dFAdB,dFAdC,dFAdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxA2);
% % B phase
% [dFBdA,dFBdB,dFBdC,dFBdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxB2);
% % C phase
% [dFCdA,dFCdB,dFCdC,dFCdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxC2);

%% Initial state of the simulator
close all
Ts = 1e-6;
tsim = 0.15;

% LUTs inversas de corrientes de Juan
load('pmsm_phi2i_21p');

sim('ee_import_fem_maxwell')

time = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.angular_velocity.series.time;
% torque_sim = load_torque.signals.values;
torque_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.electrical_torque.series.values;
thetar_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.angular_position.series.values;     % [deg] rotor angle, from wr
omegar_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.angular_velocity.series.values;     % [rpm] rotor mechanical speed wr
i_a_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.i_a.series.values;
i_b_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.i_b.series.values;
i_c_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.i_c.series.values;
% i_a_sim = phase_curr.signals.values(:,1);
% i_b_sim = phase_curr.signals.values(:,2);
% i_c_sim = phase_curr.signals.values(:,3);
v_a_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.v_a.series.values;
v_b_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.v_b.series.values;
v_c_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.v_c.series.values;
% v_a_sim = phase_volt.signals.values(:,1);
% v_b_sim = phase_volt.signals.values(:,2);
% v_c_sim = phase_volt.signals.values(:,3);

nd = length(time);

% Cambio de formato de thetar
thetar_aux = thetar_sim;
thetar_sim_nuevo = zeros(1,nd);
for i = 1:nd
    if thetar_aux(i) > 360/3/N
        thetar_aux = thetar_aux - 360/3/N;
        thetar_sim_nuevo(i) = thetar_aux(i);
    else
        thetar_sim_nuevo(i) = thetar_aux(i);
    end
end

% %% Intento de simulación discreta del sistema conocidas iabc (prueba de flujo y par)
% var_sim_iabc = zeros(3,nd);
% var_sim_iabc(:,1) = [i_a_sim(1);i_b_sim(1);i_c_sim(1)];
% var_sim_vabc = zeros(3,nd);
% var_sim_idq = zeros(3,nd);
% var_sim_fd = zeros(1,nd);
% var_sim_fq = zeros(1,nd);
% var_sim_torque = zeros(1,nd);
% fluxd = zeros(1,nd);
% fluxq = zeros(1,nd);
% 
% omega_r = RPM*2*pi/60;  % [rad/s]
% 
% method = 'linear';
% 
% % Inicialización método x=flux
% % Park
% thetar_deg = thetar_sim_nuevo(1);
% thetae = N*thetar_sim(1)*pi/180;
% 
% T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
%           sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
%           1/2            1/2                     1/2];
% invT = inv(T);
% i_dqz(:,1) = T*var_sim_iabc(:,1);
% 
% var_sim_fd(1) = interpn(idVec,iqVec',angleVec,fluxD,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);
% var_sim_fq(1) = interpn(idVec,iqVec',angleVec,fluxQ,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);
% 
% for i = 2:nd
%     thetar_deg = thetar_sim_nuevo(i-1);
%     thetae = N*thetar_sim(i-1)*pi/180;
%     
%     v_a = v_a_sim(i-1);
%     v_b = v_b_sim(i-1);
%     v_c = v_c_sim(i-1);
%     var_sim_vabc(:,i-1) = [v_a;v_b;v_c];
% 
%     var_sim_iabc(:,i-1) = [i_a_sim(i-1);i_b_sim(i-1);i_c_sim(i-1)];
% 
%     param_sim_h = time(i)-time(i-1);
%     
%     % Park
%     T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
%               sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
%               1/2            1/2                     1/2];
%     invT = inv(T);
% 
%     % MODELO EN DERIVADAS x=flux  
%     var_sim_vdq = T*var_sim_vabc(:,i-1);
%     var_sim_idq(:,i-1) = T*var_sim_iabc(:,i-1);
%     
%     fluxd(i-1) = interpn(idVec,iqVec',angleVec,fluxD,var_sim_idq(1,i-1),var_sim_idq(2,i-1),thetar_deg,method);
%     fluxq(i-1) = interpn(idVec,iqVec',angleVec,fluxQ,var_sim_idq(1,i-1),var_sim_idq(2,i-1),thetar_deg,method);
% 
%     var_sim_fd(i) = var_sim_fd(i-1) + param_sim_h*(var_sim_vdq(1) - R_s*var_sim_idq(1,i-1) - N*omega_r*var_sim_fq(i-1));
%     var_sim_fq(i) = var_sim_fq(i-1) + param_sim_h*(var_sim_vdq(2) - R_s*var_sim_idq(2,i-1) + N*omega_r*var_sim_fd(i-1));
%     
%     var_sim_torque(i-1) = interpn(idVec,iqVec',angleVec,torque,var_sim_idq(1,i-1),var_sim_idq(2,i-1),thetar_deg,method);
%     
% end
% 
% figure
% hold on
% plot(time,var_sim_idq)
% grid on
% 
% figure
% hold on
% plot(time,var_sim_fd,'--',time,var_sim_fq,'--')
% plot(time,fluxd,time,fluxq)
% grid on
% 
% figure
% hold on
% plot(time,var_sim_torque,'--')
% plot(time,torque_sim)
% grid on
%
% CONCLUSIÓN: si se inicia la iteración en un punto intermedio de la
% simulación, sí sigue bien el flujo.

% %% Intento de simulación discreta del sistema conocidas vabc
% 
% % LUTs inversas de corrientes mías
% load('currentsDQ_LUT');
% %%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% 
% % % LUTs inversas de corrientes de Juan
% % load('pmsm_phi2i_21p');
% %%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% 
% % LUTs inversas de las RN
% % Para la RN
% % Nombres de los archivos
% filename1 = 'DatosRN_completo.xlsx';
% filename2 = 'Results_completo.xlsx';
% 
% % Pesos y biases capa 1
% W1 = xlsread(filename2, 'W1');
% b1 = xlsread(filename2, 'b1');
% 
% % Pesos y biases capa 2
% W2 = xlsread(filename2, 'W2');
% b2 = xlsread(filename2, 'b2');
% 
% % Pesos y biases capa de salida
% Wout = xlsread(filename2, 'Wout');
% bout = xlsread(filename2, 'bout');
% %%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% 
% var_sim_iabc_mio = zeros(3,nd);
% var_sim_iabc_mio(:,1) = [i_a_sim(1);i_b_sim(1);i_c_sim(1)];
% var_sim_iabc_ideal = var_sim_iabc_mio;
% var_sim_iabc_rn = var_sim_iabc_mio;
% var_sim_idq_ideal = zeros(2,nd);
% var_sim_iabc_juan = var_sim_iabc_mio;
% i_dqz_mio = zeros(3,nd);
% var_sim_fd_mio = zeros(1,nd);
% var_sim_fq_mio = zeros(1,nd);
% i_dqz_rn = zeros(3,nd);
% var_sim_fdq_juan = zeros(2,nd);
% var_sim_fd_juan = zeros(1,nd);
% var_sim_fq_juan = zeros(1,nd);
% i_dqz_juan = zeros(3,nd);
% var_sim_fd_rn = zeros(1,nd);
% var_sim_fq_rn = zeros(1,nd);
% var_sim_vabc = zeros(3,nd);
% var_sim_idq = zeros(3,nd);
% var_sim_vdq = zeros(3,nd);
% 
% omega_r = RPM*2*pi/60;  % [rad/s]
% 
% method = 'linear';
% 
% % Parámetros modelo ideal
% phi_m = 0.1571;
% Ld = 0.0013;
% Lq = 0.0039;
% 
% % Inicialización método x=flux
% % Park
% thetar_deg = thetar_sim_nuevo(1);
% thetae = N*thetar_sim(1)*pi/180;
% 
% T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
%           sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
%           1/2            1/2                     1/2];
% invT = inv(T);
% i_dqz_mio(:,1) = T*var_sim_iabc_mio(:,1);
% i_dqz_rn(:,1) = i_dqz_mio(:,1);
% 
% var_sim_fd_mio(1) = interpn(idVec,iqVec',angleVec,fluxD,i_dqz_mio(1,1),i_dqz_mio(2,1),thetar_deg,method);
% var_sim_fq_mio(1) = interpn(idVec,iqVec',angleVec,fluxQ,i_dqz_mio(1,1),i_dqz_mio(2,1),thetar_deg,method);
% var_sim_fd_rn(1) = var_sim_fd_mio(1);
% var_sim_fq_rn(1) = var_sim_fq_mio(1);
% var_sim_idq_ideal(:,1) = i_dqz_mio(1:2,1);
% 
% for i = 2:nd
%     thetar_deg = thetar_sim_nuevo(i-1);
%     thetae = N*thetar_sim(i-1)*pi/180;
%     
%     v_a = v_a_sim(i-1);
%     v_b = v_b_sim(i-1);
%     v_c = v_c_sim(i-1);
% 
%     var_sim_vabc(:,i-1) = [v_a;v_b;v_c];
% 
%     param_sim_h = time(i)-time(i-1);
%     
%     % Park
%     T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
%               sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
%               1/2            1/2                     1/2];
%     invT = inv(T);
% 
%     var_sim_vdq(:,i) = T*var_sim_vabc(:,i-1);
%     var_sim_idq(:,i) = T*[i_a_sim(i);i_b_sim(i);i_c_sim(i)];
% 
%     % MODELO EN DERIVADAS x=flux (MÍO)  
%     var_sim_fd_mio(i) = var_sim_fd_mio(i-1) + param_sim_h*(var_sim_vdq(1,i) - R_s*i_dqz_mio(1,i-1) - N*omega_r*var_sim_fq_mio(i-1));
%     var_sim_fq_mio(i) = var_sim_fq_mio(i-1) + param_sim_h*(var_sim_vdq(2,i) - R_s*i_dqz_mio(2,i-1) + N*omega_r*var_sim_fd_mio(i-1));
% 
%     i_dqz_mio(1,i) = interpn(fdVec,fqVec',angleVec,iD,var_sim_fd_mio(i),var_sim_fq_mio(i),thetar_deg,method);
%     i_dqz_mio(2,i) = interpn(fdVec,fqVec',angleVec,iQ,var_sim_fd_mio(i),var_sim_fq_mio(i),thetar_deg,method);
% 
%     var_sim_iabc_mio(:,i) = invT*[i_dqz_mio(1,i);i_dqz_mio(2,i);0];
% 
%     % MODELO EN DERIVADAS x=flux (JUAN)  
%     var_sim_fdq_juan(:,i) = var_sim_fdq_juan(:,i-1) + param_sim_h*(var_sim_vdq(1:2,i) - R_s*i_dqz_juan(1:2,i-1) + N*omega_r*[0 -1;1 0]*var_sim_fdq_juan(:,i-1));
%     var_sim_fd_juan(i) = var_sim_fdq_juan(1,i);
%     var_sim_fq_juan(i) = var_sim_fdq_juan(2,i);
% %     var_sim_fd_juan(i) = var_sim_fd_juan(i-1) + param_sim_h*(var_sim_vdq(1,i) - R_s*i_dqz_juan(1,i-1) - N*omega_r*var_sim_fq_juan(i-1));
% %     var_sim_fq_juan(i) = var_sim_fq_juan(i-1) + param_sim_h*(var_sim_vdq(2,i) - R_s*i_dqz_juan(2,i-1) + N*omega_r*var_sim_fd_juan(i-1));
% 
%     i_dqz_juan(1,i) = interpn(data.phi_d,data.phi_q',data.theta,data.i_d,var_sim_fd_juan(i),var_sim_fq_juan(i),thetar_deg,method);
%     i_dqz_juan(2,i) = interpn(data.phi_d,data.phi_q',data.theta,data.i_q,var_sim_fd_juan(i),var_sim_fq_juan(i),thetar_deg,method);
% 
%     var_sim_iabc_juan(:,i) = invT*[i_dqz_juan(1,i);i_dqz_juan(2,i);0];
% 
%     % MODELO EN DERIVADAS x=flux (RN)
%     var_sim_fd_rn(i) = var_sim_fd_rn(i-1) + param_sim_h*(var_sim_vdq(1,i) - R_s*i_dqz_rn(1,i-1) - N*omega_r*var_sim_fq_rn(i-1));
%     var_sim_fq_rn(i) = var_sim_fq_rn(i-1) + param_sim_h*(var_sim_vdq(2,i) - R_s*i_dqz_rn(2,i-1) + N*omega_r*var_sim_fd_rn(i-1));
% 
%     % Aplicamos la RN
%     entrada = [var_sim_fd_rn(i),var_sim_fq_rn(i),thetar_deg];
%     % Salida capa 1
%     salida_RN_capa1 = max(entrada*W1+b1', 0);
%     % Salida capa 2
%     salida_RN_capa2 = max(salida_RN_capa1*W2+b2', 0);
%     % Salida de la RN
%     salida_RN = 300*(salida_RN_capa2*Wout+bout');
%     id_RN = salida_RN(1);
%     iq_RN = salida_RN(2);
% 
%     i_dqz_rn(1,i) = id_RN;
%     i_dqz_rn(2,i) = iq_RN;
%     
%     var_sim_iabc_rn(:,i) = invT*[i_dqz_rn(1,i);i_dqz_rn(2,i);0];
%     
%     % MODELO IDEAL
%     var_sim_idq_ideal(1,i) = var_sim_idq_ideal(1,i-1) + (param_sim_h/Ld)*(var_sim_vdq(1,i) - R_s*var_sim_idq_ideal(1,i-1) - N*omega_r*Lq*var_sim_idq_ideal(2,i-1));
%     var_sim_idq_ideal(2,i) = var_sim_idq_ideal(2,i-1) + (param_sim_h/Lq)*(var_sim_vdq(2,i) - R_s*var_sim_idq_ideal(2,i-1) + N*omega_r*(Ld*var_sim_idq_ideal(1,i-1) + phi_m));
%     var_sim_iabc_ideal(:,i) = invT*[var_sim_idq_ideal(:,i);0];
% end
% 
% figure
% hold on
% plot(time,i_a_sim,'--')
% plot(time,var_sim_iabc_mio(1,:))
% plot(time,var_sim_iabc_juan(1,:))
% plot(time,var_sim_iabc_rn(1,:))
% plot(time,var_sim_iabc_ideal(1,:))
% grid on
% legend('i_{a}^{real}','i_{a}^{mio}','i_{a}^{juan}','i_{a}^{RN}','i_{a}^{ideal}')
% 
% figure
% hold on
% plot(time,i_b_sim,'--')
% plot(time,var_sim_iabc_mio(2,:))
% plot(time,var_sim_iabc_juan(2,:))
% plot(time,var_sim_iabc_rn(2,:))
% plot(time,var_sim_iabc_ideal(2,:))
% grid on
% legend('i_{b}^{real}','i_{b}^{mio}','i_{b}^{juan}','i_{b}^{RN}','i_{b}^{ideal}')
% 
% figure
% hold on
% plot(time,i_c_sim,'--')
% plot(time,var_sim_iabc_mio(3,:))
% plot(time,var_sim_iabc_juan(3,:))
% plot(time,var_sim_iabc_rn(3,:))
% plot(time,var_sim_iabc_ideal(3,:))
% grid on
% legend('i_{c}^{real}','i_{c}^{mio}','i_{c}^{juan}','i_{c}^{RN}','i_{c}^{ideal}')
% 
% figure
% hold on
% plot(time,var_sim_idq(1,:),'--')
% plot(time,i_dqz_mio(1,:))
% plot(time,i_dqz_juan(1,:))
% plot(time,i_dqz_rn(1,:))
% plot(time,var_sim_idq_ideal(1,:))
% grid on
% legend('i_{d}^{real}','i_{d}^{mio}','i_{d}^{juan}','i_{d}^{RN}','i_{d}^{ideal}')
% 
% figure
% hold on
% plot(time,var_sim_idq(2,:),'--')
% plot(time,i_dqz_mio(2,:))
% plot(time,i_dqz_juan(2,:))
% plot(time,i_dqz_rn(2,:))
% plot(time,var_sim_idq_ideal(2,:))
% grid on
% legend('i_{q}^{real}','i_{q}^{mio}','i_{q}^{juan}','i_{q}^{RN}','i_{q}^{ideal}')
% 
% figure
% plot(time,var_sim_vdq(1:2,:))
% grid on
% legend('v_d','v_q')

%% Simulación de la PMSM estimando thetar
var_sim_iabc = zeros(3,nd);
var_sim_iabc(:,1) = [i_a_sim(1);i_b_sim(1);i_c_sim(1)];
i_dqz = zeros(3,nd);
i_dqz_MAT = zeros(3,nd);
var_sim_fdq = zeros(2,nd);
var_sim_vabc = zeros(3,nd);
var_sim_idq = zeros(3,nd);
var_sim_vdq = zeros(3,nd);

omega_r = RPM*2*pi/60;  % [rad/s]

method = 'linear';

% Inicialización método x=flux
thetar_deg = 0;
thetar = 0;
thetar_ant = 0;
thetar_deg_ant = 0;

% % Park
% thetar_deg = thetar_sim_nuevo(1);
% thetae = N*thetar_sim(1)*pi/180;
% 
% T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
%           sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
%           1/2            1/2                     1/2];
% invT = inv(T);
% i_dqz(:,1) = T*var_sim_iabc(:,1);
% 
% var_sim_fdq_mio(1,1) = interpn(idVec,iqVec',angleVec,fluxD,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);
% var_sim_fdq_mio(2,1) = interpn(idVec,iqVec',angleVec,fluxQ,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);

for i = 1:nd-1
%     thetar_deg = thetar_sim_nuevo(i-1);
%     thetae = N*thetar_sim(i-1)*pi/180;

    param_sim_h = time(i+1)-time(i);

    var_sim_vabc(:,i) = [v_a_sim(i);v_b_sim(i);v_c_sim(i)];

    % Park
    T = 2/3*[ cos(N*thetar)    cos(N*thetar - 2*pi/3)    cos(N*thetar + 2*pi/3);
              sin(N*thetar)    sin(N*thetar - 2*pi/3)    sin(N*thetar + 2*pi/3);
              1/2              1/2                       1/2];

    var_sim_vdq(:,i) = T*var_sim_vabc(:,i);
    var_sim_idq(:,i) = T*[i_a_sim(i);i_b_sim(i);i_c_sim(i)];

    % MODELO EN DERIVADAS x=flux (JUAN)  
    var_sim_fdq(:,i+1) = var_sim_fdq(:,i) + param_sim_h*(var_sim_vdq(1:2,i) - R_s*i_dqz(1:2,i) + N*omega_r*[0 -1;1 0]*var_sim_fdq(:,i));
    
    % Estimación del ángulo del rotor para h+1
    dthetar_dt = omega_r*param_sim_h;
    thetar = thetar_ant + dthetar_dt;
    thetar_ant = thetar;
    thetar_deg = thetar_deg_ant + dthetar_dt*180/pi;
    if thetar_deg > 360/3/N
        thetar_deg = thetar_deg - 360/3/N;
    end
    thetar_deg_ant = thetar_deg;

    i_dqz_MAT(1,i+1) = interpn(data.phi_d,data.phi_q',data.theta,data.i_d,var_sim_fdq(1,i+1),var_sim_fdq(2,i+1),thetar_deg,method);
    i_dqz_MAT(2,i+1) = interpn(data.phi_d,data.phi_q',data.theta,data.i_q,var_sim_fdq(1,i+1),var_sim_fdq(2,i+1),thetar_deg,method);
    
    i_dqz(1,i+1) = trilinear_interp(var_sim_fdq(1,i+1),var_sim_fdq(2,i+1),thetar_deg,data.phi_d,data.phi_q,data.theta,data.i_d);
    i_dqz(2,i+1) = trilinear_interp(var_sim_fdq(1,i+1),var_sim_fdq(2,i+1),thetar_deg,data.phi_d,data.phi_q,data.theta,data.i_q);
    
    % Park
    T = 2/3*[ cos(N*thetar)    cos(N*thetar - 2*pi/3)    cos(N*thetar + 2*pi/3);
              sin(N*thetar)    sin(N*thetar - 2*pi/3)    sin(N*thetar + 2*pi/3);
              1/2              1/2                       1/2];
    invT = inv(T);

    var_sim_iabc(:,i+1) = invT*[i_dqz(1,i+1);i_dqz(2,i+1);0];

end

figure
hold on
plot(time,i_a_sim,'--')
plot(time,var_sim_iabc(1,:))
grid on
legend('i_{a}^{real}','i_{a}^{cal}')

figure
hold on
plot(time,i_b_sim,'--')
plot(time,var_sim_iabc(2,:))
grid on
legend('i_{b}^{real}','i_{b}^{cal}')

figure
hold on
plot(time,i_c_sim,'--')
plot(time,var_sim_iabc(3,:))
grid on
legend('i_{c}^{real}','i_{c}^{cal}')

figure
hold on
plot(time,var_sim_idq(1,:),'--')
plot(time,i_dqz(1,:))
plot(time,i_dqz_MAT(1,:))
grid on
legend('i_{d}^{real}','i_{d}^{cal}','i_{d}^{cal} MAT')

figure
hold on
plot(time,var_sim_idq(2,:),'--')
plot(time,i_dqz(2,:))
plot(time,i_dqz_MAT(2,:))
grid on
legend('i_{q}^{real}','i_{q}^{cal}','i_{q}^{cal} MAT')

figure
plot(time,var_sim_vdq(1:2,:))
grid on
legend('v_d','v_q')


% % LUTs inversas de las RN
% % Para la RN
% % Nombres de los archivos
% filename1 = 'DatosRN_completo.xlsx';
% filename2 = 'Results_completo.xlsx';
% 
% % Pesos y biases capa 1
% W1 = xlsread(filename2, 'W1');
% b1 = xlsread(filename2, 'b1');
% 
% % Pesos y biases capa 2
% W2 = xlsread(filename2, 'W2');
% b2 = xlsread(filename2, 'b2');
% 
% % Pesos y biases capa de salida
% Wout = xlsread(filename2, 'Wout');
% bout = xlsread(filename2, 'bout');
% %%&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% 
% var_sim_iabc = zeros(3,nd);
% var_sim_iabc(:,1) = [i_a_sim(1);i_b_sim(1);i_c_sim(1)];
% var_sim_idq = zeros(3,nd);
% var_sim_vabc = zeros(3,nd);
% i_dqz = zeros(3,nd);
% var_sim_fd = zeros(1,nd);
% var_sim_fq = zeros(1,nd);
% 
% omega_r = RPM*2*pi/60;  % [rad/s]
% 
% method = 'linear';
% 
% % Inicialización método x=flux
% % Park
% thetar_deg = thetar_sim_nuevo(1);
% thetae = N*thetar_sim(1)*pi/180;
% 
% T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
%           sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
%           1/2            1/2                     1/2];
% invT = inv(T);
% i_dqz(:,1) = T*var_sim_iabc(:,1);
% 
% var_sim_fd(1) = interpn(idVec,iqVec',angleVec,fluxD,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);
% var_sim_fq(1) = interpn(idVec,iqVec',angleVec,fluxQ,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);
% 
% for i = 2:nd
%     thetar_deg = thetar_sim_nuevo(i-1);
%     thetae = N*thetar_sim(i-1)*pi/180;
%     
%     v_a = v_a_sim(i-1);
%     v_b = v_b_sim(i-1);
%     v_c = v_c_sim(i-1);
% 
%     var_sim_vabc(:,i-1) = [v_a;v_b;v_c];
% 
%     param_sim_h = time(i)-time(i-1);
%     
%     % Park
%     T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
%               sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
%               1/2            1/2                     1/2];
%     invT = inv(T);
% 
%     % MODELO EN DERIVADAS x=flux  
%     var_sim_vdq = T*var_sim_vabc(:,i-1);
% 
%     var_sim_fd(i) = var_sim_fd(i-1) + param_sim_h*(var_sim_vdq(1) - R_s*i_dqz(1,i-1) - N*omega_r*var_sim_fq(i-1));
%     var_sim_fq(i) = var_sim_fq(i-1) + param_sim_h*(var_sim_vdq(2) - R_s*i_dqz(2,i-1) + N*omega_r*var_sim_fd(i-1));
% 
%     % Aplicamos la RN
%     entrada = [var_sim_fd(i),var_sim_fq(i),thetar_deg];
%     % Salida capa 1
%     salida_RN_capa1 = max(entrada*W1+b1', 0);
%     % Salida capa 2
%     salida_RN_capa2 = max(salida_RN_capa1*W2+b2', 0);
%     % Salida de la RN
%     salida_RN = 300*(salida_RN_capa2*Wout+bout');
%     id_RN = salida_RN(1);
%     iq_RN = salida_RN(2);
% 
%     i_dqz(1,i) = id_RN;
%     i_dqz(2,i) = iq_RN;
%     
%     var_sim_iabc(:,i) = invT*[i_dqz(1,i);i_dqz(2,i);0];
%     var_sim_idq(:,i) = T*[i_a_sim(i);i_b_sim(i);i_c_sim(i)];
%     
% end
% 
% figure
% plot(time,var_sim_iabc)
% hold on
% plot(time,i_a_sim,'--',time,i_b_sim,'--',time,i_c_sim,'--')
% grid on
% legend('i_{a}','i_{b}','i_{c}','i_{a}^{mod}','i_{b}^{mod}','i_{c}^{mod}')
% 
% figure
% plot(time,i_dqz(1:2,:))
% hold on
% plot(time,var_sim_idq(1,:),'--',time,var_sim_idq(2,:),'--')
% grid on
% legend('i_{d}','i_{q}','i_{d}^{mod}','i_{q}^{mod}')

