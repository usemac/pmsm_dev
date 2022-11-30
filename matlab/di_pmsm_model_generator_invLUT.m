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

%% Representación de los flujos para un valor de thetar
% Parámetros de representación
tamletra = 10;
tamlinea = 1;
font = 'Times';
tammarker = 4;

nx  = length(angleVec);

for x=1:7:nx
    % FluxD
    fig1 = figure;
    set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
    hold on
    surf(iqVec,idVec,fluxD(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
    ylabel('{\iti_d}','FontSize',tamletra, 'FontName', font)
    xlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
    zlabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
    view([45 25])
%     zlim([0.05 0.25])
%     xlim([0 800])
    grid
    colormap('parula')
    set(fig1,'Position',[20 50 360 270])

    % FluxQ
    fig1 = figure;
    set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
    hold on
    surf(iqVec,idVec,fluxQ(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
    ylabel('{\iti_d}','FontSize',tamletra, 'FontName', font)
    xlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
    zlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
    view([215 25])
%     zlim([0.05 0.25])
%     xlim([0 800])
    grid
    colormap('parula')
    set(fig1,'Position',[20 50 360 270])
end

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
% xVec = [0:5:90];
% 
% % Inversión mediante minimización
% nx  = length(xVec);
% nfd = length(fdVec);
% nfq = length(fqVec);
% 
% idVec_aux = [-350:5:350];
% iqVec_aux = [-350:5:350];
% 
% nid = length(idVec_aux);
% niq = length(iqVec_aux);
% 
% iD = zeros(nfd,nfq,nx);
% iQ = zeros(nfd,nfq,nx);
% 
% metodo = 'makima';
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
%                     if id < idVec(1) || id > idVec(length(idVec)) || iq < iqVec(1) || iq > iqVec(length(iqVec))
%                         metodo = 'makima';
%                     else
%                         metodo = 'linear';
%                     end
%                     fluxd = interpn(idVec,iqVec',angleVec,fluxD,id,iq,xVec(x),metodo);
%                     fluxq = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,xVec(x),metodo);
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
% nx  = length(xVec);
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
%             fluxd = interpn(idVec,iqVec',angleVec,fluxD,id,iq,xVec(x),metodo);
%             fluxq = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,xVec(x),metodo);
%             id_LUT = interpn(fdVec,fqVec',xVec,iD,fluxd,fluxq,xVec(x),metodo);
%             iq_LUT = interpn(fdVec,fqVec',xVec,iQ,fluxd,fluxq,xVec(x),metodo);
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
%             id = interpn(fdVec,fqVec',xVec,iD,fluxd,fluxq,xVec(x),metodo);
%             iq = interpn(fdVec,fqVec',xVec,iQ,fluxd,fluxq,xVec(x),metodo);
%             fluxd_LUT = interpn(idVec,iqVec',angleVec,fluxD,id,iq,xVec(x),metodo);
%             fluxq_LUT = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,xVec(x),metodo);
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

% %% Inversión de las LUT de los flujos D y Q (MÉTODO 2)
% close all;
% 
% % Creación de la matriz de datos de entrada a la función
% nx  = length(angleVec);
% nid = length(idVec);
% niq = length(iqVec);
% 
% % Tendremos una matriz de nid*niq filas, 4 columnas (id,iq,fluxd,fluxq), y
% % nx páginas (1 por ángulo)
% InputPoints = zeros(nid*niq,4,nx);
% 
% metodo = 'linear';
% 
% for x=1:nx
%     cont = 1;
%     thetar = angleVec(x);
%     for index_d=1:nid
%         id = idVec(index_d);
%         for index_q=1:niq
%             iq = iqVec(index_q);
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
% % fdVec = [-0.2:0.01:0.4];
% % fqVec = [-0.4:0.01:0.4];
% fdVec = [round(fluxD_min,3):0.01:round(fluxD_max,3)];
% fqVec = [round(fluxQ_min,3):0.01:round(fluxQ_max,3)];
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
% 
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
tsim = 0.01;

sim('ee_import_fem_maxwell')

time = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.angular_velocity.series.time;
thetar_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.angular_position.series.values;     % [deg] rotor angle, from wr
omegar_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.angular_velocity.series.values;     % [rpm] rotor mechanical speed wr
torque_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.electrical_torque.series.values;
i_a_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.i_a.series.values;
i_b_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.i_b.series.values;
i_c_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.i_c.series.values;
i_d_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.i_d.series.values;
i_q_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.i_q.series.values;
v_a_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.v_a.series.values;
v_b_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.v_b.series.values;
v_c_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.v_c.series.values;
load_torque_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM1.torque.series.values;

nd = length(time);

% Cambio de formato de thetar
thetar_aux = thetar_sim;
thetar_sim_nuevo = zeros(1,nd);
for i = 1:nd
    if thetar_aux(i) > 360/N
        thetar_aux = thetar_aux - 360/N;
        thetar_sim_nuevo(i) = thetar_aux(i);
    else
        thetar_sim_nuevo(i) = thetar_aux(i);
    end
end

load('torqueDQ_LUT_new');
angleVec2=[0:90];

% Comprobación de corrientes abc y de par
for i = 1:nd
    thetae = N*thetar_sim(i)*pi/180;
    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
             -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];
    invT = inv(T);
    iabc(:,i) = invT(:,1:2)*[i_d_sim(i);i_q_sim(i)];
    idq(:,i) = T*[i_a_sim(i);i_b_sim(i);i_c_sim(i)];
    torque_cal(i) = interpn(idVec,iqVec',angleVec2,torque,i_d_sim(i),i_q_sim(i),thetar_sim_nuevo(i),"linear");
end

figure
plot(time,i_a_sim,time,i_b_sim,time,i_c_sim)
hold on
plot(time,iabc,'--')
% CONCLUSIÓN: Coinciden

figure
plot(time,i_d_sim,time,i_q_sim)
hold on
plot(time,idq,'--')
% CONCLUSIÓN: Coinciden

figure
plot(time,torque_sim)
hold on
plot(time,torque_cal,'--')
plot(time,load_torque_sim,':')
% CONCLUSIÓN: No coincide.

% %% Intento de simulación discreta del sistema conocidas vabc (prueba de flujos)
% var_sim_iabc = zeros(3,nd);
% var_sim_iabc(:,1) = [i_a_sim(1);i_b_sim(1);i_c_sim(1)];
% var_sim_vabc = zeros(3,nd);
% var_sim_idq = zeros(3,nd);
% i_dqz = zeros(3,nd);
% var_sim_fd = zeros(1,nd);
% var_sim_fq = zeros(1,nd);
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
% thetae = N*thetar_sim_nuevo(1)*pi/180;
% 
% T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
%          -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
%           1/2            1/2                     1/2];
% invT = inv(T);
% i_dqz(:,1) = T*var_sim_iabc(:,1);
% 
% var_sim_fd(1) = interpn(idVec,iqVec',angleVec,fluxD,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);
% var_sim_fq(1) = interpn(idVec,iqVec',angleVec,fluxQ,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);
% 
% for i = 2:nd
%     thetar_deg = thetar_sim_nuevo(i-1);
%     thetae = N*thetar_sim_nuevo(i-1)*pi/180;
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
%              -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
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
%     var_sim_fd(i) = var_sim_fd(i-1) + param_sim_h*(var_sim_vdq(1) - R_s*var_sim_idq(1,i-1) + N*omega_r*var_sim_fq(i-1));
%     var_sim_fq(i) = var_sim_fq(i-1) + param_sim_h*(var_sim_vdq(2) - R_s*var_sim_idq(2,i-1) - N*omega_r*var_sim_fd(i-1));
% 
%     i_dqz(1,i) = interpn(fdVec,fqVec',angleVec,iD,var_sim_fd(i),var_sim_fq(i),thetar_deg,method);
%     i_dqz(2,i) = interpn(fdVec,fqVec',angleVec,iQ,var_sim_fd(i),var_sim_fq(i),thetar_deg,method);
% 
% %     var_sim_iabc(:,i) = invT*[i_dqz(1,i);i_dqz(2,i);0];
%     
% end
% 
% % figure
% % plot(time,var_sim_iabc)
% % hold on
% % plot(time,i_a_sim,'--',time,i_b_sim,'--',time,i_c_sim,'--')
% 
% figure
% plot(time,i_dqz)
% hold on
% plot(time,var_sim_idq,'--')
% plot(time,i_d_sim,':',time,i_q_sim,':')
% 
% figure
% hold on
% plot(time,var_sim_fd,'--',time,var_sim_fq,'--')
% plot(time,fluxd,time,fluxq)
% grid on


%% Intento de simulación discreta del sistema conocidas vabc
% LUTs de corrientes
load('currentsDQ_LUT_new');
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&

var_sim_iabc = zeros(3,nd);
var_sim_iabc(:,1) = [i_a_sim(1);i_b_sim(1);i_c_sim(1)];
var_sim_vabc = zeros(3,nd);
i_dqz = zeros(3,nd);
var_sim_fd = zeros(1,nd);
var_sim_fq = zeros(1,nd);

omega_r = RPM*2*pi/60;  % [rad/s]

method = 'linear';

% Inicialización método x=flux
% Park
thetar_deg = thetar_sim_nuevo(1);
thetae = N*thetar_sim_nuevo(1)*pi/180;

T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
         -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
          1/2            1/2                     1/2];
invT = inv(T);
i_dqz(:,1) = T*var_sim_iabc(:,1);

var_sim_fd(1) = interpn(idVec,iqVec',angleVec,fluxD,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);
var_sim_fq(1) = interpn(idVec,iqVec',angleVec,fluxQ,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);

for i = 2:nd
    thetar_deg = thetar_sim_nuevo(i-1);
    thetae = N*thetar_sim_nuevo(i-1)*pi/180;
    
    v_a = v_a_sim(i-1);
    v_b = v_b_sim(i-1);
    v_c = v_c_sim(i-1);

    var_sim_vabc(:,i-1) = [v_a;v_b;v_c];

    param_sim_h = time(i)-time(i-1);
    
    % Park
    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
             -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];
    invT = inv(T);

    % MODELO EN DERIVADAS x=flux  
    var_sim_vdq = T*var_sim_vabc(:,i-1);

    var_sim_fd(i) = var_sim_fd(i-1) + param_sim_h*(var_sim_vdq(1) - R_s*i_dqz(1,i-1) + N*omega_r*var_sim_fq(i-1));
    var_sim_fq(i) = var_sim_fq(i-1) + param_sim_h*(var_sim_vdq(2) - R_s*i_dqz(2,i-1) - N*omega_r*var_sim_fd(i-1));

    i_dqz(1,i) = interpn(fdVec,fqVec',xVec,iD,var_sim_fd(i),var_sim_fq(i),thetar_deg,method);
    i_dqz(2,i) = interpn(fdVec,fqVec',xVec,iQ,var_sim_fd(i),var_sim_fq(i),thetar_deg,method);

    var_sim_iabc(:,i) = invT*[i_dqz(1,i);i_dqz(2,i);0];
    
end

figure
plot(time,var_sim_iabc)
hold on
plot(time,i_a_sim,'--',time,i_b_sim,'--',time,i_c_sim,'--')
grid on
legend('i_{a}','i_{b}','i_{c}','i_{a}^{mod}','i_{b}^{mod}','i_{c}^{mod}')

figure
plot(time,i_dqz(1:2,:))
hold on
plot(time,i_d_sim,'--',time,i_q_sim,'--')
grid on
legend('i_{d}','i_{q}','i_{d}^{mod}','i_{q}^{mod}')

