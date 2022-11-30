clear all
close all
clc

% %% LUT inversa obtenida a partir de RN
% % Nombres de los archivos
% filename1 = 'DatosRN.xlsx';
% filename2 = 'Results.xlsx';
% 
% % Entradas (flujos) y salidas (intensidades)
% entrada = xlsread(filename1, 'Entradas');
% salida = xlsread(filename1, 'Salidas');
% 
% % Flujos
% flujoD = entrada(:,1);
% flujoQ = entrada(:,2);
% 
% %Intensidades
% id = salida(:,1);
% iq = salida(:,2);
% 
% % Pesos y biases capa 1
% W1 = xlsread(filename2, 'W1');
% b1 = xlsread(filename2, 'b1');
% 
% %Pesos y biases capa 2
% W2 = xlsread(filename2, 'W2');
% b2 = xlsread(filename2, 'b2');
% 
% % Pesos y biases capa de salida
% Wout = xlsread(filename2, 'Wout');
% bout = xlsread(filename2, 'bout');
% 
% % Salida capa 1
% salida_RN_capa1 = max(entrada*W1+b1', 0);
% 
% % Salida capa 2
% salida_RN_capa2 = max(salida_RN_capa1*W2+b2', 0);
% 
% % Salida de la RN
% salida_RN = 300*(salida_RN_capa2*Wout+bout');
% id_RN = salida_RN(:, 1);
% iq_RN = salida_RN(:, 2);

% %% Prueba de la RN
% clear all
% close all
% clc
% 
% % Definición de la nueva matriz de entradas
% fdVec = [-0.2:0.05:0.4];
% fqVec = [-0.4:0.05:0.4];
% cont = 1;
% for i=1:length(fdVec)
%     for j=1:length(fqVec)
%         entrada(cont,1) = fdVec(i);
%         entrada(cont,2) = fqVec(j);
%         cont = cont+1;
%     end
% end
% 
% % Ejecución de la RN
% % Nombres de los archivos
% filename1 = 'DatosRN.xlsx';
% filename2 = 'Results.xlsx';
% 
% % Flujos
% flujoD = entrada(:,1);
% flujoQ = entrada(:,2);
% 
% % Pesos y biases capa 1
% W1 = xlsread(filename2, 'W1');
% b1 = xlsread(filename2, 'b1');
% 
% %Pesos y biases capa 2
% W2 = xlsread(filename2, 'W2');
% b2 = xlsread(filename2, 'b2');
% 
% % Pesos y biases capa de salida
% Wout = xlsread(filename2, 'Wout');
% bout = xlsread(filename2, 'bout');
% 
% % Salida capa 1
% salida_RN_capa1 = max(entrada*W1+b1', 0);
% 
% % Salida capa 2
% salida_RN_capa2 = max(salida_RN_capa1*W2+b2', 0);
% 
% % Salida de la RN
% salida_RN = 300*(salida_RN_capa2*Wout+bout');
% id_RN = salida_RN(:, 1);
% iq_RN = salida_RN(:, 2);
% 
% %% LUT inversa obtenida mediante optimización e iteración
% load('currentsDQ_LUT_new');
% 
% %% Comparativa de superficies para thetar = 0
% tamletra = 10;
% tamlinea = 1;
% font = 'Times';
% tammarker = 4;
% 
% % Id
% fig1 = figure;
% set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
% hold on
% surf(fqVec,fdVec,iD(:,:,1),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',0.5)
% scatter3(flujoQ, flujoD, id_RN, 'fill');
% ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
% xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
% zlabel('{\iti_d}','FontSize',tamletra, 'FontName', font)
% view([45 25])
% %     zlim([0.05 0.25])
% %     xlim([0 800])
% grid
% colormap('parula')
% set(fig1,'Position',[20 50 360 270])
% 
% % Iq
% fig1 = figure;
% set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
% hold on
% surf(fqVec,fdVec,iQ(:,:,1),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',0.5)
% scatter3(flujoQ, flujoD, iq_RN, 'fill');
% ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
% xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
% zlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
% view([215 25])
% %     zlim([0.05 0.25])
% %     xlim([0 800])
% grid
% colormap('parula')
% set(fig1,'Position',[20 50 360 270])

%% NUEVA RN INCLUYENDO EL ÁNGULO DE POSICIÓN
% % Nombres de los archivos
% filename1 = 'DatosRN_completo.xlsx';
% filename2 = 'Results_completo.xlsx';
% 
% % Entradas (flujos) y salidas (intensidades)
% entrada = xlsread(filename1, 'Entradas');
% salida = xlsread(filename1, 'Salidas');
% 
% % Flujos
% flujoD = entrada(:,1);
% flujoQ = entrada(:,2);
% thetar = entrada(:,3);
% 
% % Intensidades
% id = salida(:,1);
% iq = salida(:,2);
% 
% % Pesos y biases capa 1
% W1 = xlsread(filename2, 'W1');
% b1 = xlsread(filename2, 'b1');
% 
% %Pesos y biases capa 2
% W2 = xlsread(filename2, 'W2');
% b2 = xlsread(filename2, 'b2');
% 
% % Pesos y biases capa de salida
% Wout = xlsread(filename2, 'Wout');
% bout = xlsread(filename2, 'bout');
% 
% % Salida capa 1
% salida_RN_capa1 = max(entrada*W1+b1', 0);
% 
% % Salida capa 2
% salida_RN_capa2 = max(salida_RN_capa1*W2+b2', 0);
% 
% % Salida de la RN
% salida_RN = 300*(salida_RN_capa2*Wout+bout');
% id_RN = salida_RN(:, 1);
% iq_RN = salida_RN(:, 2);

%% Prueba de la RN
% clear all
% close all
% clc
% 
% % Definición de la nueva matriz de entradas
% fdVec = [-0.2:0.05:0.4];
% fqVec = [-0.4:0.05:0.4];
% xVec = [0:1:90];
% cont = 1;
% for k=1:length(xVec)
%     for i=1:length(fdVec)
%         for j=1:length(fqVec)
%             entrada(cont,1) = fdVec(i);
%             entrada(cont,2) = fqVec(j);
%             entrada(cont,3) = xVec(k);
%             cont = cont+1;
%         end
%     end
% end
% 
% % Ejecución de la RN
% % Nombres de los archivos
% filename1 = 'DatosRN_completo.xlsx';
% filename2 = 'Results_completo.xlsx';
% 
% % Flujos
% flujoD = entrada(:,1);
% flujoQ = entrada(:,2);
% thetar = entrada(:,3);
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
% 
% % Salida capa 1
% salida_RN_capa1 = max(entrada*W1+b1', 0);
% 
% % Salida capa 2
% salida_RN_capa2 = max(salida_RN_capa1*W2+b2', 0);
% 
% % Salida de la RN
% salida_RN = 300*(salida_RN_capa2*Wout+bout');
% id_RN = salida_RN(:, 1);
% iq_RN = salida_RN(:, 2);
% 
% %% LUT inversa obtenida mediante optimización e iteración
% load('currentsDQ_LUT_new');
% 
% %% Comparativa de superficies para thetar = 0
% tamletra = 10;
% tamlinea = 1;
% font = 'Times';
% tammarker = 4;
% 
% for x=1:5:length(xVec)
% 
%     index = find(thetar==xVec(x));
% 
%     % Id
%     fig1 = figure;
%     set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
%     hold on
%     surf(fqVec,fdVec,iD(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',0.5)
%     scatter3(flujoQ(index), flujoD(index), id_RN(index), 'fill');
%     ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\iti_d}','FontSize',tamletra, 'FontName', font)
%     view([45 25])
%     %     zlim([0.05 0.25])
%     %     xlim([0 800])
%     grid
%     colormap('parula')
%     set(fig1,'Position',[20 50 360 270])
%     
%     % Iq
%     fig1 = figure;
%     set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
%     hold on
%     surf(fqVec,fdVec,iQ(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',0.5)
%     scatter3(flujoQ(index), flujoD(index), iq_RN(index), 'fill');
%     ylabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
%     xlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
%     zlabel('{\iti_q}','FontSize',tamletra, 'FontName', font)
%     view([215 25])
%     %     zlim([0.05 0.25])
%     %     xlim([0 800])
%     grid
%     colormap('parula')
%     set(fig1,'Position',[20 50 360 270])
% end

%% PRUEBA DE LA RN EN EL SIMULADOR
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
    if thetar_aux(i) > 360/3/N
        thetar_aux = thetar_aux - 360/3/N;
        thetar_sim_nuevo(i) = thetar_aux(i);
    else
        thetar_sim_nuevo(i) = thetar_aux(i);
    end
end

%% Ejecución del simuador
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
thetar_deg = thetar_sim(1);
thetae = N*thetar_sim(1)*pi/180;

T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
          sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
          1/2            1/2                     1/2];
invT = inv(T);
i_dqz(:,1) = T*var_sim_iabc(:,1);

% Inicialización del los flujos con las LUTs originales
var_sim_fd(1) = interpn(idVec,iqVec',angleVec,fluxD,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);
var_sim_fq(1) = interpn(idVec,iqVec',angleVec,fluxQ,i_dqz(1,1),i_dqz(2,1),thetar_deg,method);

% Para la RN
% Nombres de los archivos
filename1 = 'DatosRN_completo.xlsx';
filename2 = 'Results_completo.xlsx';

% Pesos y biases capa 1
W1 = xlsread(filename2, 'W1');
b1 = xlsread(filename2, 'b1');

% Pesos y biases capa 2
W2 = xlsread(filename2, 'W2');
b2 = xlsread(filename2, 'b2');

% Pesos y biases capa de salida
Wout = xlsread(filename2, 'Wout');
bout = xlsread(filename2, 'bout');

for i = 2:nd
    thetar_deg = thetar_sim_nuevo(i-1);
    thetae = N*thetar_sim(i-1)*pi/180;
    
    v_a = v_a_sim(i-1);
    v_b = v_b_sim(i-1);
    v_c = v_c_sim(i-1);

    var_sim_vabc(:,i-1) = [v_a;v_b;v_c];

    param_sim_h = time(i)-time(i-1);
    
    % Park
    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
              sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];
    invT = inv(T);

    % MODELO EN DERIVADAS x=flux  
    var_sim_vdq = T*var_sim_vabc(:,i-1);

    var_sim_fd(i) = var_sim_fd(i-1) + param_sim_h*(var_sim_vdq(1) - R_s*i_dqz(1,i-1) - N*omega_r*var_sim_fq(i-1));
    var_sim_fq(i) = var_sim_fq(i-1) + param_sim_h*(var_sim_vdq(2) - R_s*i_dqz(2,i-1) + N*omega_r*var_sim_fd(i-1));
    
    % Aplicamos la RN
    entrada = [var_sim_fd(i),var_sim_fq(i),thetar_deg];
    % Salida capa 1
    salida_RN_capa1 = max(entrada*W1+b1', 0);
    % Salida capa 2
    salida_RN_capa2 = max(salida_RN_capa1*W2+b2', 0);
    % Salida de la RN
    salida_RN = 300*(salida_RN_capa2*Wout+bout');
    id_RN = salida_RN(1);
    iq_RN = salida_RN(2);

    i_dqz(1,i) = id_RN;
    i_dqz(2,i) = iq_RN;

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