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

%% Obtention of new LUT: fluxA, fluxB, and fluxC as functions of iA, iB, iC, and thetar

% Extension of the fluxD, fluxQ, flux0, and torque LUTs to a maximum thetar of 360/N
% 30 --> 90
for i=1:30
    fluxD(:,:,31+i) = fluxD(:,:,i+1);
    fluxQ(:,:,31+i) = fluxQ(:,:,i+1);
    flux0(:,:,31+i) = flux0(:,:,i+1);
    torque(:,:,31+i) = torque(:,:,i+1);
end
for i=1:30
    fluxD(:,:,61+i) = fluxD(:,:,i+1);
    fluxQ(:,:,61+i) = fluxQ(:,:,i+1);
    flux0(:,:,61+i) = flux0(:,:,i+1);
    torque(:,:,61+i) = torque(:,:,i+1);
end

angleVec = [0:1:90];
niD = length(idVec);
niQ = length(iqVec);
nang = length(angleVec);

% Calulation of phase fluxes LTUs as funtions of id, iq, and thetar
fluxA = zeros(niD,niQ,nang);
fluxB = zeros(niD,niQ,nang);
fluxC = zeros(niD,niQ,nang);
for x=1:nang
    for d=1:niD
        for q=1:niQ
            % Clarke-Park matrix
            theta = N*angleVec(x)*pi/180;
            T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
                     sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
                     1/2           1/2                    1/2];
            invT = inv(T);

            % fluxD, fluxQ, flux0, thetar --> fluxA, fluxB, fluxC
            fluxfase = invT*[fluxD(d,q,x);fluxQ(d,q,x);flux0(d,q,x)];
            fluxA(d,q,x) = fluxfase(1);
            fluxB(d,q,x) = fluxfase(2);
            fluxC(d,q,x) = fluxfase(3);

%             % Phase currents
%             ifase = invT(:,1:2)*[idVec(d);iqVec(q)];
%             iA(d,q,x) = ifase(1);
%             iB(d,q,x) = ifase(2);
%             iC(d,q,x) = ifase(3);
        end
    end
end

% ifa(1,1:91)=iA(11,15,:);
% ifb(1,1:91)=iB(11,15,:);
% ifc(1,1:91)=iC(11,15,:);
% figure
% plot(angleVec,ifa,angleVec,ifb,angleVec,ifc)
% grid on
% 
% fd(1,1:91) = fluxD(11,15,:);
% fq(1,1:91) = fluxQ(11,15,:);
% figure
% plot(angleVec,fd,angleVec,fq)
% legend('\phi_{d}','\phi_{q}')
% grid on
% 
% fa(1,1:91) = fluxA(11,15,:);
% fb(1,1:91) = fluxB(11,15,:);
% fc(1,1:91) = fluxC(11,15,:);
% figure
% plot(angleVec,fa,angleVec,fb,angleVec,fc)
% legend('\phi_{a}','\phi_{b}','\phi_{c}')
% grid on

% Selection of the phase currents values for the new LUTs
% iaVec = [-300 -270 -240 -210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210 240 270 300];
% ibVec = [-300 -270 -240 -210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210 240 270 300];
% icVec = [-300 -270 -240 -210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210 240 270 300];
iaVec = [-210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210];
ibVec = [-210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210];
icVec = [-210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210];

niA = length(iaVec);
niB = length(ibVec);
niC = length(icVec);

fluxA2 = zeros(niA,niB,niC,nang);
fluxB2 = zeros(niA,niB,niC,nang);
fluxC2 = zeros(niA,niB,niC,nang);

% Calculation of LTUs for the phase fluxes in terms of iA, iB, iC and
% thetar
for x=1:nang
    for a=1:niA
        for b=1:niB
            for c=1:niC
                % Clarke-Park matrix
                theta = N*angleVec(x)*pi/180;
                T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
                         sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
                         1/2           1/2                    1/2];
                % d-q currents
                idq = T*[iaVec(a);ibVec(b);icVec(c)];
                % Interpolation form the phase fluxes LTUs obtained before
                fluxA2(a,b,c,x) = interp3(iqVec,idVec',angleVec,fluxA,idq(2),idq(1),angleVec(x));
                fluxB2(a,b,c,x) = interp3(iqVec,idVec',angleVec,fluxB,idq(2),idq(1),angleVec(x));
                fluxC2(a,b,c,x) = interp3(iqVec,idVec',angleVec,fluxC,idq(2),idq(1),angleVec(x));
            end
        end
    end
end

% idq = [idVec(10);iqVec(15)];
% for i = 1:length(angleVec)
%     theta = N*angleVec(i)*pi/180;
%     T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
%              sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
%              1/2           1/2                    1/2];
%     invT = inv(T);
%     ifase(:,i) = invT(:,1:2)*idq;
% 
%     fa2(1,i) = interpn(iaVec,ibVec,icVec,angleVec,fluxA2,ifase(1),ifase(2),ifase(3),angleVec(i));
%     fb2(1,i) = interpn(iaVec,ibVec,icVec,angleVec,fluxB2,ifase(1),ifase(2),ifase(3),angleVec(i));
%     fc2(1,i) = interpn(iaVec,ibVec,icVec,angleVec,fluxC2,ifase(1),ifase(2),ifase(3),angleVec(i));
% end


% figure
% plot(angleVec,fa,angleVec,fb,angleVec,fc)
% hold on
% plot(angleVec,fa2,':',angleVec,fb2,':',angleVec,fc2,':')
% legend('\phi_{a2}','\phi_{b2}','\phi_{c2}')
% grid on

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
%     flujod(i) = interp3(iqVec,idVec',angleVec,fluxD,iq,id,angleVec(i));
%     flujoq(i) = interp3(iqVec,idVec',angleVec,fluxQ,iq,id,angleVec(i));
%     flujo0(i) = interp3(iqVec,idVec',angleVec,flux0,iq,id,angleVec(i));
%     % Pasamos a flujo abc
%     theta = N*angleVec(i)*pi/180;
%     T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
%              sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
%              1/2           1/2                    1/2];
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
%     flujoa1(i) = interp3(iqVec,idVec',angleVec,fluxA,iq,id,angleVec(i));
%     flujob1(i) = interp3(iqVec,idVec',angleVec,fluxB,iq,id,angleVec(i));
%     flujoc1(i) = interp3(iqVec,idVec',angleVec,fluxC,iq,id,angleVec(i));
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
% plot(angleVec,flujoa,'b',angleVec,flujob,'r',angleVec,flujoc,'g')
% % figure
% % plot(angleVec,flujoa1_1,'b',angleVec,flujob1_1,'r',angleVec,flujoc1_1,'g')
% % hold on
% plot(angleVec,flujoa1,'--b',angleVec,flujob1,'--r',angleVec,flujoc1,'--g')
% plot(angleVec,flujoa2,':b',angleVec,flujob2,':r',angleVec,flujoc2,':g')
% grid on

%% Flux partial derivatives calculation
% A phase
[dFAdA,dFAdB,dFAdC,dFAdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxA2);
% B phase
[dFBdA,dFBdB,dFBdC,dFBdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxB2);
% C phase
[dFCdA,dFCdB,dFCdC,dFCdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxC2);

%% Initial state of the simulator
% Speed
omega_m = 7200*2*pi/60;   % rad/s
omega_r = omega_m*N;    % rad/s
% Phase currents
var_sim_iabc = [42.68 9.167 -51.847]';
% theta: theta_e = N*theta_r
theta_ant = 0;
theta = 0;

% Total time 
param_exp_Ttotal = 0.1;
% paso de integración
param_sim_h = 1e-6; % (s)
% número de pasos a realizar para simular el experimento
param_sim_Nh = ceil(param_exp_Ttotal/param_sim_h)+1;
% Time
tray_sim_tiempo = [0:param_sim_h:param_exp_Ttotal]';

%% Recording variables
tray_sim_time       = zeros( param_sim_Nh, 1);
tray_sim_ifase      = zeros( param_sim_Nh, 3); 
tray_sim_idq0       = zeros( param_sim_Nh, 3);
% tray_sim_speed      = zeros( param_sim_Nh, 1); 
tray_sim_torque     = zeros( param_sim_Nh, 1);
tray_sim_vabc       = zeros( param_sim_Nh, 3);
tray_sim_thetar     = zeros( param_sim_Nh, 1);
tray_sim_fabc       = zeros( param_sim_Nh, 3);
tray_sim_fdq0       = zeros( param_sim_Nh, 3);
tray_sim_fdq0_lut   = zeros( param_sim_Nh, 3);

%% Simulation of derivative-based PMSM model
for h=1:param_sim_Nh

    % Clarke-Park matrix
    T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
             sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
             1/2           1/2                    1/2];
    invT = inv(T);
    
%     % Adjusting theta angle for the LUT
%     if pi/6 < theta && theta <= pi/3
%         theta_lut = theta - pi/6;
%     elseif pi/3 < theta && theta <= pi/2
%         theta_lut = theta - pi/3;
%     elseif pi/2 < theta && theta <= 2*pi/3
%         theta_lut = theta - pi/2;
%     elseif 2*pi/3 < theta && theta <= 5*pi/6
%         theta_lut = theta - 2*pi/3;
%     elseif 5*pi/6 < theta && theta <= pi
%         theta_lut = theta - 5*pi/6;
%     elseif pi < theta && theta <= 7*pi/6
%         theta_lut = theta - pi;
%     elseif 7*pi/6 < theta && theta <= 4*pi/3
%         theta_lut = theta - 7*pi/6;
%     elseif 4*pi/3 < theta && theta <= 3*pi/2
%         theta_lut = theta - 4*pi/3;
%     elseif 3*pi/2 < theta && theta <= 5*pi/3
%         theta_lut = theta - 3*pi/2;
%     elseif 5*pi/3 < theta && theta <= 11*pi/6
%         theta_lut = theta - 5*pi/3;
%     elseif 11*pi/6 < theta && theta <= 2*pi
%         theta_lut = theta - 11*pi/6;
%     end
    
    % Rotor angle for the LUTs
    thetar = theta*180/pi/N;

    % Actual abc currents
    i_a = var_sim_iabc(1);
    i_b = var_sim_iabc(2);
    i_c = var_sim_iabc(3);

    % Actual dq0 currents
    var_sim_idqz = T*var_sim_iabc;
    i_d = var_sim_idqz(1);
    i_q = var_sim_idqz(2);
    i_0 = var_sim_idqz(3);
    
    % Actual abc flux
    flux_a = interpn(iaVec,ibVec,icVec,angleVec,fluxA2,i_a,i_b,i_c,thetar);
    flux_b = interpn(iaVec,ibVec,icVec,angleVec,fluxB2,i_a,i_b,i_c,thetar);
    flux_c = interpn(iaVec,ibVec,icVec,angleVec,fluxC2,i_a,i_b,i_c,thetar);
    var_sim_fluxabc = [flux_a; flux_b; flux_c];

    % Actual dqz flux
    var_sim_fluxdq0 = T*var_sim_fluxabc;
    flux_d = var_sim_fluxdq0(1);
    flux_q = var_sim_fluxdq0(2);
    flux_0 = var_sim_fluxdq0(3);
    
    % Actual dqz flux from LUTs
    flux_d_lut = interp3(iqVec,idVec',angleVec,fluxD,i_q,i_d,thetar);
    flux_q_lut = interp3(iqVec,idVec',angleVec,fluxQ,i_q,i_d,thetar);
    flux_0_lut = interp3(iqVec,idVec',angleVec,flux0,i_q,i_d,thetar);
    var_sim_fluxdq0_lut = [flux_d_lut;flux_q_lut;flux_0_lut];

    % Actual torque
    var_sim_torque = interp3(iqVec,idVec',angleVec,torque,i_q,i_d,thetar);
    
    % Actual stator voltages
    v_a = i_a*Rload;
    v_b = i_b*Rload;
    v_c = i_c*Rload;

    % Flux derivatives
    l_aa = interpn(iaVec,ibVec,icVec,angleVec,dFAdA,i_a,i_b,i_c,thetar);
    l_ab = interpn(iaVec,ibVec,icVec,angleVec,dFAdB,i_a,i_b,i_c,thetar);
    l_ac = interpn(iaVec,ibVec,icVec,angleVec,dFAdC,i_a,i_b,i_c,thetar);
    phi_ar = interpn(iaVec,ibVec,icVec,angleVec,dFAdX,i_a,i_b,i_c,thetar);

    l_ba = interpn(iaVec,ibVec,icVec,angleVec,dFBdA,i_a,i_b,i_c,thetar);
    l_bb = interpn(iaVec,ibVec,icVec,angleVec,dFBdB,i_a,i_b,i_c,thetar);
    l_bc = interpn(iaVec,ibVec,icVec,angleVec,dFBdC,i_a,i_b,i_c,thetar);
    phi_br = interpn(iaVec,ibVec,icVec,angleVec,dFBdX,i_a,i_b,i_c,thetar);

    l_ca = interpn(iaVec,ibVec,icVec,angleVec,dFCdA,i_a,i_b,i_c,thetar);
    l_cb = interpn(iaVec,ibVec,icVec,angleVec,dFCdB,i_a,i_b,i_c,thetar);
    l_cc = interpn(iaVec,ibVec,icVec,angleVec,dFCdC,i_a,i_b,i_c,thetar);
    phi_cr = interpn(iaVec,ibVec,icVec,angleVec,dFCdX,i_a,i_b,i_c,thetar);
    
    % System evolution
    di_a = (l_ab*l_bc*v_c - l_ab*l_cc*v_b - l_ac*l_bb*v_c + l_ac*l_cb*v_b + l_bb*l_cc*v_a - l_bc*l_cb*v_a - l_ab*l_bc*omega_r*phi_cr + l_ac*l_bb*omega_r*phi_cr + l_ab*l_cc*omega_r*phi_br - l_ac*l_cb*omega_r*phi_br - l_bb*l_cc*omega_r*phi_ar + l_bc*l_cb*omega_r*phi_ar - R_s*i_a*l_bb*l_cc + R_s*i_a*l_bc*l_cb + R_s*i_b*l_ab*l_cc - R_s*i_b*l_ac*l_cb - R_s*i_c*l_ab*l_bc + R_s*i_c*l_ac*l_bb)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
    di_b = -(l_aa*l_bc*v_c - l_aa*l_cc*v_b - l_ac*l_ba*v_c + l_ac*l_ca*v_b + l_ba*l_cc*v_a - l_bc*l_ca*v_a - l_aa*l_bc*omega_r*phi_cr + l_ac*l_ba*omega_r*phi_cr + l_aa*l_cc*omega_r*phi_br - l_ac*l_ca*omega_r*phi_br - l_ba*l_cc*omega_r*phi_ar + l_bc*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cc + R_s*i_a*l_bc*l_ca + R_s*i_b*l_aa*l_cc - R_s*i_b*l_ac*l_ca - R_s*i_c*l_aa*l_bc + R_s*i_c*l_ac*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
    di_c = (l_aa*l_bb*v_c - l_aa*l_cb*v_b - l_ab*l_ba*v_c + l_ab*l_ca*v_b + l_ba*l_cb*v_a - l_bb*l_ca*v_a - l_aa*l_bb*omega_r*phi_cr + l_ab*l_ba*omega_r*phi_cr + l_aa*l_cb*omega_r*phi_br - l_ab*l_ca*omega_r*phi_br - l_ba*l_cb*omega_r*phi_ar + l_bb*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cb + R_s*i_a*l_bb*l_ca + R_s*i_b*l_aa*l_cb - R_s*i_b*l_ab*l_ca - R_s*i_c*l_aa*l_bb + R_s*i_c*l_ab*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
    
    var_sim_iabc = var_sim_iabc + param_sim_h*[di_a;di_b;di_c];
       
    % Rotor angle estimation
    theta = theta_ant + omega_m*param_sim_h;
    if theta > 2*pi
        theta = theta - 2*pi;
    elseif theta < 0
        theta = theta + 2*pi;
    end
    theta_ant = theta;
         
    % Recordings
    tray_sim_time(h)        = h*param_sim_h;
    tray_sim_ifase(h,:)     = [i_a, i_b, i_c];
    tray_sim_idq0(h,:)      = var_sim_idqz';
%     tray_sim_speed(h,:)     = omega_m;
    tray_sim_torque(h,:)    = var_sim_torque;
    tray_sim_vabc(h,:)      = [v_a, v_b, v_c];
    tray_sim_thetar(h,:)    = thetar;
    tray_sim_fabc(h,:)      = var_sim_fluxabc;
    tray_sim_fdq0(h,:)      = var_sim_fluxdq0;
    tray_sim_fdq0_lut(h,:)  = var_sim_fluxdq0_lut;
end

%% Gráficas para guardar
close all

p1 = figure;
set(p1,'OuterPosition',[1 210 470 330]);
hold on
plot(tray_sim_time,tray_sim_ifase);
fase=legend('\it{i_{a}}','\it{i_{b}}','\it{i_{c}}','Orientation','Horizontal');
set(fase,'FontSize',10,'FontName','Times New Roman')
xlabel('Time (s)','FontSize',10)
ylabel('Current (A)','FontSize',10)
% xlim([0.9 1])
% ylim([-2 2.5])
grid on
box on

p2 = figure;
set(p2,'OuterPosition',[1 210 470 330]);
hold on
plot(tray_sim_time,tray_sim_idq0);
alfabeta=legend('\it{i_{d}}','\it{i_{q}}','\it{i_{0}}','Orientation','Horizontal');
set(alfabeta,'FontSize',10,'FontName','Times New Roman')
xlabel('Time (s)','FontSize',10)
ylabel('Current (A)','FontSize',10)
% xlim([0.9 1])
% ylim([-2 2.5])
grid on
box on

p4 = figure;
set(p4,'OuterPosition',[1 210 470 330]);
hold on
plot(tray_sim_time,tray_sim_vabc);
dd=legend('\it{v_{a}}','\it{v_{b}}','\it{v_{c}}','Orientation','Horizontal');
set(dd,'FontSize',10,'FontName','Times New Roman')
xlabel('Time (s)','FontSize',10)
ylabel('Voltage (V)','FontSize',10)
% xlim([0.5 1])
% ylim([0 2.5])
grid on
box on

p6 = figure;
set(p6,'OuterPosition',[1 210 470 330]);
hold on
plot(tray_sim_time,tray_sim_torque);
legend('T','Location','NorthEast')
xlabel('Time (s)')
ylabel('Torque (N·m)')
% xlim([0.5 1])
% ylim([499.5 500.5])
grid on
box on


% %% Prueba de LUT
% 
% for i=1:length(angleVec)
%     flujo(i) = interp3(iqVec,idVec',angleVec,fluxQ,iqVec(11),idVec(11),angleVec(i));
%     flujo_p(i) = fluxQ(11,11,i);
% end
% figure
% plot(angleVec,flujo)
% hold on
% plot(angleVec,flujo_p,'r')
% 
% 
% for i=1:length(angleVec)
%     flujo(i)=interp3(iqVec,idVec',angleVec,fluxD,iqVec(11),idVec(11),angleVec(i));
%     flujo_p(i) = fluxD(11,11,i);
% end
% figure
% plot(angleVec,flujo)
% hold on
% plot(angleVec,flujo_p,'r')

