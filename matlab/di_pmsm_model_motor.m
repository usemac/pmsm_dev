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
% Rload = 10;

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

% %% TEST: representation of diferent variables for a d-q currents value.
% id = 0;
% iq = 50;
% thetar_deg = [0:1:180];
% 
% iabc = zeros(3,length(thetar_deg));
% fluxdq = zeros(3,length(thetar_deg));
% fluxabc = zeros(3,length(thetar_deg));
% torque_e = zeros(1,length(thetar_deg));
% 
% for x = 1:length(thetar_deg)
%     thetar = thetar_deg(x);
%     if thetar > 90
%         thetar = thetar - 90;
%     end
%     theta = N*thetar*pi/180;
%     
%     % Transformation matrix
%     T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
%              sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
%              1/2           1/2                    1/2];
%     invT = inv(T);
%     % Phase currents
%     iabc(:,x) = invT(:,1:2)*[id;iq];
%     % d-q fluxes
%     fluxdq(:,x) = [interp3(iqVec,idVec',angleVec,fluxD,iq,id,thetar);
%                    interp3(iqVec,idVec',angleVec,fluxQ,iq,id,thetar);
%                    interp3(iqVec,idVec',angleVec,flux0,iq,id,thetar)];
%     % a-b-c fluxes
%     fluxabc(:,x) = [interpn(iaVec,ibVec,icVec,angleVec,fluxA2,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%                     interpn(iaVec,ibVec,icVec,angleVec,fluxB2,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%                     interpn(iaVec,ibVec,icVec,angleVec,fluxC2,iabc(1,x),iabc(2,x),iabc(3,x),thetar)];
%     % torque
%     torque_e(x) = interp3(iqVec,idVec',angleVec,torque,iq,id,thetar);
% 
%     % Derivatives
%     l_aa(x) = interpn(iaVec,ibVec,icVec,angleVec,dFAdA,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%     l_ab(x) = interpn(iaVec,ibVec,icVec,angleVec,dFAdB,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%     l_ac(x) = interpn(iaVec,ibVec,icVec,angleVec,dFAdC,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%     phi_ar(x) = interpn(iaVec,ibVec,icVec,angleVec,dFAdX,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
% 
%     l_ba(x) = interpn(iaVec,ibVec,icVec,angleVec,dFBdA,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%     l_bb(x) = interpn(iaVec,ibVec,icVec,angleVec,dFBdB,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%     l_bc(x) = interpn(iaVec,ibVec,icVec,angleVec,dFBdC,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%     phi_br(x) = interpn(iaVec,ibVec,icVec,angleVec,dFBdX,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
% 
%     l_ca(x) = interpn(iaVec,ibVec,icVec,angleVec,dFCdA,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%     l_cb(x) = interpn(iaVec,ibVec,icVec,angleVec,dFCdB,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%     l_cc(x) = interpn(iaVec,ibVec,icVec,angleVec,dFCdC,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
%     phi_cr(x) = interpn(iaVec,ibVec,icVec,angleVec,dFCdX,iabc(1,x),iabc(2,x),iabc(3,x),thetar);
% 
% end
% 
% p2 = figure;
% set(p2,'OuterPosition',[1 210 470 330]);
% hold on
% plot(thetar_deg,l_aa,thetar_deg,l_bb,thetar_deg,l_cc);
% legend('di_{a}/di_{a}','di_{b}/di_{a}','di_{c}/di_{a}','Location','NorthEast')
% xlabel(['\Theta_{r}', '(degree)'])
% ylabel('di_{abc}/di_{a}')
% % xlim([0.5 1])
% % ylim([499.5 500.5])
% grid on
% box on
% 
% p3 = figure;
% set(p3,'OuterPosition',[1 210 470 330]);
% hold on
% plot(thetar_deg,iabc);
% legend('i_{a}','i_{b}','i_{c}','Location','NorthEast')
% xlabel(['\Theta_{r}', '(degree)'])
% ylabel('i_{abc} (Wb)')
% % xlim([0.5 1])
% % ylim([499.5 500.5])
% grid on
% box on
% 
% p4 = figure;
% set(p4,'OuterPosition',[1 210 470 330]);
% hold on
% plot(thetar_deg,fluxdq);
% legend('\phi_{d}','\phi_{q}','\phi_{0}','Location','NorthEast')
% xlabel(['\Theta_{r}', '(degree)'])
% ylabel('\phi_{dq0} (Wb)')
% % xlim([0.5 1])
% % ylim([499.5 500.5])
% grid on
% box on
% 
% p5 = figure;
% set(p5,'OuterPosition',[1 210 470 330]);
% hold on
% plot(thetar_deg,fluxabc);
% legend('\phi_{a}','\phi_{b}','\phi_{c}','Location','NorthEast')
% xlabel(['\Theta_{r}', '(degree)'])
% ylabel('\phi_{abc} (Wb)')
% % xlim([0.5 1])
% % ylim([499.5 500.5])
% grid on
% box on
% 
% p6 = figure;
% set(p6,'OuterPosition',[1 210 470 330]);
% hold on
% plot(thetar_deg,torque_e);
% legend('T','Location','NorthEast')
% xlabel(['\Theta_{r}', '(degree)'])
% ylabel('Torque (N·m)')
% % xlim([0.5 1])
% % ylim([499.5 500.5])
% grid on
% box on

%% Initial state of the simulator
% Speed
omega_m = 7200*2*pi/60;   % rad/s
omega_r = omega_m*N;    % rad/s
% Phase currents
% var_sim_iabc = [42.68 9.167 -51.847]';
var_sim_iabc = [0 0 0]';
% thetar
thetar_ant = 0;
thetar = 0;

% Total time 
param_exp_Ttotal = 0.01;
% paso de integración
param_sim_h = 1e-6; % (s)
% número de pasos a realizar para simular el experimento
param_sim_Nh = ceil(param_exp_Ttotal/param_sim_h);
% Time
tray_sim_tiempo = [0:param_sim_h:param_exp_Ttotal]';

% Phase voltages
Amp_vref = 50;
param_sim_fe = omega_r; 
var_sim_omega_vref = param_sim_fe*tray_sim_tiempo;
tray_sim_vfase = [Amp_vref*sin(var_sim_omega_vref),...
                  Amp_vref*sin(var_sim_omega_vref - 2*pi/3),...
                  Amp_vref*sin(var_sim_omega_vref + 2*pi/3)];


%% Recording variables
tray_sim_time       = zeros( param_sim_Nh, 1);
tray_sim_ifase      = zeros( param_sim_Nh, 3); 
tray_sim_idq0       = zeros( param_sim_Nh, 3);
% tray_sim_speed      = zeros( param_sim_Nh, 1); 
tray_sim_torque     = zeros( param_sim_Nh, 1);
tray_sim_vabc       = zeros( param_sim_Nh, 3);
tray_sim_thetar     = zeros( param_sim_Nh, 1);
tray_sim_thetar_deg = zeros( param_sim_Nh, 1);
tray_sim_fabc       = zeros( param_sim_Nh, 3);
tray_sim_fdq0       = zeros( param_sim_Nh, 3);
tray_sim_fdq0_lut   = zeros( param_sim_Nh, 3);

%% Simulation of derivative-based PMSM model
for h=1:param_sim_Nh
    
    theta = thetar*N;
    % Clarke-Park matrix
    T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
             sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
             1/2           1/2                    1/2];
    invT = inv(T);
    
    % Rotor angle for the LUTs
    thetar_deg = thetar*180/pi;

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
    flux_a = interpn(iaVec,ibVec,icVec,angleVec,fluxA2,i_a,i_b,i_c,thetar_deg);
    flux_b = interpn(iaVec,ibVec,icVec,angleVec,fluxB2,i_a,i_b,i_c,thetar_deg);
    flux_c = interpn(iaVec,ibVec,icVec,angleVec,fluxC2,i_a,i_b,i_c,thetar_deg);
    var_sim_fluxabc = [flux_a; flux_b; flux_c];

    % Actual dqz flux
    var_sim_fluxdq0 = T*var_sim_fluxabc;
    flux_d = var_sim_fluxdq0(1);
    flux_q = var_sim_fluxdq0(2);
    flux_0 = var_sim_fluxdq0(3);
    
    % Actual dqz flux from LUTs
    flux_d_lut = interp3(iqVec,idVec',angleVec,fluxD,i_q,i_d,thetar_deg);
    flux_q_lut = interp3(iqVec,idVec',angleVec,fluxQ,i_q,i_d,thetar_deg);
    flux_0_lut = interp3(iqVec,idVec',angleVec,flux0,i_q,i_d,thetar_deg);
    var_sim_fluxdq0_lut = [flux_d_lut;flux_q_lut;flux_0_lut];

    % Actual torque
    var_sim_torque = interp3(iqVec,idVec',angleVec,torque,i_q,i_d,thetar_deg);
    
    % Actual stator voltages
%     v_a = i_a*Rload;
%     v_b = i_b*Rload;
%     v_c = i_c*Rload;
    v_a = tray_sim_vfase(h,1);
    v_b = tray_sim_vfase(h,2);
    v_c = tray_sim_vfase(h,3);

    % Flux derivatives
    l_aa = interpn(iaVec,ibVec,icVec,angleVec,dFAdA,i_a,i_b,i_c,thetar_deg);
    l_ab = interpn(iaVec,ibVec,icVec,angleVec,dFAdB,i_a,i_b,i_c,thetar_deg);
    l_ac = interpn(iaVec,ibVec,icVec,angleVec,dFAdC,i_a,i_b,i_c,thetar_deg);
    phi_ar = interpn(iaVec,ibVec,icVec,angleVec,dFAdX,i_a,i_b,i_c,thetar_deg);

    l_ba = interpn(iaVec,ibVec,icVec,angleVec,dFBdA,i_a,i_b,i_c,thetar_deg);
    l_bb = interpn(iaVec,ibVec,icVec,angleVec,dFBdB,i_a,i_b,i_c,thetar_deg);
    l_bc = interpn(iaVec,ibVec,icVec,angleVec,dFBdC,i_a,i_b,i_c,thetar_deg);
    phi_br = interpn(iaVec,ibVec,icVec,angleVec,dFBdX,i_a,i_b,i_c,thetar_deg);

    l_ca = interpn(iaVec,ibVec,icVec,angleVec,dFCdA,i_a,i_b,i_c,thetar_deg);
    l_cb = interpn(iaVec,ibVec,icVec,angleVec,dFCdB,i_a,i_b,i_c,thetar_deg);
    l_cc = interpn(iaVec,ibVec,icVec,angleVec,dFCdC,i_a,i_b,i_c,thetar_deg);
    phi_cr = interpn(iaVec,ibVec,icVec,angleVec,dFCdX,i_a,i_b,i_c,thetar_deg);
    
    % System evolution
    di_a = (l_ab*l_bc*v_c - l_ab*l_cc*v_b - l_ac*l_bb*v_c + l_ac*l_cb*v_b + l_bb*l_cc*v_a - l_bc*l_cb*v_a - l_ab*l_bc*omega_r*phi_cr + l_ac*l_bb*omega_r*phi_cr + l_ab*l_cc*omega_r*phi_br - l_ac*l_cb*omega_r*phi_br - l_bb*l_cc*omega_r*phi_ar + l_bc*l_cb*omega_r*phi_ar - R_s*i_a*l_bb*l_cc + R_s*i_a*l_bc*l_cb + R_s*i_b*l_ab*l_cc - R_s*i_b*l_ac*l_cb - R_s*i_c*l_ab*l_bc + R_s*i_c*l_ac*l_bb)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
    di_b = -(l_aa*l_bc*v_c - l_aa*l_cc*v_b - l_ac*l_ba*v_c + l_ac*l_ca*v_b + l_ba*l_cc*v_a - l_bc*l_ca*v_a - l_aa*l_bc*omega_r*phi_cr + l_ac*l_ba*omega_r*phi_cr + l_aa*l_cc*omega_r*phi_br - l_ac*l_ca*omega_r*phi_br - l_ba*l_cc*omega_r*phi_ar + l_bc*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cc + R_s*i_a*l_bc*l_ca + R_s*i_b*l_aa*l_cc - R_s*i_b*l_ac*l_ca - R_s*i_c*l_aa*l_bc + R_s*i_c*l_ac*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
    di_c = (l_aa*l_bb*v_c - l_aa*l_cb*v_b - l_ab*l_ba*v_c + l_ab*l_ca*v_b + l_ba*l_cb*v_a - l_bb*l_ca*v_a - l_aa*l_bb*omega_r*phi_cr + l_ab*l_ba*omega_r*phi_cr + l_aa*l_cb*omega_r*phi_br - l_ab*l_ca*omega_r*phi_br - l_ba*l_cb*omega_r*phi_ar + l_bb*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cb + R_s*i_a*l_bb*l_ca + R_s*i_b*l_aa*l_cb - R_s*i_b*l_ab*l_ca - R_s*i_c*l_aa*l_bb + R_s*i_c*l_ab*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
    
    var_sim_iabc = var_sim_iabc + param_sim_h*[di_a;di_b;di_c];
       
    % Rotor angle estimation
    thetar = thetar_ant + omega_r*param_sim_h;
    if thetar > 2*pi/N
        thetar = thetar - 2*pi/N;
    elseif thetar < 0
        thetar = thetar + 2*pi/N;
    end
    thetar_ant = thetar;
         
    % Recordings
    tray_sim_time(h)        = h*param_sim_h;
    tray_sim_ifase(h,:)     = [i_a, i_b, i_c];
    tray_sim_idq0(h,:)      = var_sim_idqz';
%     tray_sim_speed(h,:)     = omega_m;
    tray_sim_torque(h,:)    = var_sim_torque;
    tray_sim_vabc(h,:)      = [v_a, v_b, v_c];
    tray_sim_thetar(h,:)    = thetar;
    tray_sim_thetar_deg(h,:)= thetar_deg;
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

p7 = figure;
set(p7,'OuterPosition',[1 210 470 330]);
hold on
plot(tray_sim_time,tray_sim_fdq0);
plot(tray_sim_time,tray_sim_fdq0_lut,':');
legend('\phi_{d}','\phi_{q}','\phi_{z}','\phi_{dLUT}','\phi_{dLUT}','\phi_{dLUT}','Location','NorthEast')
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

