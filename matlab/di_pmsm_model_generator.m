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

%% Obtention of new LUT: fluxA, fluxB, and fluxC as functions of iA, iB, iC, and thetar

% Extension of the fluxD, fluxQ, flux0, and torque LUTs to a maximum thetar of 360/N
% 30 (360/3/N) --> 90 (360/N)
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

% Calulation of phase fluxes LUTs as funtions of id, iq, and thetar
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

            % Phase currents
            ifase = invT(:,1:2)*[idVec(d);iqVec(q)];
            iA(d,q,x) = ifase(1);
            iB(d,q,x) = ifase(2);
            iC(d,q,x) = ifase(3);
        end
    end
end

ifa(1,1:91)=iA(11,15,:);
ifb(1,1:91)=iB(11,15,:);
ifc(1,1:91)=iC(11,15,:);
figure
plot(angleVec,ifa,angleVec,ifb,angleVec,ifc)
grid on

fd(1,1:91) = fluxD(11,15,:);
fq(1,1:91) = fluxQ(11,15,:);
figure
plot(angleVec,fd,angleVec,fq)
legend('\phi_{d}','\phi_{q}')
grid on

fa(1,1:91) = fluxA(11,15,:);
fb(1,1:91) = fluxB(11,15,:);
fc(1,1:91) = fluxC(11,15,:);
figure
plot(angleVec,fa,angleVec,fb,angleVec,fc)
legend('\phi_{a}','\phi_{b}','\phi_{c}')
grid on

%%
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

method = 'linear';

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
%                 T = 2/3*[ cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
%                          -sin(theta)   -sin(theta - 2*pi/3)   -sin(theta + 2*pi/3)
%                           1/2           1/2                    1/2];
                % d-q currents
                idq = T*[iaVec(a);ibVec(b);icVec(c)];
                % Interpolation form the phase fluxes LTUs obtained before
                fluxA2(a,b,c,x) = interpn(idVec,iqVec',angleVec,fluxA,idq(1),idq(2),angleVec(x),method);
                fluxB2(a,b,c,x) = interpn(idVec,iqVec',angleVec,fluxB,idq(1),idq(2),angleVec(x),method);
                fluxC2(a,b,c,x) = interpn(idVec,iqVec',angleVec,fluxC,idq(1),idq(2),angleVec(x),method);
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

%% Comprobamos que los flujos son iguales
id = idVec(11);
iq = iqVec(15);

for i=1:length(angleVec)
%     % Flujos extraídos directamente de la LUT sin interpolar
%     fd(i) = fluxD(11,15,i);
%     fq(i) = fluxQ(11,15,i);
%     f0(i) = flux0(11,15,i);
    % Flujos de fase sacados de fluxD y fluxQ
    flujod(i) = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(i));
    flujoq(i) = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(i));
    flujo0(i) = interpn(idVec,iqVec',angleVec,flux0,id,iq,angleVec(i));
    % Pasamos a flujo abc
    theta = N*angleVec(i)*pi/180;
    T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
             sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
             1/2           1/2                    1/2];
%     T = 2/3*[ cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
%              -sin(theta)   -sin(theta - 2*pi/3)   -sin(theta + 2*pi/3)
%               1/2           1/2                    1/2];
    invT = inv(T);
    
%     fabc = invT*[fd(i);fq(i);f0(i)];
%     fa(i) = fabc(1); 
%     fb(i) = fabc(2); 
%     fc(i) = fabc(3);

    flujosabc = invT*[flujod(i);flujoq(i);flujo0(i)];
    flujoa(i) = flujosabc(1); 
    flujob(i) = flujosabc(2); 
    flujoc(i) = flujosabc(3);
    
    % Flujos de fase sacados de fluxA y fluxB y fluxC
%     flujoa1_1(i) = fluxA(11,15,i);
%     flujob1_1(i) = fluxB(11,15,i);
%     flujoc1_1(i) = fluxC(11,15,i);
    flujoa1(i) = interpn(idVec,iqVec',angleVec,fluxA,id,iq,angleVec(i));
    flujob1(i) = interpn(idVec,iqVec',angleVec,fluxB,id,iq,angleVec(i));
    flujoc1(i) = interpn(idVec,iqVec',angleVec,fluxC,id,iq,angleVec(i));

    % Flujos de fase sacados de fluxA2 y fluxB2 y fluxC2
    iabc = invT(:,1:2)*[id;iq];
    flujoa2(i) = interpn(iaVec,ibVec,icVec,angleVec,fluxA2,iabc(1),iabc(2),iabc(3),angleVec(i));
    flujob2(i) = interpn(iaVec,ibVec,icVec,angleVec,fluxB2,iabc(1),iabc(2),iabc(3),angleVec(i));
    flujoc2(i) = interpn(iaVec,ibVec,icVec,angleVec,fluxC2,iabc(1),iabc(2),iabc(3),angleVec(i));
end

% Flujos de fase 
figure
% plot(angleVec,fa,'b',angleVec,fb,'r',angleVec,fc,'g')
hold on
plot(angleVec,flujoa,angleVec,flujob,angleVec,flujoc)
% figure
% plot(angleVec,flujoa1_1,'b',angleVec,flujob1_1,'r',angleVec,flujoc1_1,'g')
% hold on
plot(angleVec,flujoa1,'--',angleVec,flujob1,'--',angleVec,flujoc1,'--')
plot(angleVec,flujoa2,':',angleVec,flujob2,':',angleVec,flujoc2,':')
grid on

%% Flux partial derivatives calculation
% A phase
[dFAdA,dFAdB,dFAdC,dFAdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxA2);
% B phase
[dFBdA,dFBdB,dFBdC,dFBdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxB2);
% C phase
[dFCdA,dFCdB,dFCdC,dFCdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxC2);

%% Initial state of the simulator
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

% Comprobación de corrientes abc y de par
for i = 1:nd
    thetae = N*thetar_sim(i)*pi/180;
    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
              sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];
    invT = inv(T);
    iabc(:,i) = invT(:,1:2)*[i_d_sim(i);i_q_sim(i)];
    idq(:,i) = T*[i_a_sim(i);i_b_sim(i);i_c_sim(i)];
    torque_cal(i) = interpn(idVec,iqVec',angleVec,torque,idq(1,i),idq(2,i),thetar_sim_nuevo(i),"makima");
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


% %% FORMA 2: suponiendo iabc conocido para todo t
% var_sim_iabc = [i_a_sim';i_b_sim';i_c_sim'];
% var_sim_vabc = zeros(3,nd);
% % var_sim_fluxabc_lut_sim = zeros(3,nd);
% % var_sim_fluxabc_lut = zeros(3,nd);
% % var_sim_fluxdq0 = zeros(3,nd);
% % var_sim_fluxdq0_lut = zeros(3,nd);
% % var_sim_fluxdq0_lut_sim = zeros(3,nd);
% % var_sim_idqz = zeros(3,nd);
% 
% % l_aa = zeros(1,nd);
% % l_ab = zeros(1,nd);
% % l_ac = zeros(1,nd);
% % l_ba = zeros(1,nd);
% % l_bb = zeros(1,nd);
% % l_bc = zeros(1,nd);
% % l_ca = zeros(1,nd);
% % l_cb = zeros(1,nd);
% % l_cc = zeros(1,nd);
% % phi_ar = zeros(1,nd);
% % phi_br = zeros(1,nd);
% % phi_cr = zeros(1,nd);
% 
% divi = zeros(1,nd);
% detMatrix = zeros(1,nd);
% 
% omega_r = RPM*2*pi/60;  % [rad/s]
% 
% method = 'linear';
% 
% for i=2:nd
%     
%     thetar_deg = thetar_sim_nuevo(i-1);
% 
%     i_a = var_sim_iabc(1,i-1);
%     i_b = var_sim_iabc(2,i-1);
%     i_c = var_sim_iabc(3,i-1);
%     
% %     % COMPROBACIÓN DE FLUJOS
% %     % Flujos dq de LUT con idq recogidas de la simulación
% %     flux_d_lut_sim = interpn(idVec,iqVec',angleVec,fluxD,i_d_sim(i-1),i_q_sim(i-1),thetar_deg);
% %     flux_q_lut_sim = interpn(idVec,iqVec',angleVec,fluxQ,i_d_sim(i-1),i_q_sim(i-1),thetar_deg);
% %     flux_0_lut_sim = interpn(idVec,iqVec',angleVec,flux0,i_d_sim(i-1),i_q_sim(i-1),thetar_deg);
% %     var_sim_fluxdq0_lut_sim(:,i-1) = [flux_d_lut_sim;flux_q_lut_sim;flux_0_lut_sim];
% % 
% %     % Flujos dq de LUT con idq calculadas a partir de iabc
% %     thetae = N*thetar_sim(i-1)*pi/180;
% %     T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
% %              -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
% %               1/2            1/2                     1/2];
% %     var_sim_idqz(:,i-1) = T*var_sim_iabc(:,i-1);
% %     i_d = var_sim_idqz(1,i-1);
% %     i_q = var_sim_idqz(2,i-1);
% %     flux_d_lut = interpn(idVec,iqVec',angleVec,fluxD,i_d,i_q,thetar_deg);
% %     flux_q_lut = interpn(idVec,iqVec',angleVec,fluxQ,i_d,i_q,thetar_deg);
% %     flux_0_lut = interpn(idVec,iqVec',angleVec,flux0,i_d,i_q,thetar_deg);
% %     var_sim_fluxdq0_lut(:,i-1) = [flux_d_lut;flux_q_lut;flux_0_lut];
% % 
% %     % Flujos abc de LUT con iabc recogidas de la simulación
% %     flux_a_lut_sim = interpn(iaVec,ibVec,icVec,angleVec,fluxA2,i_a,i_b,i_c,thetar_deg);
% %     flux_b_lut_sim = interpn(iaVec,ibVec,icVec,angleVec,fluxB2,i_a,i_b,i_c,thetar_deg);
% %     flux_c_lut_sim = interpn(iaVec,ibVec,icVec,angleVec,fluxC2,i_a,i_b,i_c,thetar_deg);
% %     var_sim_fluxabc_lut_sim(:,i-1) = [flux_a_lut_sim; flux_b_lut_sim; flux_c_lut_sim];
% % 
% %     % Flujos abc de LUT con idq calculadas a partir de idq
% %     var_sim_iabc_cal = inv(T)*[i_d_sim(i-1);i_q_sim(i-1);0];
% %     flux_a_lut = interpn(iaVec,ibVec,icVec,angleVec,fluxA2,var_sim_iabc_cal(1),var_sim_iabc_cal(2),var_sim_iabc_cal(3),thetar_deg);
% %     flux_b_lut = interpn(iaVec,ibVec,icVec,angleVec,fluxB2,var_sim_iabc_cal(1),var_sim_iabc_cal(2),var_sim_iabc_cal(3),thetar_deg);
% %     flux_c_lut = interpn(iaVec,ibVec,icVec,angleVec,fluxC2,var_sim_iabc_cal(1),var_sim_iabc_cal(2),var_sim_iabc_cal(3),thetar_deg);
% %     var_sim_fluxabc_lut(:,i-1) = [flux_a_lut; flux_b_lut; flux_c_lut];
% % 
% %     % COMPROBACIÓN DE DERIVADAS DE FLUJOS
% %     l_aa(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFAdA,i_a,i_b,i_c,thetar_deg,method);
% %     l_ab(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFAdB,i_a,i_b,i_c,thetar_deg,method);
% %     l_ac(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFAdC,i_a,i_b,i_c,thetar_deg,method);
% %     phi_ar(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFAdX,i_a,i_b,i_c,thetar_deg,method);
% %     
% %     l_ba(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFBdA,i_a,i_b,i_c,thetar_deg,method);
% %     l_bb(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFBdB,i_a,i_b,i_c,thetar_deg,method);
% %     l_bc(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFBdC,i_a,i_b,i_c,thetar_deg,method);
% %     phi_br(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFBdX,i_a,i_b,i_c,thetar_deg,method);
% %     
% %     l_ca(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFCdA,i_a,i_b,i_c,thetar_deg,method);
% %     l_cb(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFCdB,i_a,i_b,i_c,thetar_deg,method);
% %     l_cc(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFCdC,i_a,i_b,i_c,thetar_deg,method);
% %     phi_cr(i-1) = interpn(iaVec,ibVec,icVec,angleVec,dFCdX,i_a,i_b,i_c,thetar_deg,method);    
% 
% 
%     % Flux derivatives
%     l_aa = interpn(iaVec,ibVec,icVec,angleVec,dFAdA,i_a,i_b,i_c,thetar_deg,method);
%     l_ab = interpn(iaVec,ibVec,icVec,angleVec,dFAdB,i_a,i_b,i_c,thetar_deg,method);
%     l_ac = interpn(iaVec,ibVec,icVec,angleVec,dFAdC,i_a,i_b,i_c,thetar_deg,method);
%     phi_ar = interpn(iaVec,ibVec,icVec,angleVec,dFAdX,i_a,i_b,i_c,thetar_deg,method);
%     
%     l_ba = interpn(iaVec,ibVec,icVec,angleVec,dFBdA,i_a,i_b,i_c,thetar_deg,method);
%     l_bb = interpn(iaVec,ibVec,icVec,angleVec,dFBdB,i_a,i_b,i_c,thetar_deg,method);
%     l_bc = interpn(iaVec,ibVec,icVec,angleVec,dFBdC,i_a,i_b,i_c,thetar_deg,method);
%     phi_br = interpn(iaVec,ibVec,icVec,angleVec,dFBdX,i_a,i_b,i_c,thetar_deg,method);
%     
%     l_ca = interpn(iaVec,ibVec,icVec,angleVec,dFCdA,i_a,i_b,i_c,thetar_deg,method);
%     l_cb = interpn(iaVec,ibVec,icVec,angleVec,dFCdB,i_a,i_b,i_c,thetar_deg,method);
%     l_cc = interpn(iaVec,ibVec,icVec,angleVec,dFCdC,i_a,i_b,i_c,thetar_deg,method);
%     phi_cr = interpn(iaVec,ibVec,icVec,angleVec,dFCdX,i_a,i_b,i_c,thetar_deg,method);
%     
%     % System evolution
%     param_sim_h = time(i)-time(i-1);
%     
%     i_a_p1 = var_sim_iabc(1,i);
%     i_b_p1 = var_sim_iabc(2,i);
%     i_c_p1 = var_sim_iabc(3,i);
% 
%     v_a =  l_aa*(i_a_p1-i_a)/param_sim_h + l_ab*(i_b_p1-i_b)/param_sim_h + l_ac*(i_c_p1-i_c)/param_sim_h + R_s * i_a + phi_ar * omega_r;    
%     v_b =  l_ba*(i_a_p1-i_a)/param_sim_h + l_bb*(i_b_p1-i_b)/param_sim_h + l_bc*(i_c_p1-i_c)/param_sim_h + R_s * i_b + phi_br * omega_r;
%     v_c =  l_ca*(i_a_p1-i_a)/param_sim_h + l_cb*(i_b_p1-i_b)/param_sim_h + l_cc*(i_c_p1-i_c)/param_sim_h + R_s * i_c + phi_cr * omega_r;
%     
%     var_sim_vabc(:,i-1) = [v_a;v_b;v_c];
% 
%     divi(i-1) = l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca;
%     Matrix_dF_di = [l_aa, l_ab, l_ac; l_ba, l_bb, l_bc; l_ca, l_cb, l_cc];
%     detMatrix(i-1) = det(Matrix_dF_di);
% 
% end
% 
% % figure
% % plot(time,i_d_sim,time,i_q_sim)
% % hold on
% % plot(time,var_sim_idqz(1:2,:),'--')
% % grid on
% % 
% % figure
% % plot(time,var_sim_fluxdq0_lut)
% % hold on
% % plot(time,var_sim_fluxdq0_lut_sim,'--')
% % grid on
% % 
% % figure
% % plot(time,var_sim_fluxabc_lut)
% % hold on
% % plot(time,var_sim_fluxabc_lut_sim,'--')
% % grid on
% 
% 
% % figure
% % plot(time,l_aa,time,l_bb,time,l_cc)
% % legend('l_{aa}','l_{bb}','l_{cc}')
% % grid on
% % 
% % figure
% % plot(time,l_aa,time,l_ab,time,l_ac)
% % legend('l_{aa}','l_{ab}','l_{ac}')
% % grid on
% % 
% % figure
% % plot(time,l_ba,time,l_bb,time,l_bc)
% % legend('l_{bc}','l_{bb}','l_{bc}')
% % grid on
% % 
% % figure
% % plot(time,l_ca,time,l_cb,time,l_cc)
% % legend('l_{ca}','l_{cb}','l_{cc}')
% % grid on
% % 
% % figure
% % plot(time,phi_ar,time,phi_br,time,phi_cr)
% % legend('d\phi_{a}/d\theta_r','d\phi_{b}/d\theta_r','d\phi_{c}/d\theta_r')
% % grid on
% % 
% % figure
% % plot(time,l_ab,time,l_ac,time,l_bc)
% % hold on
% % plot(time,l_ba,'--',time,l_ca,'--',time,l_cb,'--')
% % legend('l_{ab}','l_{ac}','l_{bc}','l_{ba}','l_{ca}','l_{cb}')
% % grid on
% 
% 
% figure
% plot(time,var_sim_vabc)
% hold on
% plot(time,v_a_sim,'--',time,v_b_sim,'--',time,v_c_sim,'--')
% ylim([-600,600])
% grid on


%% Intento de simulación discreta del sistema conocidas vabc
var_sim_iabc = zeros(3,nd);
var_sim_iabc(:,1) = [i_a_sim(1);i_b_sim(1);i_c_sim(1)];
divi = zeros(1,nd);
var_sim_vabc = zeros(3,nd);

omega_r = RPM*2*pi/60;  % [rad/s]

method = 'linear';

for i = 2:nd
    
    thetar_deg = thetar_sim_nuevo(i-1);
    
    v_a = v_a_sim(i-1);
    v_b = v_b_sim(i-1);
    v_c = v_c_sim(i-1);

    i_a = var_sim_iabc(1,i-1);
    i_b = var_sim_iabc(2,i-1);
    i_c = var_sim_iabc(3,i-1);

    % Flux derivatives
    l_aa = interpn(iaVec,ibVec,icVec,angleVec,dFAdA,i_a,i_b,i_c,thetar_deg,method);
    l_ab = interpn(iaVec,ibVec,icVec,angleVec,dFAdB,i_a,i_b,i_c,thetar_deg,method);
    l_ac = interpn(iaVec,ibVec,icVec,angleVec,dFAdC,i_a,i_b,i_c,thetar_deg,method);
    phi_ar = interpn(iaVec,ibVec,icVec,angleVec,dFAdX,i_a,i_b,i_c,thetar_deg,method);
    
    l_ba = interpn(iaVec,ibVec,icVec,angleVec,dFBdA,i_a,i_b,i_c,thetar_deg,method);
    l_bb = interpn(iaVec,ibVec,icVec,angleVec,dFBdB,i_a,i_b,i_c,thetar_deg,method);
    l_bc = interpn(iaVec,ibVec,icVec,angleVec,dFBdC,i_a,i_b,i_c,thetar_deg,method);
    phi_br = interpn(iaVec,ibVec,icVec,angleVec,dFBdX,i_a,i_b,i_c,thetar_deg,method);
    
    l_ca = interpn(iaVec,ibVec,icVec,angleVec,dFCdA,i_a,i_b,i_c,thetar_deg,method);
    l_cb = interpn(iaVec,ibVec,icVec,angleVec,dFCdB,i_a,i_b,i_c,thetar_deg,method);
    l_cc = interpn(iaVec,ibVec,icVec,angleVec,dFCdC,i_a,i_b,i_c,thetar_deg,method);
    phi_cr = interpn(iaVec,ibVec,icVec,angleVec,dFCdX,i_a,i_b,i_c,thetar_deg,method);
    
    
%     % System evolution
%     divi(i) = l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca;
%     di_a = (l_ab*l_bc*v_c - l_ab*l_cc*v_b - l_ac*l_bb*v_c + l_ac*l_cb*v_b + l_bb*l_cc*v_a - l_bc*l_cb*v_a - l_ab*l_bc*omega_r*phi_cr + l_ac*l_bb*omega_r*phi_cr + l_ab*l_cc*omega_r*phi_br - l_ac*l_cb*omega_r*phi_br - l_bb*l_cc*omega_r*phi_ar + l_bc*l_cb*omega_r*phi_ar - R_s*i_a*l_bb*l_cc + R_s*i_a*l_bc*l_cb + R_s*i_b*l_ab*l_cc - R_s*i_b*l_ac*l_cb - R_s*i_c*l_ab*l_bc + R_s*i_c*l_ac*l_bb)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
%     di_b = -(l_aa*l_bc*v_c - l_aa*l_cc*v_b - l_ac*l_ba*v_c + l_ac*l_ca*v_b + l_ba*l_cc*v_a - l_bc*l_ca*v_a - l_aa*l_bc*omega_r*phi_cr + l_ac*l_ba*omega_r*phi_cr + l_aa*l_cc*omega_r*phi_br - l_ac*l_ca*omega_r*phi_br - l_ba*l_cc*omega_r*phi_ar + l_bc*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cc + R_s*i_a*l_bc*l_ca + R_s*i_b*l_aa*l_cc - R_s*i_b*l_ac*l_ca - R_s*i_c*l_aa*l_bc + R_s*i_c*l_ac*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
%     di_c = (l_aa*l_bb*v_c - l_aa*l_cb*v_b - l_ab*l_ba*v_c + l_ab*l_ca*v_b + l_ba*l_cb*v_a - l_bb*l_ca*v_a - l_aa*l_bb*omega_r*phi_cr + l_ab*l_ba*omega_r*phi_cr + l_aa*l_cb*omega_r*phi_br - l_ab*l_ca*omega_r*phi_br - l_ba*l_cb*omega_r*phi_ar + l_bb*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cb + R_s*i_a*l_bb*l_ca + R_s*i_b*l_aa*l_cb - R_s*i_b*l_ab*l_ca - R_s*i_c*l_aa*l_bb + R_s*i_c*l_ab*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
%     
%     param_sim_h = time(i)-time(i-1);
%     var_sim_iabc(:,i) = var_sim_iabc(:,i-1) + param_sim_h*[di_a;di_b;di_c];
    
    % System evolution (MATRICIAL)
    Matrix_dF_di = [l_aa, l_ab, l_ac; l_ba, l_bb, l_bc; l_ca, l_cb, l_cc];
    Matrix_dF_dx = [phi_ar; phi_br; phi_cr];
    
    var_sim_vabc(:,i-1) = [v_a;v_b;v_c];
    param_sim_h = time(i)-time(i-1);
    di_dt = Matrix_dF_di\(var_sim_vabc(:,i-1) - R_s*var_sim_iabc(:,i-1) - omega_r*Matrix_dF_dx);
    var_sim_iabc(:,i) = var_sim_iabc(:,i-1) + param_sim_h*di_dt; 
    
end

figure
plot(time,var_sim_iabc)
hold on
plot(time,i_a_sim,'--',time,i_b_sim,'--',time,i_c_sim,'--')

% CONCLUSIÓN: No converje porde divi es muy próximo a cero.

% %% Simulation of derivative-based PMSM model
% 
% for h=1:param_sim_Nh
% 
%     % Clarke-Park matrix
%     T = 2/3*[cos(theta)    cos(theta - 2*pi/3)    cos(theta + 2*pi/3);
%              sin(theta)    sin(theta - 2*pi/3)    sin(theta + 2*pi/3)
%              1/2           1/2                    1/2];
%     invT = inv(T);
%     
%     % Rotor angle for the LUTs
%     thetar = theta*180/pi/N;
% 
%     % Actual abc currents
%     i_a = var_sim_iabc(1);
%     i_b = var_sim_iabc(2);
%     i_c = var_sim_iabc(3);
% 
%     % Actual dq0 currents
%     var_sim_idqz = T*var_sim_iabc;
%     i_d = var_sim_idqz(1);
%     i_q = var_sim_idqz(2);
%     i_0 = var_sim_idqz(3);
%     
%     % Actual abc flux
%     flux_a = interpn(iaVec,ibVec,icVec,angleVec,fluxA2,i_a,i_b,i_c,thetar);
%     flux_b = interpn(iaVec,ibVec,icVec,angleVec,fluxB2,i_a,i_b,i_c,thetar);
%     flux_c = interpn(iaVec,ibVec,icVec,angleVec,fluxC2,i_a,i_b,i_c,thetar);
%     var_sim_fluxabc = [flux_a; flux_b; flux_c];
% 
%     % Actual dqz flux
%     var_sim_fluxdq0 = T*var_sim_fluxabc;
%     flux_d = var_sim_fluxdq0(1);
%     flux_q = var_sim_fluxdq0(2);
%     flux_0 = var_sim_fluxdq0(3);
%     
%     % Actual dqz flux from LUTs
%     flux_d_lut = interp3(iqVec,idVec',angleVec,fluxD,i_q,i_d,thetar);
%     flux_q_lut = interp3(iqVec,idVec',angleVec,fluxQ,i_q,i_d,thetar);
%     flux_0_lut = interp3(iqVec,idVec',angleVec,flux0,i_q,i_d,thetar);
%     var_sim_fluxdq0_lut = [flux_d_lut;flux_q_lut;flux_0_lut];
% 
%     % Actual torque
%     var_sim_torque = interp3(iqVec,idVec',angleVec,torque,i_q,i_d,thetar);
%     
%     % Actual stator voltages
%     v_a = i_a*Rload;
%     v_b = i_b*Rload;
%     v_c = i_c*Rload;
% 
%     % Flux derivatives
%     l_aa = interpn(iaVec,ibVec,icVec,angleVec,dFAdA,i_a,i_b,i_c,thetar);
%     l_ab = interpn(iaVec,ibVec,icVec,angleVec,dFAdB,i_a,i_b,i_c,thetar);
%     l_ac = interpn(iaVec,ibVec,icVec,angleVec,dFAdC,i_a,i_b,i_c,thetar);
%     phi_ar = interpn(iaVec,ibVec,icVec,angleVec,dFAdX,i_a,i_b,i_c,thetar);
% 
%     l_ba = interpn(iaVec,ibVec,icVec,angleVec,dFBdA,i_a,i_b,i_c,thetar);
%     l_bb = interpn(iaVec,ibVec,icVec,angleVec,dFBdB,i_a,i_b,i_c,thetar);
%     l_bc = interpn(iaVec,ibVec,icVec,angleVec,dFBdC,i_a,i_b,i_c,thetar);
%     phi_br = interpn(iaVec,ibVec,icVec,angleVec,dFBdX,i_a,i_b,i_c,thetar);
% 
%     l_ca = interpn(iaVec,ibVec,icVec,angleVec,dFCdA,i_a,i_b,i_c,thetar);
%     l_cb = interpn(iaVec,ibVec,icVec,angleVec,dFCdB,i_a,i_b,i_c,thetar);
%     l_cc = interpn(iaVec,ibVec,icVec,angleVec,dFCdC,i_a,i_b,i_c,thetar);
%     phi_cr = interpn(iaVec,ibVec,icVec,angleVec,dFCdX,i_a,i_b,i_c,thetar);
%     
%     % System evolution
%     di_a = (l_ab*l_bc*v_c - l_ab*l_cc*v_b - l_ac*l_bb*v_c + l_ac*l_cb*v_b + l_bb*l_cc*v_a - l_bc*l_cb*v_a - l_ab*l_bc*omega_r*phi_cr + l_ac*l_bb*omega_r*phi_cr + l_ab*l_cc*omega_r*phi_br - l_ac*l_cb*omega_r*phi_br - l_bb*l_cc*omega_r*phi_ar + l_bc*l_cb*omega_r*phi_ar - R_s*i_a*l_bb*l_cc + R_s*i_a*l_bc*l_cb + R_s*i_b*l_ab*l_cc - R_s*i_b*l_ac*l_cb - R_s*i_c*l_ab*l_bc + R_s*i_c*l_ac*l_bb)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
%     di_b = -(l_aa*l_bc*v_c - l_aa*l_cc*v_b - l_ac*l_ba*v_c + l_ac*l_ca*v_b + l_ba*l_cc*v_a - l_bc*l_ca*v_a - l_aa*l_bc*omega_r*phi_cr + l_ac*l_ba*omega_r*phi_cr + l_aa*l_cc*omega_r*phi_br - l_ac*l_ca*omega_r*phi_br - l_ba*l_cc*omega_r*phi_ar + l_bc*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cc + R_s*i_a*l_bc*l_ca + R_s*i_b*l_aa*l_cc - R_s*i_b*l_ac*l_ca - R_s*i_c*l_aa*l_bc + R_s*i_c*l_ac*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
%     di_c = (l_aa*l_bb*v_c - l_aa*l_cb*v_b - l_ab*l_ba*v_c + l_ab*l_ca*v_b + l_ba*l_cb*v_a - l_bb*l_ca*v_a - l_aa*l_bb*omega_r*phi_cr + l_ab*l_ba*omega_r*phi_cr + l_aa*l_cb*omega_r*phi_br - l_ab*l_ca*omega_r*phi_br - l_ba*l_cb*omega_r*phi_ar + l_bb*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cb + R_s*i_a*l_bb*l_ca + R_s*i_b*l_aa*l_cb - R_s*i_b*l_ab*l_ca - R_s*i_c*l_aa*l_bb + R_s*i_c*l_ab*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
%     
%     var_sim_iabc = var_sim_iabc + param_sim_h*[di_a;di_b;di_c];
%        
%     % Rotor angle estimation
%     theta = theta_ant + omega_m*param_sim_h;
%     if theta > 2*pi
%         theta = theta - 2*pi;
%     elseif theta < 0
%         theta = theta + 2*pi;
%     end
%     theta_ant = theta;
%          
%     % Recordings
%     tray_sim_time(h)        = h*param_sim_h;
%     tray_sim_ifase(h,:)     = [i_a, i_b, i_c];
%     tray_sim_idq0(h,:)      = var_sim_idqz';
% %     tray_sim_speed(h,:)     = omega_m;
%     tray_sim_torque(h,:)    = var_sim_torque;
%     tray_sim_vabc(h,:)      = [v_a, v_b, v_c];
%     tray_sim_thetar(h,:)    = thetar;
%     tray_sim_fabc(h,:)      = var_sim_fluxabc;
%     tray_sim_fdq0(h,:)      = var_sim_fluxdq0;
%     tray_sim_fdq0_lut(h,:)  = var_sim_fluxdq0_lut;
% end
% 
% %% Gráficas para guardar
% close all
% 
% p1 = figure;
% set(p1,'OuterPosition',[1 210 470 330]);
% hold on
% plot(tray_sim_time,tray_sim_ifase);
% fase=legend('\it{i_{a}}','\it{i_{b}}','\it{i_{c}}','Orientation','Horizontal');
% set(fase,'FontSize',10,'FontName','Times New Roman')
% xlabel('Time (s)','FontSize',10)
% ylabel('Current (A)','FontSize',10)
% % xlim([0.9 1])
% % ylim([-2 2.5])
% grid on
% box on
% 
% p2 = figure;
% set(p2,'OuterPosition',[1 210 470 330]);
% hold on
% plot(tray_sim_time,tray_sim_idq0);
% alfabeta=legend('\it{i_{d}}','\it{i_{q}}','\it{i_{0}}','Orientation','Horizontal');
% set(alfabeta,'FontSize',10,'FontName','Times New Roman')
% xlabel('Time (s)','FontSize',10)
% ylabel('Current (A)','FontSize',10)
% % xlim([0.9 1])
% % ylim([-2 2.5])
% grid on
% box on
% 
% p4 = figure;
% set(p4,'OuterPosition',[1 210 470 330]);
% hold on
% plot(tray_sim_time,tray_sim_vabc);
% dd=legend('\it{v_{a}}','\it{v_{b}}','\it{v_{c}}','Orientation','Horizontal');
% set(dd,'FontSize',10,'FontName','Times New Roman')
% xlabel('Time (s)','FontSize',10)
% ylabel('Voltage (V)','FontSize',10)
% % xlim([0.5 1])
% % ylim([0 2.5])
% grid on
% box on
% 
% p6 = figure;
% set(p6,'OuterPosition',[1 210 470 330]);
% hold on
% plot(tray_sim_time,tray_sim_torque);
% legend('T','Location','NorthEast')
% xlabel('Time (s)')
% ylabel('Torque (N·m)')
% % xlim([0.5 1])
% % ylim([499.5 500.5])
% grid on
% box on
% 
% 
% % %% Prueba de LUT
% % 
% % for i=1:length(angleVec)
% %     flujo(i) = interp3(iqVec,idVec',angleVec,fluxQ,iqVec(11),idVec(11),angleVec(i));
% %     flujo_p(i) = fluxQ(11,11,i);
% % end
% % figure
% % plot(angleVec,flujo)
% % hold on
% % plot(angleVec,flujo_p,'r')
% % 
% % 
% % for i=1:length(angleVec)
% %     flujo(i)=interp3(iqVec,idVec',angleVec,fluxD,iqVec(11),idVec(11),angleVec(i));
% %     flujo_p(i) = fluxD(11,11,i);
% % end
% % figure
% % plot(angleVec,flujo)
% % hold on
% % plot(angleVec,flujo_p,'r')

