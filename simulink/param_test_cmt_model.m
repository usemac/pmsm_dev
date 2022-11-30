
clear all;
clc;

%% Motor parameters
ee_ece_table
R_s = 0.07;

Rload = 10;     % [Ohm]
RPM = 7200;     % [rpm] mechanical speed wr

%% Ejecución del modelo de Simulink
sim('params_test_cmt')

%% Load data
time = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.angular_velocity.series.time;
thetar_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.angular_position.series.values;     % [deg] rotor angle, from wr
omegar_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.angular_velocity.series.values;     % [rpm] rotor mechanical speed wr
torque_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.electrical_torque.series.values;
ia_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.i_a.series.values;
ib_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.i_b.series.values;
ic_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.i_c.series.values;
id_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.i_d.series.values;
iq_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.i_q.series.values;
va_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.v_a.series.values;
vb_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.v_b.series.values;
vc_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM.v_c.series.values;

nd = length(time);

%% Extension of the fluxD, fluxQ, flux0, and torque LUTs to a maximum thetar of 360/N
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

%% Cambio de formato de thetar
thetar_aux = thetar_sim;
thetar_sim_nuevo = zeros(1,nd);
for i = 1:nd
    if thetar_aux(i) > 90
        thetar_aux = thetar_aux-90;
        thetar_sim_nuevo(i) = thetar_aux(i);
    else
        thetar_sim_nuevo(i) = thetar_aux(i);
    end
end

%% Tests
% Cálculo de thetar a partir de la velocidad de simulink
thetar = zeros(1,length(time));

thetar(1) = 0;
thetar_ant = thetar(1);

for i = 2:length(time)
    omegar = omegar_sim(i)*2*pi/60;
    thetar(i) = thetar_ant + omegar*(time(i)-time(i-1));
    thetar_ant = thetar(i);
end

figure
plot(time,thetar_sim)
hold on
plot(time,thetar*180/pi,'--')
grid on

% Comprobación de corrientes abc y de par
for i = 1:nd
    thetae = N*thetar_sim(i)*pi/180;
    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
              sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];
    invT = inv(T);
    iabc(:,i) = invT(:,1:2)*[id_sim(i);iq_sim(i)];
    idq(:,i) = T*[ia_sim(i);ib_sim(i);ic_sim(i)];
    torque_cal(i) = interp3(iqVec,idVec',angleVec,torque,iq_sim(i),id_sim(i),thetar_sim_nuevo(i),"linear");
end

figure
plot(time,ia_sim,time,ib_sim,time,ic_sim)
hold on
plot(time,iabc,'--')

figure
plot(time,id_sim,time,iq_sim)
hold on
plot(time,idq,'--')

figure
plot(time,torque_sim)
hold on
plot(time,torque_cal,'--')


%% Obtention of new LUT: fluxA, fluxB, and fluxC as functions of iA, iB, iC, and thetar

% Calulation of phase fluxes LTUs as funtions of id, iq, and thetar
fluxA = zeros(niD,niQ,nang);
fluxB = zeros(niD,niQ,nang);
fluxC = zeros(niD,niQ,nang);
for x=1:nang
    for d=1:niD
        for q=1:niQ
            % Clarke-Park matrix
            theta = N*angleVec(x)*pi/180;
%             T = 2/3*[ cos(theta)     cos(theta - 2*pi/3)     cos(theta + 2*pi/3);
%                      -sin(theta)    -sin(theta - 2*pi/3)    -sin(theta + 2*pi/3)
%                       1/2            1/2                     1/2];
            T = 2/3*[ cos(theta)     cos(theta - 2*pi/3)     cos(theta + 2*pi/3);
                      sin(theta)     sin(theta - 2*pi/3)     sin(theta + 2*pi/3)
                      1/2            1/2                     1/2];
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
% iaVec = [-210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210];
% ibVec = [-210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210];
% icVec = [-210 -180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180 210];
iaVec = [-350 -300 -250 -200 -150 -100 -50 -0 50 100 150 200 250 300 350];
ibVec = [-350 -300 -250 -200 -150 -100 -50 -0 50 100 150 200 250 300 350];
icVec = [-350 -300 -250 -200 -150 -100 -50 -0 50 100 150 200 250 300 350];

niA = length(iaVec);
niB = length(ibVec);
niC = length(icVec);

fluxA2 = zeros(niA,niB,niC,nang);
fluxB2 = zeros(niA,niB,niC,nang);
fluxC2 = zeros(niA,niB,niC,nang);
torque2 = zeros(niA,niB,niC,nang);

% Calculation of LTUs for the phase fluxes in terms of iA, iB, iC and
% thetar
for x=1:nang
    for a=1:niA
        for b=1:niB
            for c=1:niC
                % Clarke-Park matrix
                theta = N*angleVec(x)*pi/180;
%                 T = 2/3*[ cos(theta)     cos(theta - 2*pi/3)     cos(theta + 2*pi/3);
%                          -sin(theta)    -sin(theta - 2*pi/3)    -sin(theta + 2*pi/3)
%                           1/2            1/2                     1/2];
                T = 2/3*[ cos(theta)     cos(theta - 2*pi/3)     cos(theta + 2*pi/3);
                          sin(theta)     sin(theta - 2*pi/3)     sin(theta + 2*pi/3)
                          1/2            1/2                     1/2];
                invT = inv(T);
                % d-q currents
                idq = T*[iaVec(a);ibVec(b);icVec(c)];
                
                if idq(1) < idVec(1) || idq(1) > idVec(niD) || idq(2) < iqVec(1) || idq(2) > iqVec(niQ) || angleVec(x) < angleVec(1) || angleVec(x) > angleVec(nang)
                    method = 'makima';
                else
                    method = 'linear';
                end

                % Interpolation form the phase fluxes LTUs obtained before
                fluxA2(a,b,c,x) = interp3(iqVec,idVec',angleVec,fluxA,idq(2),idq(1),angleVec(x),method);
                fluxB2(a,b,c,x) = interp3(iqVec,idVec',angleVec,fluxB,idq(2),idq(1),angleVec(x),method);
                fluxC2(a,b,c,x) = interp3(iqVec,idVec',angleVec,fluxC,idq(2),idq(1),angleVec(x),method);
                torque2(a,b,c,x) = interp3(iqVec,idVec',angleVec,torque,idq(2),idq(1),angleVec(x),method);
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

% Flux partial derivatives calculation
% A phase
[dFAdA,dFAdB,dFAdC,dFAdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxA2);
% B phase
[dFBdA,dFBdB,dFBdC,dFBdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxB2);
% C phase
[dFCdA,dFCdB,dFCdC,dFCdX] = ee_calculateFluxPartialDerivatives(iaVec,ibVec,icVec,angleVec*pi/180,fluxC2);

%% Representation of the partial derivatives
% Operational point
id = idVec(11);
iq = iqVec(15);

% T = 2/3*[ cos(theta)     cos(theta - 2*pi/3)     cos(theta + 2*pi/3);
%          -sin(theta)    -sin(theta - 2*pi/3)    -sin(theta + 2*pi/3)
%           1/2            1/2                     1/2];
T = 2/3*[ cos(theta)     cos(theta - 2*pi/3)     cos(theta + 2*pi/3);
          sin(theta)     sin(theta - 2*pi/3)     sin(theta + 2*pi/3)
          1/2            1/2                     1/2];
invT = inv(T);
iabc = invT*[id;iq;0];

laa = zeros(1,length(angleVec));
lab = zeros(1,length(angleVec));
lac = zeros(1,length(angleVec));
lba = zeros(1,length(angleVec));
lbb = zeros(1,length(angleVec));
lbc = zeros(1,length(angleVec));
lca = zeros(1,length(angleVec));
lcb = zeros(1,length(angleVec));
lcc = zeros(1,length(angleVec));

for i = 1:length(angleVec)
    % Flux derivatives
    laa(i) = interpn(iaVec,ibVec,icVec,angleVec,dFAdA,iabc(1),iabc(2),iabc(3),angleVec(i));
    lab(i) = interpn(iaVec,ibVec,icVec,angleVec,dFAdB,iabc(1),iabc(2),iabc(3),angleVec(i));
    lac(i) = interpn(iaVec,ibVec,icVec,angleVec,dFAdC,iabc(1),iabc(2),iabc(3),angleVec(i));
    
    lba(i) = interpn(iaVec,ibVec,icVec,angleVec,dFBdA,iabc(1),iabc(2),iabc(3),angleVec(i));
    lbb(i) = interpn(iaVec,ibVec,icVec,angleVec,dFBdB,iabc(1),iabc(2),iabc(3),angleVec(i));
    lbc(i) = interpn(iaVec,ibVec,icVec,angleVec,dFBdC,iabc(1),iabc(2),iabc(3),angleVec(i));
    
    lca(i) = interpn(iaVec,ibVec,icVec,angleVec,dFCdA,iabc(1),iabc(2),iabc(3),angleVec(i));
    lcb(i) = interpn(iaVec,ibVec,icVec,angleVec,dFCdB,iabc(1),iabc(2),iabc(3),angleVec(i));
    lcc(i) = interpn(iaVec,ibVec,icVec,angleVec,dFCdC,iabc(1),iabc(2),iabc(3),angleVec(i));
end

% laa = dFA/diA; lbb = dFB/diB; lcc = dFC/diC
figure
plot(angleVec,laa,angleVec,lbb,angleVec,lcc)
grid on
legend('l_{aa}','l_{bb}','l_{cc}')

% lab = dFA/diB; lba = dFB/diA
figure
plot(angleVec,lab)
hold on
plot(angleVec,lba,'--')
legend('l_{ab}','l_{ba}')

% lac = dFA/diC; lca = dFC/diA
figure
plot(angleVec,lac)
hold on
plot(angleVec,lca,'--')
legend('l_{ac}','l_{ca}')

% lbc = dFB/diC; lcb = dFC/diB
figure
plot(angleVec,lbc)
hold on
plot(angleVec,lcb,'--')
legend('l_{bc}','l_{cb}')

% lab = dFA/diB; lbc = dFB/diC; lca = dFC/diA
figure
plot(angleVec,lab,angleVec,lbc,angleVec,lca)
grid on
legend('l_{ab}','l_{bc}','l_{ca}')

%% Intento de simulación del sistema conocidas vabc

var_sim_iabc = zeros(3,nd);
var_sim_iabc(:,1) = [ia_sim(1);ib_sim(1);ic_sim(1)];

omega_r = RPM*2*pi/60;

for i=2:nd
    thetar_deg = thetar_sim_nuevo(i-1);

    v_a = va_sim(i-1);
    v_b = vb_sim(i-1);
    v_c = vc_sim(i-1);

    i_a = var_sim_iabc(1,i-1);
    i_b = var_sim_iabc(2,i-1);
    i_c = var_sim_iabc(3,i-1);

    if i_a < iaVec(1) || i_a > iaVec(length(iaVec)) || i_b < ibVec(1) || i_b > ibVec(length(ibVec)) || i_c < icVec(1) || i_c > icVec(length(icVec)) || thetar_deg < angleVec(1) || thetar_deg > angleVec(nang)
        method = 'makima';
    else
        method = 'linear';
    end

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
    
    % System evolution
    di_a = (l_ab*l_bc*v_c - l_ab*l_cc*v_b - l_ac*l_bb*v_c + l_ac*l_cb*v_b + l_bb*l_cc*v_a - l_bc*l_cb*v_a - l_ab*l_bc*omega_r*phi_cr + l_ac*l_bb*omega_r*phi_cr + l_ab*l_cc*omega_r*phi_br - l_ac*l_cb*omega_r*phi_br - l_bb*l_cc*omega_r*phi_ar + l_bc*l_cb*omega_r*phi_ar - R_s*i_a*l_bb*l_cc + R_s*i_a*l_bc*l_cb + R_s*i_b*l_ab*l_cc - R_s*i_b*l_ac*l_cb - R_s*i_c*l_ab*l_bc + R_s*i_c*l_ac*l_bb)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
    di_b = -(l_aa*l_bc*v_c - l_aa*l_cc*v_b - l_ac*l_ba*v_c + l_ac*l_ca*v_b + l_ba*l_cc*v_a - l_bc*l_ca*v_a - l_aa*l_bc*omega_r*phi_cr + l_ac*l_ba*omega_r*phi_cr + l_aa*l_cc*omega_r*phi_br - l_ac*l_ca*omega_r*phi_br - l_ba*l_cc*omega_r*phi_ar + l_bc*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cc + R_s*i_a*l_bc*l_ca + R_s*i_b*l_aa*l_cc - R_s*i_b*l_ac*l_ca - R_s*i_c*l_aa*l_bc + R_s*i_c*l_ac*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
    di_c = (l_aa*l_bb*v_c - l_aa*l_cb*v_b - l_ab*l_ba*v_c + l_ab*l_ca*v_b + l_ba*l_cb*v_a - l_bb*l_ca*v_a - l_aa*l_bb*omega_r*phi_cr + l_ab*l_ba*omega_r*phi_cr + l_aa*l_cb*omega_r*phi_br - l_ab*l_ca*omega_r*phi_br - l_ba*l_cb*omega_r*phi_ar + l_bb*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cb + R_s*i_a*l_bb*l_ca + R_s*i_b*l_aa*l_cb - R_s*i_b*l_ab*l_ca - R_s*i_c*l_aa*l_bb + R_s*i_c*l_ab*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
    
    param_sim_h = time(i)-time(i-1);
    var_sim_iabc(:,i) = var_sim_iabc(:,i-1) + param_sim_h*[di_a;di_b;di_c];
end

figure
plot(time,var_sim_iabc)
hold on
plot(time,ia_sim,'--',time,ib_sim,'--',time,ic_sim,'--')

%% FORMA 2: suponiendo iabs conocido para todo t
var_sim_iabc = [ia_sim,ib_sim,ic_sim]';
var_sim_vabc = zeros(3,nd);

omega_r = RPM*2*pi/60;

method = 'linear';

for i=2:nd
    thetar_deg = thetar_sim_nuevo(i-1);

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
    
    % System evolution
    param_sim_h = time(i)-time(i-1);
    
    i_a_p1 = var_sim_iabc(1,i);
    i_b_p1 = var_sim_iabc(2,i);
    i_c_p1 = var_sim_iabc(3,i);

    v_a =  l_aa*(i_a_p1-i_a)/param_sim_h + l_ab*(i_b_p1-i_b)/param_sim_h + l_ac*(i_c_p1-i_c)/param_sim_h + R_s * i_a + phi_ar * omega_r;    
    v_b =  l_ba*(i_a_p1-i_a)/param_sim_h + l_bb*(i_b_p1-i_b)/param_sim_h + l_bc*(i_c_p1-i_c)/param_sim_h + R_s * i_b + phi_br * omega_r;
    v_c =  l_ca*(i_a_p1-i_a)/param_sim_h + l_cb*(i_b_p1-i_b)/param_sim_h + l_cc*(i_c_p1-i_c)/param_sim_h + R_s * i_c + phi_cr * omega_r;
    
    var_sim_vabc(:,i-1) = [v_a;v_b;v_c];
end

figure
plot(time,var_sim_vabc)
hold on
plot(time,va_sim,'--',time,vb_sim,'--',time,vc_sim,'--')
grid on