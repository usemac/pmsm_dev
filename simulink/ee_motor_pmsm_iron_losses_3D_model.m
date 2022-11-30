%% PMSM Drive Data - Create Example Flux Linkage Data
close all;
clear all;
clc;

%% Motor parameters
PM = 0.1;     % Permanent magnet flux
N = 6;        % Number of pole pairs
Ld = 0.0002;  % D-axis inductance
Lq = 0.0005;  % Q-axis inductance
L0 = 0.00018; % Zero-sequence inductance
R_s = 0.013;   % Stator resistance

%% Tabulate flux linkage partial derivatives and torque in terms of id, iq and rotor angle
idVec = linspace(-250,250,5);
iqVec = idVec;
xVec = pi/180*linspace(0,360/N,180/N+1); % Rotor angle vector
[~,TorqueMatrix,dfluxAdiaMatrix,dfluxAdibMatrix,dfluxAdicMatrix,dfluxAdxMatrix] = ...
                        ee_generateIdealPMSMfluxData(PM,Ld,Lq,L0,idVec,iqVec,xVec);


% %% Prueba de derivadas del flujo
% xVec = [0:1:360]*pi/180;
% Ls = (L0+Ld+Lq)/3;
% Lm = (Ld-Lq)/3;
% Ms = (Ld+Lq)/6-L0/3;
% 
% Laa = zeros(1,length(xVec));
% Lab = zeros(1,length(xVec));
% Lac = zeros(1,length(xVec));
% Lba = zeros(1,length(xVec));
% Lbb = zeros(1,length(xVec));
% Lbc = zeros(1,length(xVec));
% Lca = zeros(1,length(xVec));
% Lcb = zeros(1,length(xVec));
% Lcc = zeros(1,length(xVec));
% phi_ra = zeros(1,length(xVec));
% phi_rb = zeros(1,length(xVec));
% phi_rc = zeros(1,length(xVec));
% 
% for i =1:length(xVec)
%     Laa(i) = Ls + Lm*cos(2*xVec(i));
%     Lbb(i) = Ls + Lm*cos(2*(xVec(i)-2*pi/3));
%     Lcc(i) = Ls + Lm*cos(2*(xVec(i)+2*pi/3));
%     Lab(i) = - Ms - Lm*cos(xVec(i)+pi/6);
%     Lbc(i) = - Ms - Lm*cos(xVec(i)+pi/6-2*pi/3);
%     Lca(i) = - Ms - Lm*cos(xVec(i)+pi/6+2*pi/3);
%     Lba(i) = Lab(i);
%     Lcb(i) = Lbc(i);
%     Lac(i) = Lca(i);
%     phi_ra(i) = PM*cos(N*xVec(i));
%     phi_rb(i) = PM*cos(N*xVec(i)-2*pi/3);
%     phi_rc(i) = PM*cos(N*xVec(i)+2*pi/3);
% end
% 
% % Intento de sacar las derivadas en b y c a partir de las derivadas en a
% for i = 1:length(xVec)
%     thetar_bb = xVec(i) - 2*pi/3;
%     if thetar_bb < 0
%         thetar_bb = thetar_bb + 2*pi;
%     end
%     Lbb_p(i) = interp1(xVec,Laa,thetar_bb);
%     thetar_cc = xVec(i) + 2*pi/3;
%     if thetar_cc > 2*pi
%         thetar_cc = thetar_cc - 2*pi;
%     end
%     Lcc_p(i) = interp1(xVec,Laa,thetar_cc);
%     
%     thetar_bc = xVec(i) - 2*pi/3;
%     if thetar_bc < 0
%         thetar_bc = thetar_bc + 2*pi;
%     end
%     Lbc_p(i) = interp1(xVec,Lab,thetar_bc);
%     thetar_ca = xVec(i) + 2*pi/3;
%     if thetar_ca > 2*pi
%         thetar_ca = thetar_ca - 2*pi;
%     end
%     Lca_p(i) = interp1(xVec,Lab,thetar_ca);
% 
%     thetar_rc = xVec(i) - 2*pi/3;
%     if thetar_rc < 0
%         thetar_rc = thetar_rc + 2*pi;
%     end
%     phi_rc_p(i) = interp1(xVec,phi_rb,thetar_rc);
%     thetar_rb = xVec(i) + 2*pi/3;
%     if thetar_rb > 2*pi
%         thetar_rb = thetar_rb - 2*pi;
%     end
%     phi_rb_p(i) = interp1(xVec,phi_ra,thetar_rb);
% end
% 
% % laa = dFA/diA; lbb = dFB/diB; lcc = dFC/diC
% figure
% plot(xVec,Laa,xVec,Lbb,xVec,Lcc)
% hold on
% plot(xVec,Lbb_p,'--',xVec,Lcc_p,'--')
% grid on
% legend('l_{aa}','l_{bb}','l_{cc}','l_{bbp}','l_{ccp}')
% 
% % lab = dFA/diB; lba = dFB/diA
% figure
% plot(xVec,Lab)
% hold on
% plot(xVec,Lba,'--')
% legend('l_{ab}','l_{ba}')
% 
% % lac = dFA/diC; lca = dFC/diA
% figure
% plot(xVec,Lac)
% hold on
% plot(xVec,Lca,'--',xVec,Lca_p,'o')
% legend('l_{ac}','l_{ca}')
% 
% % lbc = dFB/diC; lcb = dFC/diB
% figure
% plot(xVec,Lbc)
% hold on
% plot(xVec,Lcb,'--',xVec,Lbc_p,'o')
% legend('l_{bc}','l_{cb}')
% 
% figure
% plot(xVec,phi_ra,xVec,phi_rb,xVec,phi_rc)
% hold on
% plot(xVec,phi_rb_p,'--',xVec,phi_rc_p,'--')
% grid on
% legend('\phi_{ra}','\phi_{rb}','\phi_{rc}','\phi_{rbp}','\phi_{rcp}')

%% Obtención de las derivadas del flujo para b y c a partir de a
id = 0;
iq = 100;

Ls = (L0+Ld+Lq)/3;
Lm = (Ld-Lq)/3;
Ms = (Ld+Lq)/6-L0/3;

Laa = zeros(1,length(xVec));
Lab = zeros(1,length(xVec));
Lac = zeros(1,length(xVec));
Lba = zeros(1,length(xVec));
Lbb = zeros(1,length(xVec));
Lbc = zeros(1,length(xVec));
Lca = zeros(1,length(xVec));
Lcb = zeros(1,length(xVec));
Lcc = zeros(1,length(xVec));
phi_ra = zeros(1,length(xVec));
phi_rb = zeros(1,length(xVec));
phi_rc = zeros(1,length(xVec));

for i = 1:length(xVec)

    thetar_deg = xVec(i)*180/pi;
    thetara_deg = thetar_deg;


    thetarbb_deg = thetar_deg-360/3/N;
    if thetarbb_deg < 0
        thetarbb_deg = thetarbb_deg + 60;
    elseif thetarbb_deg > 60
        thetarbb_deg = thetarbb_deg - 60;
    end
    thetarcc_deg = thetar_deg+360/3/N;
    if thetarcc_deg < 0
        thetarcc_deg = thetarcc_deg + 60;
    elseif thetarcc_deg > 60
        thetarcc_deg = thetarcc_deg - 60;
    end
    thetarbc_deg = thetar_deg-360/3/N;
    if thetarbc_deg < 0
        thetarbc_deg = thetarbc_deg + 60;
    elseif thetarbc_deg > 60
        thetarbc_deg = thetarbc_deg - 60;
    end
    thetarca_deg = thetar_deg+360/3/N;
    if thetarca_deg < 0
        thetarca_deg = thetarca_deg + 60;
    elseif thetarca_deg > 60
        thetarca_deg = thetarca_deg - 60;
    end
    thetarrb_deg = thetar_deg-360/3/N;
    if thetarrb_deg < 0
        thetarrb_deg = thetarrb_deg + 60;
    elseif thetarrb_deg > 60
        thetarrb_deg = thetarrb_deg - 60;
    end
    thetarrc_deg = thetar_deg+360/3/N;
    if thetarrc_deg < 0
        thetarrc_deg = thetarrc_deg + 60;
    elseif thetarrc_deg > 60
        thetarrc_deg = thetarrc_deg - 60;
    end

    thetara_rad = thetara_deg*pi/180;
    thetarbb_rad = thetarbb_deg*pi/180;
    thetarcc_rad = thetarcc_deg*pi/180;
    thetarbc_rad = thetarbc_deg*pi/180;
    thetarca_rad = thetarca_deg*pi/180;
    thetarrb_rad = thetarrb_deg*pi/180;
    thetarrc_rad = thetarrc_deg*pi/180;

    % Flux derivatives
    laa(i) = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,iq,id,thetara_rad);
    lab(i) = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,iq,id,thetara_rad);
    lac(i) = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,iq,id,thetara_rad);

    lbb(i) = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,iq,id,thetarbb_rad);
    lcc(i) = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,iq,id,thetarcc_rad);
    
    lbc(i) = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,iq,id,thetarbc_rad);
%     lca(i) = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,iq,id,thetarca_rad);
    lca(i) = lac(i);

    lba(i) = lab(i);
    lcb(i) = lbc(i);
    
    phi_ra(i) = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,iq,id,thetara_rad);
    phi_rb(i) = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,iq,id,thetarrb_rad);
    phi_rc(i) = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,iq,id,thetarrc_rad);
end

% laa = dFA/diA; lbb = dFB/diB; lcc = dFC/diC
figure
plot(xVec,laa,xVec,lbb,xVec,lcc)
grid on
legend('l_{aa}','l_{bb}','l_{cc}')

% lab = dFA/diB; lba = dFB/diA
figure
plot(xVec,lab)
hold on
plot(xVec,lba,'--')
legend('l_{ab}','l_{ba}')

% lac = dFA/diC; lca = dFC/diA
figure
plot(xVec,lac)
hold on
plot(xVec,lca,'--')
legend('l_{ac}','l_{ca}')

% lbc = dFB/diC; lcb = dFC/diB
figure
plot(xVec,lbc)
hold on
plot(xVec,lcb,'--')
legend('l_{bc}','l_{cb}')

figure
plot(xVec,phi_ra,xVec,phi_rb,xVec,phi_rc)
grid on
legend('\phi_{ra}','\phi_{rb}','\phi_{rc}')

%% Controller parameters

Ts = 2e-6;                % Fundamental sample time
fsw = 2e3;                % Switching frequency (Hz)
fc = fsw*10;              % Control loop frequency (Hz)
Tsc = 1/fc;               % Control loop sample time
Imax = 230;               % Assumed max stator current (peak value)
Tmax = 1.5*N*PM*Imax;     % Maximum electromagnetic torque
fnom = 140;               % Nominal frequency (Hz)                  % Esta sería fe (la relacionada con we)
rpm_nom = 60*fnom/N;      % Nominal rotor speed in rpm              % Esta es wr, la mecánica (wm para mí)
omegam_nom = 2*pi*fnom/N; % Nominal mechanical rotor speed (rad/s)  % Esta es wr, la mecánica (wm para mí)
Pmax = omegam_nom*Tmax;   % Maximum power


f_test = 140;
Rload = 10;

% PRIMERO SE EJECUTA EL MODELO DE SIMULINK PARA OBTENER LOS RESULTADOS
sim('ee_motor_pmsm_iron_losses_3D')

%% Prueba de simulación del sistema en modo generador modelo en derivadas
% Lectura de los resultados
time = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.angular_velocity.series.time;
thetar_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.angular_position.series.values;     % [deg] rotor angle, from wr
omegar_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.angular_velocity.series.values;     % [rad/s] rotor mechanical speed wr
torque_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.electrical_torque.series.values;
ia_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.i_a.series.values;
ib_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.i_b.series.values;
ic_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.i_c.series.values;
id_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.i_d.series.values;
iq_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.i_q.series.values;
va_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.v_a.series.values;
vb_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.v_b.series.values;
vc_sim = simlog_ee_motor_pmsm_iron_losses.FEM_Parameterized_PMSM_Short_Circuit_Test.v_c.series.values;

% thetar_sim = thetam_sim/N;
nd = length(time);

% Cálculo de thetar a partir de la velocidad de simulink
thetar = zeros(1,length(time));

thetar(1) = 0;
thetar_ant = thetar(1);

for i = 2:length(time)
    omegar = omegar_sim(i);
    thetar(i) = thetar_ant + omegar*(time(i)-time(i-1));
    thetar_ant = thetar(i);
end

figure
plot(time,thetar_sim)
hold on
plot(time,thetar*180/pi,'--')
grid on

% Cambio de formato de thetar
thetar_aux = thetar_sim;
thetar_sim_nuevo = zeros(1,nd);
for i = 1:nd
    if thetar_aux(i) > 60
        thetar_aux = thetar_aux - 60;
        thetar_sim_nuevo(i) = thetar_aux(i);
    else
        thetar_sim_nuevo(i) = thetar_aux(i);
    end
end

% Comprobación de corrientes abc y de par
for i = 1:nd
    thetae = N*thetar_sim(i)*pi/180;
    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
             -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];
    invT = inv(T);
    iabc(:,i) = invT(:,1:2)*[id_sim(i);iq_sim(i)];
    idq(:,i) = T*[ia_sim(i);ib_sim(i);ic_sim(i)];
    torque_cal(i) = interp3(iqVec,idVec',xVec,TorqueMatrix,iq_sim(i),id_sim(i),thetar_sim_nuevo(i)*pi/180,"linear");
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

% %% Representation of the partial derivatives
% % Operational point
% id = 0;
% iq = 100;
% 
% % T = 2/3*[ cos(theta)     cos(theta - 2*pi/3)     cos(theta + 2*pi/3);
% %          -sin(theta)    -sin(theta - 2*pi/3)    -sin(theta + 2*pi/3)
% %           1/2            1/2                     1/2];
% % 
% % invT = inv(T);
% % iabc = invT*[id;iq;0];
% 
% laa = zeros(1,length(xVec));
% lab = zeros(1,length(xVec));
% lac = zeros(1,length(xVec));
% lba = zeros(1,length(xVec));
% lbb = zeros(1,length(xVec));
% lbc = zeros(1,length(xVec));
% lca = zeros(1,length(xVec));
% lcb = zeros(1,length(xVec));
% lcc = zeros(1,length(xVec));
% 
% for i = 1:length(xVec)
%     thetar_deg = xVec(i)*180/pi;
%     
%     thetara_deg = thetar_deg;
%     thetarb_deg = thetar_deg-360/3/N;
%     if thetarb_deg < 0
%         thetarb_deg = thetarb_deg + 60;
%     elseif thetarb_deg > 60
%         thetarb_deg = thetarb_deg - 60;
%     end
%     thetarc_deg = thetar_deg-2*360/3/N;
%     if thetarc_deg < 0
%         thetarc_deg = thetarc_deg + 60;
%     elseif thetarc_deg > 60
%         thetarc_deg = thetarc_deg - 60;
%     end
%     
%     thetara_rad = thetara_deg*pi/180;
%     thetarb_rad = thetarb_deg*pi/180;
%     thetarc_rad = thetarc_deg*pi/180;
% 
%     % Flux derivatives
%     laa(i) = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,iq,id,thetara_rad);
%     lab(i) = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,iq,id,thetara_rad);
%     lac(i) = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,iq,id,thetara_rad);
%     
%     lba(i) = lab(i);
%     lbb(i) = laa(i);
%     lbc(i) = lab(i);
%     
%     lca(i) = lac(i);
%     lcb(i) = lac(i);
%     lcc(i) = laa(i);
% end
% 
% % laa = dFA/diA; lbb = dFB/diB; lcc = dFC/diC
% figure
% plot(xVec,laa,xVec,lbb,xVec,lcc)
% grid on
% legend('l_{aa}','l_{bb}','l_{cc}')
% 
% % lab = dFA/diB; lba = dFB/diA
% figure
% plot(xVec,lab)
% hold on
% plot(xVec,lba,'--')
% legend('l_{ab}','l_{ba}')
% 
% % lac = dFA/diC; lca = dFC/diA
% figure
% plot(xVec,lac)
% hold on
% plot(xVec,lca,'--')
% legend('l_{ac}','l_{ca}')
% 
% % lbc = dFB/diC; lcb = dFC/diB
% figure
% plot(xVec,lbc)
% hold on
% plot(xVec,lcb,'--')
% legend('l_{bc}','l_{cb}')
% 
% % lab = dFA/diB; lbc = dFB/diC; lca = dFC/diA
% figure
% plot(xVec,lab,xVec,lbc,xVec,lca)
% grid on
% legend('l_{ab}','l_{bc}','l_{ca}')

%% Intento de simulación del sistema conocidas vabc
var_sim_iabc = zeros(3,nd);
var_sim_iabc(:,1) = [ia_sim(1);ib_sim(1);ic_sim(1)];

omega_r = omegam_nom;   % wr es la velocidad mecánica del rotor

method = 'linear';

for i = 2:nd
%     thetar_deg = thetar_sim_nuevo(i-1);
%     thetae_rad = thetar_deg*N*pi/180;
%     
%     thetara_deg = thetar_deg;
%     thetarb_deg = thetar_deg-360/3/N;
%     if thetarb_deg < 0
%         thetarb_deg = thetarb_deg + 60;
%     elseif thetarb_deg > 60
%         thetarb_deg = thetarb_deg - 60;
%     end
%     thetarc_deg = thetar_deg-2*360/3/N;
%     if thetarc_deg < 0
%         thetarc_deg = thetarc_deg + 60;
%     elseif thetarc_deg > 60
%         thetarc_deg = thetarc_deg - 60;
%     end
%     
%     thetara_rad = thetara_deg*pi/180;
%     thetarb_rad = thetarb_deg*pi/180;
%     thetarc_rad = thetarc_deg*pi/180;
% 
%     v_a = va_sim(i-1);
%     v_b = vb_sim(i-1);
%     v_c = vc_sim(i-1);
% 
%     i_a = var_sim_iabc(1,i-1);
%     i_b = var_sim_iabc(2,i-1);
%     i_c = var_sim_iabc(3,i-1);
% 
%     T = 2/3*[ cos(thetae_rad)    cos(thetae_rad - 2*pi/3)    cos(thetae_rad + 2*pi/3);
%              -sin(thetae_rad)   -sin(thetae_rad - 2*pi/3)   -sin(thetae_rad + 2*pi/3);
%               1/2                1/2                         1/2];
%     invT = inv(T);
% 
%     i_dqz = T*var_sim_iabc(:,i-1);
% 
%     if i_dqz(1) < idVec(1) || i_dqz(1) > idVec(length(idVec)) || i_dqz(2) < iqVec(1) || i_dqz(2) > iqVec(length(iqVec)) || thetar_deg < xVec(1) || thetar_deg > xVec(length(xVec))
%         method = 'makima';
%     else
%         method = 'linear';
%     end
% 
%     % Flux derivatives
%     l_aa = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
%     l_ab = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
%     l_ac = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
%     phi_ar = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
%     
% %     l_ba = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
% %     l_bb = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
% %     l_bc = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
%     l_ba = l_ab;
%     l_bb = l_aa;
%     l_bc = l_ba;
%     phi_br = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
%     
% %     l_ca = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
% %     l_cb = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
% %     l_cc = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
%     l_ca = l_ac;
%     l_cb = l_ac;
%     l_cc = l_aa;
%     phi_cr = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
    
    thetar_deg = thetar_sim_nuevo(i-1);
    thetara_deg = thetar_deg;
    thetae_rad = thetar_deg*N*pi/180;
    thetarbb_deg = thetar_deg-360/3/N;
    if thetarbb_deg < 0
        thetarbb_deg = thetarbb_deg + 60;
    elseif thetarbb_deg > 60
        thetarbb_deg = thetarbb_deg - 60;
    end
    thetarcc_deg = thetar_deg+360/3/N;
    if thetarcc_deg < 0
        thetarcc_deg = thetarcc_deg + 60;
    elseif thetarcc_deg > 60
        thetarcc_deg = thetarcc_deg - 60;
    end
    thetarbc_deg = thetar_deg-360/3/N;
    if thetarbc_deg < 0
        thetarbc_deg = thetarbc_deg + 60;
    elseif thetarbc_deg > 60
        thetarbc_deg = thetarbc_deg - 60;
    end
    thetarca_deg = thetar_deg+360/3/N;
    if thetarca_deg < 0
        thetarca_deg = thetarca_deg + 60;
    elseif thetarca_deg > 60
        thetarca_deg = thetarca_deg - 60;
    end
    thetarrb_deg = thetar_deg-360/3/N;
    if thetarrb_deg < 0
        thetarrb_deg = thetarrb_deg + 60;
    elseif thetarrb_deg > 60
        thetarrb_deg = thetarrb_deg - 60;
    end
    thetarrc_deg = thetar_deg+360/3/N;
    if thetarrc_deg < 0
        thetarrc_deg = thetarrc_deg + 60;
    elseif thetarrc_deg > 60
        thetarrc_deg = thetarrc_deg - 60;
    end

    thetara_rad = thetara_deg*pi/180;
    thetarbb_rad = thetarbb_deg*pi/180;
    thetarcc_rad = thetarcc_deg*pi/180;
    thetarbc_rad = thetarbc_deg*pi/180;
    thetarca_rad = thetarca_deg*pi/180;
    thetarrb_rad = thetarrb_deg*pi/180;
    thetarrc_rad = thetarrc_deg*pi/180;

    v_a = va_sim(i-1);
    v_b = vb_sim(i-1);
    v_c = vc_sim(i-1);

    i_a = var_sim_iabc(1,i-1);
    i_b = var_sim_iabc(2,i-1);
    i_c = var_sim_iabc(3,i-1);

    T = 2/3*[ cos(thetae_rad)    cos(thetae_rad - 2*pi/3)    cos(thetae_rad + 2*pi/3);
             -sin(thetae_rad)   -sin(thetae_rad - 2*pi/3)   -sin(thetae_rad + 2*pi/3);
              1/2                1/2                         1/2];
    invT = inv(T);

    i_dqz = T*var_sim_iabc(:,i-1);

    if i_dqz(1) < idVec(1) || i_dqz(1) > idVec(length(idVec)) || i_dqz(2) < iqVec(1) || i_dqz(2) > iqVec(length(iqVec)) || thetar_deg < xVec(1) || thetar_deg > xVec(length(xVec))
        method = 'makima';
    else
        method = 'linear';
    end

    % Flux derivatives
    l_aa = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    l_ab = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    l_ac = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);

    l_bb = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetarbb_rad,method);
    l_cc = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetarcc_rad,method);
    
    l_bc = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetarbc_rad,method);
%     lca(i) = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,iq,id,thetarca_rad);
    l_ca = l_ac;

    l_ba = l_ab;
    l_cb = l_bc;
    
    phi_ar = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    phi_br = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetarrb_rad,method);
    phi_cr = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetarrc_rad,method);
    
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
var_sim_iabc = zeros(3,nd);
var_sim_iabc(1,:) = ia_sim';
var_sim_iabc(2,:) = ib_sim';
var_sim_iabc(3,:) = ic_sim';
var_sim_vabc = zeros(3,nd);
i_dqz = zeros(3,nd);

omega_r = omegam_nom;   % wr es la velocidad mecánica del rotor

method = 'linear';

for i = 2:nd

%     thetar_deg = thetar_sim_nuevo(i-1);
%     thetae_rad = thetar_deg*N*pi/180;
%     
%     thetara_deg = thetar_deg;
%     thetarb_deg = thetar_deg-360/3/N;
%     if thetarb_deg < 0
%         thetarb_deg = thetarb_deg + 60;
%     elseif thetarb_deg > 60
%         thetarb_deg = thetarb_deg - 60;
%     end
%     thetarc_deg = thetar_deg-2*360/3/N;
%     if thetarc_deg < 0
%         thetarc_deg = thetarc_deg + 60;
%     elseif thetarc_deg > 60
%         thetarc_deg = thetarc_deg - 60;
%     end
%     
%     thetara_rad = thetara_deg*pi/180;
%     thetarb_rad = thetarb_deg*pi/180;
%     thetarc_rad = thetarc_deg*pi/180;
    
    thetar_deg = thetar_sim_nuevo(i-1);
    thetara_deg = thetar_deg;
    thetae_rad = thetar_deg*N*pi/180;
    thetarbb_deg = thetar_deg-360/3/N;
    if thetarbb_deg < 0
        thetarbb_deg = thetarbb_deg + 60;
    elseif thetarbb_deg > 60
        thetarbb_deg = thetarbb_deg - 60;
    end
    thetarcc_deg = thetar_deg+360/3/N;
    if thetarcc_deg < 0
        thetarcc_deg = thetarcc_deg + 60;
    elseif thetarcc_deg > 60
        thetarcc_deg = thetarcc_deg - 60;
    end
    thetarbc_deg = thetar_deg-360/3/N;
    if thetarbc_deg < 0
        thetarbc_deg = thetarbc_deg + 60;
    elseif thetarbc_deg > 60
        thetarbc_deg = thetarbc_deg - 60;
    end
    thetarca_deg = thetar_deg+360/3/N;
    if thetarca_deg < 0
        thetarca_deg = thetarca_deg + 60;
    elseif thetarca_deg > 60
        thetarca_deg = thetarca_deg - 60;
    end
    thetarrb_deg = thetar_deg-360/3/N;
    if thetarrb_deg < 0
        thetarrb_deg = thetarrb_deg + 60;
    elseif thetarrb_deg > 60
        thetarrb_deg = thetarrb_deg - 60;
    end
    thetarrc_deg = thetar_deg+360/3/N;
    if thetarrc_deg < 0
        thetarrc_deg = thetarrc_deg + 60;
    elseif thetarrc_deg > 60
        thetarrc_deg = thetarrc_deg - 60;
    end

    thetara_rad = thetara_deg*pi/180;
    thetarbb_rad = thetarbb_deg*pi/180;
    thetarcc_rad = thetarcc_deg*pi/180;
    thetarbc_rad = thetarbc_deg*pi/180;
    thetarca_rad = thetarca_deg*pi/180;
    thetarrb_rad = thetarrb_deg*pi/180;
    thetarrc_rad = thetarrc_deg*pi/180;

    i_a = var_sim_iabc(1,i-1);
    i_b = var_sim_iabc(2,i-1);
    i_c = var_sim_iabc(3,i-1);

    T = 2/3*[ cos(thetae_rad)    cos(thetae_rad - 2*pi/3)    cos(thetae_rad + 2*pi/3);
             -sin(thetae_rad)   -sin(thetae_rad - 2*pi/3)   -sin(thetae_rad + 2*pi/3);
              1/2                1/2                         1/2];
    invT = inv(T);

    i_dqz(:,i-1) = T*var_sim_iabc(:,i-1);

    % Flux derivatives
    l_aa = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetara_rad,method);
    l_ab = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetara_rad,method);
    l_ac = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetara_rad,method);

    l_bb = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarbb_rad,method);
    l_cc = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarcc_rad,method);
    
    l_bc = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarbc_rad,method);
%     lca(i) = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,iq,id,thetarca_rad);
    l_ca = l_ac;

    l_ba = l_ab;
    l_cb = l_bc;
    
    phi_ar = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetara_rad,method);
    phi_br = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarrb_rad,method);
    phi_cr = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarrc_rad,method);

%     % Flux derivatives
%     l_aa = interp3(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetara_rad,method);
%     l_ab = interp3(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetara_rad,method);
%     l_ac = interp3(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetara_rad,method);
%     phi_ar = interp3(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetara_rad,method);
%     
% %     l_ba = interp3(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarb_rad,method);
% %     l_bb = interp3(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarb_rad,method);
% %     l_bc = interp3(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarb_rad,method);
%     l_ba = l_ab;
%     l_bb = l_aa;
%     l_bc = l_ba;
%     phi_br = interp3(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarb_rad,method);
%     
% %     l_ca = interp3(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarc_rad,method);
% %     l_cb = interp3(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarc_rad,method);
% %     l_cc = interp3(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarc_rad,method);
%     l_ca = l_ac;
%     l_cb = l_ac;
%     l_cc = l_aa;
%     phi_cr = interp3(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2,i-1),i_dqz(1,i-1),thetarc_rad,method);
    
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

figure
plot(time,i_dqz(1,:),time,i_dqz(2,:),time,i_dqz(3,:))