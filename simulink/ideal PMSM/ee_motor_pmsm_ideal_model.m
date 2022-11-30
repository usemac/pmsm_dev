%% PMSM Drive Data - Create Example Flux Linkage Data
close all;
clear all;
clc;

%% Motor parameters
% PM = 0.1;     % Permanent magnet flux
% N = 6;        % Number of pole pairs
% Ld = 0.0002;  % D-axis inductance
% Lq = 0.0005;  % Q-axis inductance
% L0 = 0.00018; % Zero-sequence inductance
% R_s = 0.013;   % Stator resistance
% % NOTA: con estos parámetros no estima bien usando el modelo de corrientes

% % Caso Ms=0, Lm=0
% PM = 0.114;     % Permanent magnet flux
% N = 4;        % Number of pole pairs
% Ld = 0.00363;  % D-axis inductance
% Lq = 0.00363;  % Q-axis inductance
% L0 = 0.00363; % Zero-sequence inductance
% R_s = 2.6;   % Stator resistance

% % Caso Lm=0
% PM = 0.114;     % Permanent magnet flux
% N = 4;        % Number of pole pairs
% Ld = 0.00363;  % D-axis inductance
% Lq = 0.00363;  % Q-axis inductance
% L0 = 0.00163; % Zero-sequence inductance
% R_s = 2.6;   % Stator resistance

% Caso L diferente de 0
PM = 0.114;     % Permanent magnet flux
N = 4;        % Number of pole pairs
Ld = 0.00363;  % D-axis inductance
Lq = 0.00263;  % Q-axis inductance
L0 = 0.00163; % Zero-sequence inductance
R_s = 2.6;   % Stator resistance

%% Tabulate flux linkage partial derivatives and torque in terms of id, iq and rotor angle
idVec = linspace(-250,250,5);
iqVec = idVec;
xVec = pi/180*linspace(0,360/N,180/N+1); % Rotor angle vector
[FluxMatrix,TorqueMatrix,dfluxAdiaMatrix,dfluxAdibMatrix,dfluxAdicMatrix,dfluxAdxMatrix] = ...
                        ee_generateIdealPMSMfluxData(PM,Ld,Lq,L0,idVec,iqVec,xVec);


%% Controller parameters

% Ts = 2e-6;                % Fundamental sample time
% fsw = 2e3;                % Switching frequency (Hz)
% fc = fsw*10;              % Control loop frequency (Hz)
% Tsc = 1/fc;               % Control loop sample time
Imax = 230;               % Assumed max stator current (peak value)
Tmax = 1.5*N*PM*Imax;     % Maximum electromagnetic torque
fnom = 140;               % Nominal frequency (Hz)                  % Esta sería fe (la relacionada con we)
rpm_nom = 60*fnom/N;      % Nominal rotor speed in rpm              % Esta es wr, la mecánica (wm para mí)
omegam_nom = 2*pi*fnom/N; % Nominal mechanical rotor speed (rad/s)  % Esta es wr, la mecánica (wm para mí)
Pmax = omegam_nom*Tmax;   % Maximum power
tsim = 0.5;
Ts = 1e-4;

f_test = 140;
Rload = 10;

% PRIMERO SE EJECUTA EL MODELO DE SIMULINK PARA OBTENER LOS RESULTADOS
sim('ee_motor_pmsm_ideal')

%% Prueba de simulación del sistema en modo generador modelo en derivadas
% Lectura de los resultados
time = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.angular_velocity.series.time;
thetar_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.angular_position.series.values;     % [deg] rotor angle, from wr
omegar_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.angular_velocity.series.values;     % [rad/s] rotor mechanical speed wr
torque_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.electrical_torque.series.values;
ia_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.i_a.series.values;
ib_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.i_b.series.values;
ic_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.i_c.series.values;
id_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.i_d.series.values;
iq_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.i_q.series.values;
va_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.v_a.series.values;
vb_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.v_b.series.values;
vc_sim = simlog_ee_motor_pmsm_ideal.FEM_Parameterized_PMSM_Short_Circuit_Test.v_c.series.values;

nd = length(time);

% % Cálculo de thetar a partir de la velocidad de simulink
% thetar = zeros(1,length(time));
% 
% thetar(1) = 0;
% thetar_ant = thetar(1);
% 
% for i = 2:nd
%     omegar = omegar_sim(i);
%     thetar(i) = thetar_ant + omegar*(time(i)-time(i-1));
%     thetar_ant = thetar(i);
% end
% 
% figure
% plot(time,thetar_sim)
% hold on
% plot(time,thetar*180/pi,'--')
% grid on
% % CONCLUSIÓN: El ángulo thetar de simlog es el ángulo thertar, y la velocidad medida es
% % la omegar (la mecánica)

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
             -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];
    invT = inv(T);
    iabc(:,i) = invT(:,1:2)*[id_sim(i);iq_sim(i)];
    idq(:,i) = T*[ia_sim(i);ib_sim(i);ic_sim(i)];
    torque_cal(i) = interpn(idVec,iqVec',xVec,TorqueMatrix,id_sim(i),iq_sim(i),thetar_sim_nuevo(i)*pi/180,"makima");
end

figure
plot(time,ia_sim,time,ib_sim,time,ic_sim)
hold on
plot(time,iabc,'--')
% CONCLUSIÓN: Coinciden

figure
plot(time,id_sim,time,iq_sim)
hold on
plot(time,idq,'--')
% CONCLUSIÓN: Coinciden

figure
plot(time,torque_sim)
hold on
plot(time,torque_cal,'--')
% CONCLUSIÓN: Error de 0.01 en el valor con método de interpolación 'linear', pero misma forma de la curva.
% Coinciden con el método 'makima'.


% %% Intento de simulación discreta del sistema conocidas vabc
% var_sim_iabc = zeros(3,nd);
% var_sim_iabc(:,1) = [ia_sim(1);ib_sim(1);ic_sim(1)];
% 
% omega_r = omegam_nom;   % wr es la velocidad mecánica del rotor
% 
% method = 'linear';
% 
% for i = 2:nd
% 
%     thetar_deg = thetar_sim_nuevo(i-1);
%     thetara_deg = thetar_deg;
%     thetae_rad = thetar_deg*N*pi/180;
%     thetarbb_deg = thetar_deg-360/3/N;
%     if thetarbb_deg < 0
%         thetarbb_deg = thetarbb_deg + 60;
%     elseif thetarbb_deg > 60
%         thetarbb_deg = thetarbb_deg - 60;
%     end
%     thetarcc_deg = thetar_deg+360/3/N;
%     if thetarcc_deg < 0
%         thetarcc_deg = thetarcc_deg + 60;
%     elseif thetarcc_deg > 60
%         thetarcc_deg = thetarcc_deg - 60;
%     end
% 
%     thetara_rad = thetara_deg*pi/180;
%     thetarbb_rad = thetarbb_deg*pi/180;
%     thetarcc_rad = thetarcc_deg*pi/180;
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
% %     if i_dqz(1) < idVec(1) || i_dqz(1) > idVec(length(idVec)) || i_dqz(2) < iqVec(1) || i_dqz(2) > iqVec(length(iqVec)) || thetar_deg < xVec(1) || thetar_deg > xVec(length(xVec))
% %         method = 'makima';
% %     else
% %         method = 'linear';
% %     end
% 
%     % Flux derivatives
%     l_aa = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz(1),i_dqz(2),thetara_rad,method);
%     l_ab = interpn(idVec,iqVec',xVec,dfluxAdibMatrix,i_dqz(1),i_dqz(2),thetara_rad,method);
%     l_ac = interpn(idVec,iqVec',xVec,dfluxAdicMatrix,i_dqz(1),i_dqz(2),thetara_rad,method);
% 
%     l_ba = l_ab;
%     l_bb = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz(1),i_dqz(2),thetarbb_rad,method);
%     l_bc = interpn(idVec,iqVec',xVec,dfluxAdibMatrix,i_dqz(1),i_dqz(2),thetarcc_rad,method);
%     
%     l_ca = l_ac;
%     l_cb = l_bc;
%     l_cc = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz(1),i_dqz(2),thetarcc_rad,method);
%     
%     phi_ar = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz(1),i_dqz(2),thetara_rad,method);
%     phi_br = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz(1),i_dqz(2),thetarbb_rad,method);
%     phi_cr = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz(1),i_dqz(2),thetarcc_rad,method);
%     
%     % System evolution
%     divi(i) = l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca;
%     di_a(i) = (l_ab*l_bc*v_c - l_ab*l_cc*v_b - l_ac*l_bb*v_c + l_ac*l_cb*v_b + l_bb*l_cc*v_a - l_bc*l_cb*v_a - l_ab*l_bc*omega_r*phi_cr + l_ac*l_bb*omega_r*phi_cr + l_ab*l_cc*omega_r*phi_br - l_ac*l_cb*omega_r*phi_br - l_bb*l_cc*omega_r*phi_ar + l_bc*l_cb*omega_r*phi_ar - R_s*i_a*l_bb*l_cc + R_s*i_a*l_bc*l_cb + R_s*i_b*l_ab*l_cc - R_s*i_b*l_ac*l_cb - R_s*i_c*l_ab*l_bc + R_s*i_c*l_ac*l_bb)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
%     di_b(i) = -(l_aa*l_bc*v_c - l_aa*l_cc*v_b - l_ac*l_ba*v_c + l_ac*l_ca*v_b + l_ba*l_cc*v_a - l_bc*l_ca*v_a - l_aa*l_bc*omega_r*phi_cr + l_ac*l_ba*omega_r*phi_cr + l_aa*l_cc*omega_r*phi_br - l_ac*l_ca*omega_r*phi_br - l_ba*l_cc*omega_r*phi_ar + l_bc*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cc + R_s*i_a*l_bc*l_ca + R_s*i_b*l_aa*l_cc - R_s*i_b*l_ac*l_ca - R_s*i_c*l_aa*l_bc + R_s*i_c*l_ac*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
%     di_c(i) = (l_aa*l_bb*v_c - l_aa*l_cb*v_b - l_ab*l_ba*v_c + l_ab*l_ca*v_b + l_ba*l_cb*v_a - l_bb*l_ca*v_a - l_aa*l_bb*omega_r*phi_cr + l_ab*l_ba*omega_r*phi_cr + l_aa*l_cb*omega_r*phi_br - l_ab*l_ca*omega_r*phi_br - l_ba*l_cb*omega_r*phi_ar + l_bb*l_ca*omega_r*phi_ar - R_s*i_a*l_ba*l_cb + R_s*i_a*l_bb*l_ca + R_s*i_b*l_aa*l_cb - R_s*i_b*l_ab*l_ca - R_s*i_c*l_aa*l_bb + R_s*i_c*l_ab*l_ba)/(l_aa*l_bb*l_cc - l_aa*l_bc*l_cb - l_ab*l_ba*l_cc + l_ab*l_bc*l_ca + l_ac*l_ba*l_cb - l_ac*l_bb*l_ca);
%     
%     param_sim_h = time(i)-time(i-1);
%     var_sim_iabc(:,i) = var_sim_iabc(:,i-1) + param_sim_h*[di_a(i);di_b(i);di_c(i)];
% end
% 
% figure
% plot(time,var_sim_iabc)
% hold on
% plot(time,ia_sim,'--',time,ib_sim,'--',time,ic_sim,'--')
% 
% % CONCLUSIÓN: No converje porde divi es muy próximo a cero.
% 
% %% FORMA 2: suponiendo iabs conocido para todo t
% var_sim_iabc = zeros(3,nd);
% var_sim_iabc(1,:) = ia_sim';
% var_sim_iabc(2,:) = ib_sim';
% var_sim_iabc(3,:) = ic_sim';
% var_sim_vabc = zeros(3,nd);
% i_dqz = zeros(3,nd);
% 
% Ls = (L0+Ld+Lq)/3;
% Lm = (Ld-Lq)/3;
% Ms = (Ld+Lq)/6-L0/3;
% 
% omega_r = omegam_nom;   % wr es la velocidad mecánica del rotor
% 
% method = 'makima';
% 
% Laa = zeros(1,nd);
% Lab = zeros(1,nd);
% Lac = zeros(1,nd);
% Lba = zeros(1,nd);
% Lbb = zeros(1,nd);
% Lbc = zeros(1,nd);
% Lca = zeros(1,nd);
% Lcb = zeros(1,nd);
% Lcc = zeros(1,nd);
% phi_ma = zeros(1,nd);
% phi_mb = zeros(1,nd);
% phi_mc = zeros(1,nd);
% 
% l_aa = zeros(1,nd);
% l_ab = zeros(1,nd);
% l_ac = zeros(1,nd);
% l_ba = zeros(1,nd);
% l_bb = zeros(1,nd);
% l_bc = zeros(1,nd);
% l_ca = zeros(1,nd);
% l_cb = zeros(1,nd);
% l_cc = zeros(1,nd);
% phi_ar = zeros(1,nd);
% phi_br = zeros(1,nd);
% phi_cr = zeros(1,nd);
% 
% Fa = zeros(1,nd);
% f_a = zeros(1,nd);
% der_Fa_dx = zeros(1,nd);
% der_phim_dx = zeros(1,nd);
% divi = zeros(1,nd);
% 
% for i = 1:nd-1
%     
%     thetar_deg = thetar_sim_nuevo(i);
%     thetara_deg = thetar_deg;
%     thetae_rad = thetar_deg*N*pi/180;
%     
%     thetarbb_deg = thetar_deg-360/3/N;
%     if thetarbb_deg < 0
%         thetarbb_deg = thetarbb_deg + 360/N;
%     elseif thetarbb_deg > 360/N
%         thetarbb_deg = thetarbb_deg - 360/N;
%     end
%     thetarcc_deg = thetar_deg+360/3/N;
%     if thetarcc_deg < 0
%         thetarcc_deg = thetarcc_deg + 360/N;
%     elseif thetarcc_deg > 360/N
%         thetarcc_deg = thetarcc_deg - 360/N;
%     end
% 
%     thetara_rad = thetara_deg*pi/180;
%     thetarbb_rad = thetarbb_deg*pi/180;
%     thetarcc_rad = thetarcc_deg*pi/180;
% 
%     i_a = var_sim_iabc(1,i);
%     i_b = var_sim_iabc(2,i);
%     i_c = var_sim_iabc(3,i);
% 
%     T = 2/3*[ cos(thetae_rad)    cos(thetae_rad - 2*pi/3)    cos(thetae_rad + 2*pi/3);
%              -sin(thetae_rad)   -sin(thetae_rad - 2*pi/3)   -sin(thetae_rad + 2*pi/3);
%               1/2                1/2                         1/2];
%     invT = inv(T);
% 
%     i_dqz(:,i) = T*var_sim_iabc(:,i);
%     
%     angulo_sim = thetae_rad;
% 
%     % Inductancias y flujo de magnetización
%     Laa(i) = Ls + Lm*cos(2*angulo_sim);
%     Lbb(i) = Ls + Lm*cos(2*(angulo_sim-2*pi/3));
%     Lcc(i) = Ls + Lm*cos(2*(angulo_sim+2*pi/3));
%     Lab(i) = - Ms - Lm*cos(2*(angulo_sim+pi/6));
%     Lbc(i) = - Ms - Lm*cos(2*(angulo_sim+pi/6-2*pi/3));
%     Lca(i) = - Ms - Lm*cos(2*(angulo_sim+pi/6+2*pi/3));
%     Lba(i) = Lab(i);
%     Lcb(i) = Lbc(i);
%     Lac(i) = Lca(i);
%     phi_ma(i) = PM*cos(angulo_sim);
%     phi_mb(i) = PM*cos(angulo_sim-2*pi/3);
%     phi_mc(i) = PM*cos(angulo_sim+2*pi/3);
% 
%     % Flux derivatives
%     l_aa(i) = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz(1,i),i_dqz(2,i),thetara_rad,method);
%     l_ab(i) = interpn(idVec,iqVec',xVec,dfluxAdibMatrix,i_dqz(1,i),i_dqz(2,i),thetara_rad,method);
%     l_ac(i) = interpn(idVec,iqVec',xVec,dfluxAdicMatrix,i_dqz(1,i),i_dqz(2,i),thetara_rad,method);
%     
%     l_ba(i) = l_ab(i);
%     l_bb(i) = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz(1,i),i_dqz(2,i),thetarbb_rad,method);
%     l_bc(i) = interpn(idVec,iqVec',xVec,dfluxAdibMatrix,i_dqz(1,i),i_dqz(2,i),thetarcc_rad,method);
%     
%     l_ca(i) = l_ac(i);
%     l_cb(i) = l_bc(i);
%     l_cc(i) = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz(1,i),i_dqz(2,i),thetarcc_rad,method);
%     
%     phi_ar(i) = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz(1,i),i_dqz(2,i),thetara_rad,method);
%     phi_br(i) = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz(1,i),i_dqz(2,i),thetarbb_rad,method);
%     phi_cr(i) = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz(1,i),i_dqz(2,i),thetarcc_rad,method);
%     
%     % System evolution
%     param_sim_h = time(i+1)-time(i);
%     
%     i_a_p1 = var_sim_iabc(1,i+1);
%     i_b_p1 = var_sim_iabc(2,i+1);
%     i_c_p1 = var_sim_iabc(3,i+1);
% 
%     v_a =  l_aa(i)*(i_a_p1-i_a)/param_sim_h + l_ab(i)*(i_b_p1-i_b)/param_sim_h + l_ac(i)*(i_c_p1-i_c)/param_sim_h + R_s * i_a + phi_ar(i) * omega_r;
%     v_b =  l_ba(i)*(i_a_p1-i_a)/param_sim_h + l_bb(i)*(i_b_p1-i_b)/param_sim_h + l_bc(i)*(i_c_p1-i_c)/param_sim_h + R_s * i_b + phi_br(i) * omega_r;
%     v_c =  l_ca(i)*(i_a_p1-i_a)/param_sim_h + l_cb(i)*(i_b_p1-i_b)/param_sim_h + l_cc(i)*(i_c_p1-i_c)/param_sim_h + R_s * i_c + phi_cr(i) * omega_r;
% 
%     % Flux A
%     Fa(i) = interpn(idVec,iqVec',xVec,FluxMatrix,i_dqz(1,i),i_dqz(2,i),thetara_rad,method);
%     f_a(i) = Laa(i)*i_a+Lab(i)*i_b+Lac(i)*i_c+phi_ma(i);
% 
%     % Derivatives
%     der_Fa_dx(i) = -2*N*Lm*sin(2*angulo_sim)*i_a+2*N*Lm*sin(2*(angulo_sim+pi/6))*i_b+2*N*Lm*sin(2*(angulo_sim+pi/6+2*pi/3))*i_c-PM*N*sin(angulo_sim);
%     der_phim_dx(i) = -PM*N*sin(angulo_sim);
% 
%     divi(i) = l_aa(i)*l_bb(i)*l_cc(i) - l_aa(i)*l_bc(i)*l_cb(i) - l_ab(i)*l_ba(i)*l_cc(i) + l_ab(i)*l_bc(i)*l_ca(i) + l_ac(i)*l_ba(i)*l_cb(i) - l_ac(i)*l_bb(i)*l_ca(i);
%     
%     var_sim_vabc(:,i) = [v_a;v_b;v_c];
% end
% 
% figure
% plot(time,var_sim_vabc)
% hold on
% plot(time,va_sim,'--',time,vb_sim,'--',time,vc_sim,'--')
% grid on
% 
% % figure
% % plot(time,i_dqz(1,:),time,i_dqz(2,:),time,i_dqz(3,:))
% 
% % CONCLUSIÓN: Parece que va bien.
% 
% % laa = dFA/diA
% figure
% plot(time,Laa)
% hold on
% plot(time,l_aa,'--')
% grid on
% legend('L_{aa}','dF_a/di_a')
% 
% % lab = dFA/diB
% figure
% plot(time,Lab)
% hold on
% plot(time,l_ab,'--')
% legend('L_{ab}','dF_a/di_b')
% 
% % lac = dFA/diC
% figure
% plot(time,Lac)
% hold on
% plot(time,l_ac,'--')
% legend('L_{ac}','dF_a/di_c')
% % CONCLUSIÓN: Las derivadas respecto a las corrientes coinciden con las inductancias.
% 
% % dFA/dx
% figure
% plot(time,phi_ar)
% hold on
% plot(time,der_Fa_dx,'--')
% plot(time,der_phim_dx)
% legend('dF_a/d\theta_r','dF_a/d\theta_r cal','d\phi_{ma}/d\theta_r cal')
% % CONCLUSIÓN: Las derivadas respecto de thetar coinciden.
% 
% % Fa = f_a
% figure
% plot(time,Fa)
% hold on
% plot(time,f_a,'--')
% grid on
% legend('\phi_{a}','\phi_{acal}')
% % CONCLUSIÓN: El flujo coincide

%% FORMA 3: matricial
var_sim_iabc = zeros(3,nd);
var_sim_iabc_2 = zeros(3,nd);
var_sim_iabc(:,1) = [ia_sim(1);ib_sim(1);ic_sim(1)];
var_sim_iabc_2(:,1) =  [ia_sim(1);ib_sim(1);ic_sim(1)];
var_sim_vabc = zeros(3,nd);
var_sim_vabc(1,:) = va_sim';
var_sim_vabc(2,:) = vb_sim';
var_sim_vabc(3,:) = vc_sim';
i_dqz = zeros(3,nd);

omega_r = omegam_nom;   % wr es la velocidad mecánica del rotor

% Parámetros modelo ideal
Ls = (L0+Ld+Lq)/3;
Lm = (Ld-Lq)/3;
Ms = (Ld+Lq)/6-L0/3;

method = 'makima';

for i = 1:nd-1
    
    thetar_deg = thetar_sim_nuevo(i);
    thetara_deg = thetar_deg;
    thetae_rad = thetar_deg*N*pi/180;
    
    thetarbb_deg = thetar_deg-360/3/N;
    if thetarbb_deg < 0
        thetarbb_deg = thetarbb_deg + 360/N;
    elseif thetarbb_deg > 360/N
        thetarbb_deg = thetarbb_deg - 360/N;
    end
    thetarcc_deg = thetar_deg+360/3/N;
    if thetarcc_deg < 0
        thetarcc_deg = thetarcc_deg + 360/N;
    elseif thetarcc_deg > 360/N
        thetarcc_deg = thetarcc_deg - 360/N;
    end

    thetara_rad = thetara_deg*pi/180;
    thetarbb_rad = thetarbb_deg*pi/180;
    thetarcc_rad = thetarcc_deg*pi/180;

    angulo_sim = thetae_rad;
    
%     % MODELO IDEAL (Lm=0, Ms=0)
%     % Inductancias y flujo de magnetización
%     Laa = Ls + Lm*cos(2*angulo_sim);
%     Lbb = Ls + Lm*cos(2*(angulo_sim-2*pi/3));
%     Lcc = Ls + Lm*cos(2*(angulo_sim+2*pi/3));
%     Lab = - Ms - Lm*cos(2*(angulo_sim+pi/6));
%     Lbc = - Ms - Lm*cos(2*(angulo_sim+pi/6-2*pi/3));
%     Lca = - Ms - Lm*cos(2*(angulo_sim+pi/6+2*pi/3));
%     Lba = Lab;
%     Lcb = Lbc;
%     Lac = Lca;
%     phi_ma = PM*cos(angulo_sim);
%     phi_mb = PM*cos(angulo_sim-2*pi/3);
%     phi_mc = PM*cos(angulo_sim+2*pi/3);
% 
%     Matrix_dF_di = [Laa, Lab, Lac; Lba, Lbb, Lbc; Lca, Lcb, Lcc];
%     Matrix_dF_dx = [-N*PM*sin(angulo_sim); -N*PM*sin(angulo_sim-2*pi/3); -N*PM*sin(angulo_sim+2*pi/3)];
% 
%     % MODELO IDEAL (Lm=0)
%     % Inductancias y flujo de magnetización
%     Laa = Ls + Lm*cos(2*angulo_sim);
%     Lbb = Ls + Lm*cos(2*(angulo_sim-2*pi/3));
%     Lcc = Ls + Lm*cos(2*(angulo_sim+2*pi/3));
%     Lab = - Ms - Lm*cos(2*(angulo_sim+pi/6));
%     Lbc = - Ms - Lm*cos(2*(angulo_sim+pi/6-2*pi/3));
%     Lca = - Ms - Lm*cos(2*(angulo_sim+pi/6+2*pi/3));
%     Lba = Lab;
%     Lcb = Lbc;
%     Lac = Lca;
%     phi_ma = PM*cos(angulo_sim);
%     phi_mb = PM*cos(angulo_sim-2*pi/3);
%     phi_mc = PM*cos(angulo_sim+2*pi/3);
% 
%     Matrix_dF_di = [Laa, Lab, Lac; Lba, Lbb, Lbc; Lca, Lcb, Lcc];
%     Matrix_dF_dx = [-N*PM*sin(angulo_sim); -N*PM*sin(angulo_sim-2*pi/3); -N*PM*sin(angulo_sim+2*pi/3)];

    % MODELO IDEAL (L diferente de 0)
    % Inductancias y flujo de magnetización
    Laa = Ls + Lm*cos(2*angulo_sim);
    Lbb = Ls + Lm*cos(2*(angulo_sim-2*pi/3));
    Lcc = Ls + Lm*cos(2*(angulo_sim+2*pi/3));
    Lab = - Ms - Lm*cos(2*(angulo_sim+pi/6));
    Lbc = - Ms - Lm*cos(2*(angulo_sim+pi/6-2*pi/3));
    Lca = - Ms - Lm*cos(2*(angulo_sim+pi/6+2*pi/3));
    Lba = Lab;
    Lcb = Lbc;
    Lac = Lca;
    phi_ma = PM*cos(angulo_sim);
    phi_mb = PM*cos(angulo_sim-2*pi/3);
    phi_mc = PM*cos(angulo_sim+2*pi/3);

    Matrix_dF_di = [Laa, Lab, Lac; Lba, Lbb, Lbc; Lca, Lcb, Lcc];
    Matrix_dF_dx = [-2*N*Lm*sin(2*angulo_sim)*var_sim_iabc(1,i)+2*N*Lm*sin(2*(angulo_sim+pi/6))*var_sim_iabc(2,i)+2*N*Lm*sin(2*(angulo_sim+pi/6+2*pi/3))*var_sim_iabc(3,i)-N*PM*sin(angulo_sim); 
                     2*N*Lm*sin(2*(angulo_sim+pi/6))*var_sim_iabc(1,i)-2*N*Lm*sin(2*(angulo_sim-2*pi/3))*var_sim_iabc(2,i)+2*N*Lm*sin(2*(angulo_sim+pi/6-2*pi/3))*var_sim_iabc(3,i)-N*PM*sin(angulo_sim-2*pi/3);
                     2*N*Lm*sin(2*(angulo_sim+pi/6+2*pi/3))*var_sim_iabc(1,i)+2*N*Lm*sin(2*(angulo_sim+pi/6-2*pi/3))*var_sim_iabc(2,i)-2*N*Lm*sin(2*(angulo_sim+2*pi/3))*var_sim_iabc(3,i)-N*PM*sin(angulo_sim+2*pi/3)];

    % Evolución del sistema
    Ts = time(i+1)-time(i);
    di_dt = Matrix_dF_di\(var_sim_vabc(:,i) - R_s*var_sim_iabc(:,i) - omega_r*Matrix_dF_dx);
    var_sim_iabc(:,i+1) = var_sim_iabc(:,i) + Ts*di_dt;
    
    % MODELO EN DERIVADAS
    % Park
    T = 2/3*[ cos(thetae_rad)    cos(thetae_rad - 2*pi/3)    cos(thetae_rad + 2*pi/3);
             -sin(thetae_rad)   -sin(thetae_rad - 2*pi/3)   -sin(thetae_rad + 2*pi/3);
              1/2                1/2                         1/2];
    invT = inv(T);

    i_dqz(:,i) = T*var_sim_iabc_2(:,i);

    % Flux derivatives
    l_aa = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz(1,i),i_dqz(2,i),thetara_rad,method);
    l_ab = interpn(idVec,iqVec',xVec,dfluxAdibMatrix,i_dqz(1,i),i_dqz(2,i),thetara_rad,method);
    l_ac = interpn(idVec,iqVec',xVec,dfluxAdicMatrix,i_dqz(1,i),i_dqz(2,i),thetara_rad,method);
    
    l_ba = l_ab;
    l_bb = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz(1,i),i_dqz(2,i),thetarbb_rad,method);
    l_bc = interpn(idVec,iqVec',xVec,dfluxAdibMatrix,i_dqz(1,i),i_dqz(2,i),thetarcc_rad,method);
    
    l_ca = l_ac;
    l_cb = l_bc;
    l_cc = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz(1,i),i_dqz(2,i),thetarcc_rad,method);
    
    phi_ar = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz(1,i),i_dqz(2,i),thetara_rad,method);
    phi_br = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz(1,i),i_dqz(2,i),thetarbb_rad,method);
    phi_cr = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz(1,i),i_dqz(2,i),thetarcc_rad,method);

    % Matriz de derivadas del flujo respecto a la corriente
    Matrix_dF_di = [l_aa, l_ab, l_ac; l_ba, l_bb, l_bc; l_ca, l_cb, l_cc];
    Matrix_dF_dx = [phi_ar; phi_br; phi_cr];

    % Evolución del sistema
    Ts = time(i+1)-time(i);
    di_dt = Matrix_dF_di\(var_sim_vabc(:,i) - R_s*var_sim_iabc_2(:,i) - omega_r*Matrix_dF_dx);
    var_sim_iabc_2(:,i+1) = var_sim_iabc_2(:,i) + Ts*di_dt;    

end

figure
plot(time,var_sim_iabc,'*')
hold on
plot(time,var_sim_iabc_2)
plot(time,ia_sim,'--',time,ib_sim,'--',time,ic_sim,'--')
grid on

% CONCLUSIÓN: Hay que revisar el modelo en derivadas parciales porque no
% termina de converger al valor final por algún motivo
