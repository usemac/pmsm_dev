%% PMSM Drive Data - Create Example Flux Linkage Data
clear all;
clc;

%% Motor parameters
PM = 0.1;     % Permanent magnet flux
N = 6;        % Number of pole pairs
Ld = 0.0002;  % D-axis inductance
Lq = 0.0002;  % Q-axis inductance
L0 = 0.00018; % Zero-sequence inductance
R_s = 0.013;   % Stator resistance

%% Tabulate flux linkage partial derivatives and torque in terms of id, iq and rotor angle
idVec = linspace(-250,250,5);
iqVec = idVec;
xVec = pi/180*linspace(0,360/N,180/N+1); % Rotor angle vector
[~,TorqueMatrix,dfluxAdiaMatrix,dfluxAdibMatrix,dfluxAdicMatrix,dfluxAdxMatrix] = ...
                        ee_generateIdealPMSMfluxData(PM,Ld,Lq,L0,idVec,iqVec,xVec);

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
time = ang_speed.time;
thetar_sim = ang_position.signals.values;   % [deg] 
omegar_sim = ang_speed.signals.values;      % [rad/s]
torque_sim = torque.signals.values;
vabc_sim = phase_voltages.signals.values;
iabc_sim = phase_currents.signals.values;
id_sim = dq_currents.signals(1).values;
iq_sim = dq_currents.signals(2).values;

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
    iabc(:,i) = invT(:,1:2)*[-id_sim(i);-iq_sim(i)];
    idq(:,i) = T*(-iabc_sim(i,:))';
    torque_cal(i) = interp3(iqVec,idVec',xVec,TorqueMatrix,-iq_sim(i),-id_sim(i),thetar_sim_nuevo(i)*pi/180,"linear");
end

figure
plot(time,iabc_sim)
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


%% Intento de simulación del sistema conocidas vabc

var_sim_iabc = zeros(3,nd);
var_sim_iabc(:,1) = iabc_sim(1,:)';

omega_r = omegam_nom;   % wr es la velocidad mecánica del rotor

method = 'linear';

for i = 2:nd
    thetar_deg = thetar_sim_nuevo(i-1);
    thetae_rad = thetar_deg*N*pi/180;
    
    thetara_deg = thetar_deg;
    thetarb_deg = thetar_deg-360/3/N;
    if thetarb_deg < 0
        thetarb_deg = thetarb_deg + 60;
    elseif thetarb_deg > 60
        thetarb_deg = thetarb_deg - 60;
    end
    thetarc_deg = thetar_deg-2*360/3/N;
    if thetarc_deg < 0
        thetarc_deg = thetarc_deg + 60;
    elseif thetarc_deg > 60
        thetarc_deg = thetarc_deg - 60;
    end
    
    thetara_rad = thetara_deg*pi/180;
    thetarb_rad = thetarb_deg*pi/180;
    thetarc_rad = thetarc_deg*pi/180;

    v_a = vabc_sim(i-1,1);
    v_b = vabc_sim(i-1,2);
    v_c = vabc_sim(i-1,3);

    i_a = var_sim_iabc(1,i-1);
    i_b = var_sim_iabc(2,i-1);
    i_c = var_sim_iabc(3,i-1);

    T = 2/3*[ cos(thetae_rad)    cos(thetae_rad - 2*pi/3)    cos(thetae_rad + 2*pi/3);
             -sin(thetae_rad)   -sin(thetae_rad - 2*pi/3)   -sin(thetae_rad + 2*pi/3);
              1/2                1/2                         1/2];
    invT = inv(T);

    i_dqz = T*var_sim_iabc(:,i-1);

    % Flux derivatives
    l_aa = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    l_ab = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    l_ac = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    phi_ar = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    
    l_ba = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
    l_bb = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
    l_bc = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
    phi_br = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
    
    l_ca = interpn(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
    l_cb = interpn(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
    l_cc = interpn(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
    phi_cr = interpn(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
    
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
plot(time,iabc_sim,'--')

%% FORMA 2: suponiendo iabs conocido para todo t
var_sim_iabc = iabc_sim';
var_sim_vabc = zeros(3,nd);
i_dqz = zeros(3,nd);

omega_r = omegam_nom;   % wr es la velocidad mecánica del rotor

method = 'linear';

for i = 2:nd

    thetar_deg = thetar_sim_nuevo(i-1);
    thetae_rad = thetar_deg*N*pi/180;
    
    thetara_deg = thetar_deg;
    thetarb_deg = thetar_deg-360/3/N;
    if thetarb_deg < 0
        thetarb_deg = thetarb_deg + 60;
    elseif thetarb_deg > 60
        thetarb_deg = thetarb_deg - 60;
    end
    thetarc_deg = thetar_deg-2*360/3/N;
    if thetarc_deg < 0
        thetarc_deg = thetarc_deg + 60;
    elseif thetarc_deg > 60
        thetarc_deg = thetarc_deg - 60;
    end
    
    thetara_rad = thetara_deg*pi/180;
    thetarb_rad = thetarb_deg*pi/180;
    thetarc_rad = thetarc_deg*pi/180;

    i_a = var_sim_iabc(1,i-1);
    i_b = var_sim_iabc(2,i-1);
    i_c = var_sim_iabc(3,i-1);

    T = 2/3*[ cos(thetae_rad)    cos(thetae_rad - 2*pi/3)    cos(thetae_rad + 2*pi/3);
             -sin(thetae_rad)   -sin(thetae_rad - 2*pi/3)   -sin(thetae_rad + 2*pi/3);
              1/2                1/2                         1/2];
    invT = inv(T);

    i_dqz(:,i-1) = T*var_sim_iabc(:,i-1);

    % Flux derivatives
    l_aa = interp3(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    l_ab = interp3(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    l_ac = interp3(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    phi_ar = interp3(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetara_rad,method);
    
    l_ba = interp3(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
    l_bb = interp3(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
    l_bc = interp3(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
    phi_br = interp3(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetarb_rad,method);
    
    l_ca = interp3(iqVec,idVec',xVec,dfluxAdiaMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
    l_cb = interp3(iqVec,idVec',xVec,dfluxAdibMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
    l_cc = interp3(iqVec,idVec',xVec,dfluxAdicMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
    phi_cr = interp3(iqVec,idVec',xVec,dfluxAdxMatrix,i_dqz(2),i_dqz(1),thetarc_rad,method);
    
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
plot(time,vabc_sim,'--')
grid on

figure
plot(time,i_dqz)