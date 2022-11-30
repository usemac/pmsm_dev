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
xVec = pi/180*linspace(0,360/N,180/N+1); % Rotor angle vector [rad]
[Fa,TorqueMatrix,dfluxAdiaMatrix,dfluxAdibMatrix,dfluxAdicMatrix,dfluxAdxMatrix] = ...
                        ee_generateIdealPMSMfluxData(PM,Ld,Lq,L0,idVec,iqVec,xVec);

%% Create new LUT for Fd and Fq

% Calculate Fb and Fc
nid = length(idVec);
niq = length(iqVec);
nang = length(xVec);

Fb = zeros(nid,niq,nang);
Fc = zeros(nid,niq,nang);
Fd = zeros(nid,niq,nang);
Fq = zeros(nid,niq,nang);

method = 'makima';

for id = 1:nid
    for iq = 1:niq
        for ang = 1:nang
            % Angles
            anga = xVec(ang);
            angb = xVec(ang)-2*pi/3/N;
            if angb < 0
                angb = angb + 2*pi/N;
            elseif angb > 2*pi/N
                angb = angb - 2*pi/N;
            end
            angc= xVec(ang)+2*pi/3/N;
            if angc < 0
                angc = angc + 2*pi/N;
            elseif angc > 2*pi/N
                angc = angc - 2*pi/N;
            end

            % Fluxes Fb and Fc
            fluxa(ang) = interpn(idVec,iqVec',xVec,Fa,idVec(id),iqVec(iq),anga,method);
            Fb(id,iq,ang) = interpn(idVec,iqVec',xVec,Fa,idVec(id),iqVec(iq),angb,method);
            fluxb(ang) = Fb(id,iq,ang);
            Fc(id,iq,ang) = interpn(idVec,iqVec',xVec,Fa,idVec(id),iqVec(iq),angc,method);
            fluxc(ang) = Fc(id,iq,ang);
            
            % Park's transformation
            thetae = N*xVec(ang);
            T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
                     -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
                      1/2            1/2                     1/2];

            % Fluxes Fd and Fq 
            fluxes_dq = T*[fluxa(ang);fluxb(ang);fluxc(ang)];
            Fd(id,iq,ang) = fluxes_dq(1);
            Fq(id,iq,ang) = fluxes_dq(2);
        end
%         figure
%         hold on
%         plot(xVec,fluxa,xVec,fluxb,xVec,fluxc)
%         legend('\phi_a','\phi_b','\phi_c')
%         grid on
    end
end

% Comprobación de los flujos
id = idVec(3);
iq = iqVec(4);

Ls = (L0+Ld+Lq)/3;
Lm = (Ld-Lq)/3;
Ms = (Ld+Lq)/6-L0/3;

fluxA = zeros(1,nang);
fluxB = zeros(1,nang);
fluxC = zeros(1,nang);
fluxD = zeros(1,nang);
fluxQ = zeros(1,nang);
fluxA_ideal = zeros(1,nang);
fluxB_ideal = zeros(1,nang);
fluxC_ideal = zeros(1,nang);
fluxD_ideal = zeros(1,nang);
fluxQ_ideal = zeros(1,nang);

for x=1:nang
    % Park's transformation
    thetae = N*xVec(x);
    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
             -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];
    iabc(:,x) = inv(T)*[id;iq;0];
    
    % Flujos por LUT
    fluxA(x) = interpn(idVec,iqVec',xVec,Fa,id,iq,xVec(x),method);
    fluxB(x) = interpn(idVec,iqVec',xVec,Fb,id,iq,xVec(x),method);
    fluxC(x) = interpn(idVec,iqVec',xVec,Fc,id,iq,xVec(x),method);
    
    fluxD(x) = interpn(idVec,iqVec',xVec,Fd,id,iq,xVec(x),method);
    fluxQ(x) = interpn(idVec,iqVec',xVec,Fq,id,iq,xVec(x),method);
    
    % Flujos por modelo ideal
    Laa = Ls + Lm*cos(2*thetae);
    Lbb = Ls + Lm*cos(2*(thetae-2*pi/3));
    Lcc = Ls + Lm*cos(2*(thetae+2*pi/3));
    Lab = - Ms - Lm*cos(2*(thetae+pi/6));
    Lbc = - Ms - Lm*cos(2*(thetae+pi/6-2*pi/3));
    Lca = - Ms - Lm*cos(2*(thetae+pi/6+2*pi/3));
    Lba = Lab;
    Lcb = Lbc;
    Lac = Lca;
    phi_am = PM*cos(thetae);
    phi_bm = PM*cos(thetae-2*pi/3);
    phi_cm = PM*cos(thetae+2*pi/3);

    fluxA_ideal(x) = Laa*iabc(1,x)+Lab*iabc(2,x)+Lac*iabc(3,x)+phi_am;
    fluxB_ideal(x) = Lba*iabc(1,x)+Lbb*iabc(2,x)+Lbc*iabc(3,x)+phi_bm;
    fluxC_ideal(x) = Lca*iabc(1,x)+Lcb*iabc(2,x)+Lcc*iabc(3,x)+phi_cm;

    fluxdq = T*[fluxA_ideal(x);fluxB_ideal(x);fluxC_ideal(x)];
    fluxD_ideal(x) = fluxdq(1);
    fluxQ_ideal(x) = fluxdq(2);
end

figure
hold on
plot(xVec,fluxA)
plot(xVec,fluxB)
plot(xVec,fluxC)
plot(xVec,fluxA_ideal,'--')
plot(xVec,fluxB_ideal,'--')
plot(xVec,fluxC_ideal,'--')
grid on
legend('\phi_a','\phi_b','\phi_c','\phi_a^{ideal}','\phi_b^{ideal}','\phi_c^{ideal}')

figure
hold on
plot(xVec,fluxD)
plot(xVec,fluxQ)
plot(xVec,fluxD_ideal,'--')
plot(xVec,fluxQ_ideal,'--')
grid on
legend('\phi_d','\phi_q','\phi_d^{ideal}','\phi_q^{ideal}')


% %% Comprobar si las LUT son o no biyectivas
% 
% % Redondeamos a 4 decimales las LUTs para facilitar la comprobación
% Fd_round = round(Fd,6);
% Fq_round = round(Fq,6);
% 
% % Comprobamos punto por punto que no haya coincidencias (sin considerar theta=0)
% for index_d = 1:nid
%     for index_q = 1:niq
%         for index_x = 2:nang-1
%             i_d = idVec(index_d);
%             i_q = iqVec(index_q);
%             theta_r = xVec(index_x);
%             phi_d = Fd_round(index_d,index_q,index_x);
%             phi_q = Fq_round(index_d,index_q,index_x);
% 
%             [row_d,col_d] = find(Fd_round(:,:,2:nang-1)==phi_d);
%             [row_q,col_q] = find(Fq_round(:,:,2:nang-1)==phi_q);
%             
%             if length(row_d)>1 || length(row_q)>1
%                 break;
%             end
%         end
%         if length(row_d)>1 || length(row_q)>1
%             break;
%         end
%     end
%     if length(row_d)>1 || length(row_q)>1
%         break;
%     end
% end
% 
% % CONCLUSIÓN: Fd y Fq, en el caso ideal, no varían con theta ni con una de
% % las corrientes. Por lo tanto, se pueden sacar las LUT de id e iq sin
% % problemas.

%% LUT para id e iq en función del flujo (no dependen del ángulo)
Fd_round = round(Fd,6);
Fq_round = round(Fq,6);

% Como Fd=f(id), id=g(Fd)
% Como Fq=f(iq), iq=g(Fq)
fdVec = sort(Fd_round(:,1,1))';
fqVec = sort(Fq_round(1,:,1));
Id = idVec;
Iq = iqVec;

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
tsim = 0.7;
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

%% Evolución del sistema: x = i, (matricial)
var_sim_iabc_ideal = zeros(3,nd);
var_sim_iabc_i = zeros(3,nd);
var_sim_iabc_f = zeros(3,nd);
var_sim_iabc_ideal(:,1) = [ia_sim(1);ib_sim(1);ic_sim(1)];
var_sim_iabc_i(:,1) =  [ia_sim(1);ib_sim(1);ic_sim(1)];
var_sim_iabc_f(:,1) =  [ia_sim(1);ib_sim(1);ic_sim(1)];
var_sim_vabc = zeros(3,nd);
var_sim_vabc(1,:) = va_sim';
var_sim_vabc(2,:) = vb_sim';
var_sim_vabc(3,:) = vc_sim';
i_dqz_i = zeros(3,nd);
i_dqz_f = zeros(3,nd);
var_sim_fd = zeros(1,nd);
var_sim_fq = zeros(1,nd);

omega_r = omegam_nom;   % wr es la velocidad mecánica del rotor

% Parámetros modelo ideal
Ls = (L0+Ld+Lq)/3;
Lm = (Ld-Lq)/3;
Ms = (Ld+Lq)/6-L0/3;

method = 'makima';

% Inicialización método x=flux
% Park
thetar_deg = thetar_sim_nuevo(1);
thetae_rad = thetar_deg*N*pi/180;
thetara_rad = thetar_deg*pi/180;

T = 2/3*[ cos(thetae_rad)    cos(thetae_rad - 2*pi/3)    cos(thetae_rad + 2*pi/3);
         -sin(thetae_rad)   -sin(thetae_rad - 2*pi/3)   -sin(thetae_rad + 2*pi/3);
          1/2                1/2                         1/2];
invT = inv(T);
i_dqz_f(:,1) = T*var_sim_iabc_f(:,1);

var_sim_fd(1) = interpn(idVec,iqVec',xVec,Fd,i_dqz_f(1,1),i_dqz_f(2,1),thetara_rad,method);
var_sim_fq(1) = interpn(idVec,iqVec',xVec,Fq,i_dqz_f(1,1),i_dqz_f(2,1),thetara_rad,method);

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

    Ts = time(i+1)-time(i);

    % Park
    T = 2/3*[ cos(thetae_rad)    cos(thetae_rad - 2*pi/3)    cos(thetae_rad + 2*pi/3);
             -sin(thetae_rad)   -sin(thetae_rad - 2*pi/3)   -sin(thetae_rad + 2*pi/3);
              1/2                1/2                         1/2];
    invT = inv(T);
    
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
    Matrix_dF_dx = [-2*N*Lm*sin(2*angulo_sim)*var_sim_iabc_ideal(1,i)+2*N*Lm*sin(2*(angulo_sim+pi/6))*var_sim_iabc_ideal(2,i)+2*N*Lm*sin(2*(angulo_sim+pi/6+2*pi/3))*var_sim_iabc_ideal(3,i)-N*PM*sin(angulo_sim); 
                     2*N*Lm*sin(2*(angulo_sim+pi/6))*var_sim_iabc_ideal(1,i)-2*N*Lm*sin(2*(angulo_sim-2*pi/3))*var_sim_iabc_ideal(2,i)+2*N*Lm*sin(2*(angulo_sim+pi/6-2*pi/3))*var_sim_iabc_ideal(3,i)-N*PM*sin(angulo_sim-2*pi/3);
                     2*N*Lm*sin(2*(angulo_sim+pi/6+2*pi/3))*var_sim_iabc_ideal(1,i)+2*N*Lm*sin(2*(angulo_sim+pi/6-2*pi/3))*var_sim_iabc_ideal(2,i)-2*N*Lm*sin(2*(angulo_sim+2*pi/3))*var_sim_iabc_ideal(3,i)-N*PM*sin(angulo_sim+2*pi/3)];

    % Evolución del sistema
    di_dt = Matrix_dF_di\(var_sim_vabc(:,i) - R_s*var_sim_iabc_ideal(:,i) - omega_r*Matrix_dF_dx);
    var_sim_iabc_ideal(:,i+1) = var_sim_iabc_ideal(:,i) + Ts*di_dt;
    
    % MODELO EN DERIVADAS x=i
    i_dqz_i(:,i) = T*var_sim_iabc_i(:,i);

    % Flux derivatives
    l_aa = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz_i(1,i),i_dqz_i(2,i),thetara_rad,method);
    l_ab = interpn(idVec,iqVec',xVec,dfluxAdibMatrix,i_dqz_i(1,i),i_dqz_i(2,i),thetara_rad,method);
    l_ac = interpn(idVec,iqVec',xVec,dfluxAdicMatrix,i_dqz_i(1,i),i_dqz_i(2,i),thetara_rad,method);
    
    l_ba = l_ab;
    l_bb = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz_i(1,i),i_dqz_i(2,i),thetarbb_rad,method);
    l_bc = interpn(idVec,iqVec',xVec,dfluxAdibMatrix,i_dqz_i(1,i),i_dqz_i(2,i),thetarcc_rad,method);
    
    l_ca = l_ac;
    l_cb = l_bc;
    l_cc = interpn(idVec,iqVec',xVec,dfluxAdiaMatrix,i_dqz_i(1,i),i_dqz_i(2,i),thetarcc_rad,method);
    
    phi_ar = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz_i(1,i),i_dqz_i(2,i),thetara_rad,method);
    phi_br = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz_i(1,i),i_dqz_i(2,i),thetarbb_rad,method);
    phi_cr = interpn(idVec,iqVec',xVec,dfluxAdxMatrix,i_dqz_i(1,i),i_dqz_i(2,i),thetarcc_rad,method);

    % Matriz de derivadas del flujo respecto a la corriente
    Matrix_dF_di = [l_aa, l_ab, l_ac; l_ba, l_bb, l_bc; l_ca, l_cb, l_cc];
    Matrix_dF_dx = [phi_ar; phi_br; phi_cr];

    % Evolución del sistema
    di_dt = Matrix_dF_di\(var_sim_vabc(:,i) - R_s*var_sim_iabc_i(:,i) - omega_r*Matrix_dF_dx);
    var_sim_iabc_i(:,i+1) = var_sim_iabc_i(:,i) + Ts*di_dt;    

    % MODELO EN DERIVADAS x=flux

%     i_d_k = interp1(fdVec,Id,var_sim_fd(i));
%     i_q_k = interp1(fqVec,Iq,var_sim_fq(i));
    
    var_sim_vdq = T*var_sim_vabc(:,i);

    var_sim_fd(i+1) = var_sim_fd(i) + Ts*(var_sim_vdq(1) - R_s*i_dqz_f(1,i) + N*omega_r*var_sim_fq(i));
    var_sim_fq(i+1) = var_sim_fq(i) + Ts*(var_sim_vdq(2) - R_s*i_dqz_f(2,i) - N*omega_r*var_sim_fd(i));

    i_dqz_f(1,i+1) = interp1(fdVec,Id,var_sim_fd(i+1),method);
    i_dqz_f(2,i+1) = interp1(fqVec,Iq,var_sim_fq(i+1),method);

    var_sim_iabc_f(:,i+1) = invT*[i_dqz_f(1,i+1);i_dqz_f(2,i+1);0];

end

figure
hold on
plot(time,var_sim_iabc_ideal,'*')
plot(time,var_sim_iabc_i)
plot(time,var_sim_iabc_f,'o')
plot(time,ia_sim,'--',time,ib_sim,'--',time,ic_sim,'--')
grid on


