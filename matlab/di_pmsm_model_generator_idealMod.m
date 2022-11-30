%% EL OBJETIVO DE ESTE SCRIPT ES MODELAR LA PMSM DE FORMA IDEAL
% Consideraremos las ecuaciones ideales, para lo cual hay que obtener Ld,
% Lq y PM de las LUTs.
% El objetivo es poder montar un control FOC.

clear all;
close all;
clc;

%% Parameters initialization

% LUT from ANSYS Maxwell:
% id, iq, thetar --> fluxD, fluxQ, flux0, torque 
ee_ece_table;

% Stator resistance [Ohm]
R_s = 0.07;

%% Búsqueda de Ld, Lq y PM
% PM: flujoD para id=iq=0
id = 0;
iq = 0;

metodo = 'linear';

for x=1:length(angleVec)
    fd_0(x) = interpn(idVec,iqVec',angleVec,fluxD,id,iq,angleVec(x),metodo);
    fq_0(x) = interpn(idVec,iqVec',angleVec,fluxQ,id,iq,angleVec(x),metodo);
end

fd_0mean = mean(fd_0);
fq_0mean = mean(fq_0);

PM = fd_0;
PM_mean = fd_0mean;

figure
plot(angleVec,PM)
xlabel('\theta_r [deg]')
ylabel('\phi_m')
grid on


% Para obtener Ld y Lq vamos a usar las ecuaciones:
% fluxd = Ld*id + PM
% fluxq = Lq*iq

metodo = 'linear';
tamletra = 10;
tamlinea = 1;
font = 'Times';
tammarker = 4;

% Lq: fq(id=0)/iq
% Representamos el flujo respecto a la corriente
id = 0;

figure
hold on
for x=1:length(angleVec)
    for q=1:length(iqVec)
        fq(q,x) = interpn(idVec,iqVec',angleVec,fluxQ,id,iqVec(q),angleVec(x),metodo);
%         Lq_matrix(q,x) = fq(q,x)/iqVec(q); %Con q=1:10
    end
%     Lq_mean(x) = mean(Lq_matrix(:,x)); %Con q=1:10
    plot(iqVec,fq(:,x))
end
% Lq = mean(Lq_mean); %Con q=1:10
xlabel('i_q')
ylabel('\phi_q')
grid on

% Podemos considerar que Lq es igual a la pendiente de la cruva en la zona
% lineal
m30_idx = find(iqVec==-30);
p30_idx = find(iqVec==30);
Lq = (mean(fq(p30_idx,:))-mean(fq(m30_idx,:)))/(30-(-30))


% Ld: (fd(iq=0)-PM)/id
% Tenemos dos opciones:
% Opción 1: Restamos el flujo de los imanes al flujod y calculamos la
% pendiente
iq = 0;

figure
hold on
for x=1:length(angleVec)
    fM = PM(x);
    for d=1:length(idVec)
        fd(d,x) = interpn(idVec,iqVec',angleVec,fluxD,idVec(d),iq,angleVec(x),metodo);
        dflux(d,x) = fd(d,x)-fM;
%         Ld_matrix(d,x) = (fd-fM)/idVec(d); %Con d=1:10
    end
%     Ld_mean(x) = mean(Ld_matrix(:,x)); %Con d=1:10
    plot(idVec,dflux(:,x))
end
% Ld = mean(Ld_mean) %Con d=1:10
xlabel('i_d')
ylabel('\phi_d-\phi_m')
grid on

% Calculamos la pendiente en la zona lineal
m90_idx = find(idVec==-90);
p30_idx = find(idVec==30);
Ld = (mean(dflux(p30_idx,:))-mean(dflux(m90_idx,:)))/(30-(-90))

% Opción 2: directamente representamos el flujo respecto a la corriente
iq = 0;

figure
hold on
for x=1:length(angleVec)
    fM = PM(x);
    for d=1:length(idVec)
        fd(d,x) = interpn(idVec,iqVec',angleVec,fluxD,idVec(d),iq,angleVec(x),metodo);
%         Ld_matrix(d,x) = (fd-fM)/idVec(d); %Con d=1:10
    end
%     Ld_mean(x) = mean(Ld_matrix(:,x)); %Con d=1:10
    plot(idVec,fd(:,x))
end
% Ld = mean(Ld_mean) %Con d=1:10
xlabel('i_d')
ylabel('\phi_d')
grid on

% Calculamos la pendiente en la zona lineal
m90_idx = find(idVec==-90);
p60_idx = find(idVec==60);
Ld = (mean(dflux(p60_idx,:))-mean(dflux(m90_idx,:)))/(60-(-90))

% SALE LO MISMO CON AMBAS OPCIONES

%% Simulación de la máquina ideal con los parámetros calculados (SIMULINK)
L0 = Ld;

xVec = pi/180*linspace(0,360/N,180/N+1); % Rotor angle vector [rad]
[F,T,dFdA,dFdB,dFdC,dFdX] = ee_generateIdealPMSMfluxData(PM_mean,Ld,Lq,L0,idVec,iqVec,xVec);

Ts = 1e-6;
tsim = 0.01;

Rload = 5;
RPM = 7200;

sim('ee_import_fem_maxwell_idealMod')

%% SIMULACIÓN DE LAS MÁQUINAS EN MATLAB DE FORMA DISCRETA
% Comprobación de las LUTs del modelo no-sinusoidal
time = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.angular_velocity.series.time;
thetar_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.angular_position.series.values;     % [deg] rotor angle, from wr
omegar_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.angular_velocity.series.values;     % [rpm] rotor mechanical speed wr
torque_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.electrical_torque.series.values;
i_a_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.i_a.series.values;
i_b_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.i_b.series.values;
i_c_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.i_c.series.values;
i_d_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.i_d.series.values;
i_q_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.i_q.series.values;
v_a_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.v_a.series.values;
v_b_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.v_b.series.values;
v_c_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.v_c.series.values;
load_torque_sim = simlog_ee_import_fem_maxwell.FEM_Parameterized_PMSM2.torque.series.values;

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

% Comprobación de corrientes abc y de par
for i = 1:nd
    thetae = N*thetar_sim(i)*pi/180;
    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
              sin(thetae)    sin(thetae - 2*pi/3)    sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];
    invT = inv(T);
    iabc(:,i) = invT(:,1:2)*[i_d_sim(i);i_q_sim(i)];
    idq(:,i) = T*[i_a_sim(i);i_b_sim(i);i_c_sim(i)];
    torque_cal(i) = interpn(idVec,iqVec',angleVec,torque,idq(1,i),idq(2,i),thetar_sim_nuevo(i),"linear");
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
