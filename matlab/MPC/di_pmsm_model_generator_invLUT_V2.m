clear all;
close all;
clc;

%% Parameters initialization

% LUT from ANSYS Maxwell:
% id, iq, thetar --> fluxD, fluxQ, flux0, torque 
ee_ece_table;

% Speed
RPM = 1000;

% Voltage frecuency 
omegae = N*RPM*2*pi/60;

Ts = 1e-6;
tsim = 1;

% LUTs inversas de corrientes de Juan
load('pmsm_phi2i_21p');
% LUTs para Park original
Torque_orig = fliplr(torque);
Id_orig = fliplr(data.i_d);
Iq_orig = fliplr(-data.i_q);

R_s = 0.07;     % Stator resistance [Ohm]
PM = 0.1571;    % Permanent magnet flux
Ld = 0.0013;    % D-axis inductance
Lq = 0.0039;    % Q-axis inductance
L0 = 0.0013;    % Zero-sequence inductance
J = 0.02;       % Inertia

% Pasar de x-y-z a z-x-y en las LUT
Torque_perm = permute(Torque_orig,[3,1,2]);
Id_perm = permute(Id_orig,[3,1,2]);
Iq_perm = permute(Iq_orig,[3,1,2]);
Phid_Vec = data.phi_d;
Phiq_Vec = data.phi_q;
Id_Vec = idVec;
Iq_Vec = iqVec;

% % Nos quemadmos con una parte del la LUT
% Id_Vec = Id_Vec(9:13);
% Iq_Vec = Iq_Vec(9:13);
% 
% Torque_perm = Torque_perm(:,9:13,9:13);

%% Nuevos vectores
index_q = [1, 3, 5, 7, 8, 11, 14, 15, 17, 19, 21];
index_d = [1, 2, 3, 4, 5, 9, 13, 15, 17, 19, 21];

fdVec = Phid_Vec(index_d);
fqVec = Phiq_Vec(index_q);

% Nos quemadmos con una parte del la LUT
Id_LUT = Id_perm(:,index_d,index_q);
Iq_LUT = Iq_perm(:,index_d,index_q);

% save('matlab.mat', '-v7')

%% Nos quedamos únicamente con valores destacados
% Representamos el flujo respecto a la corriente
x = 1; % la dependencia con el ángulo no es importante para definir los punto de la LUT
figure
hold on
for i=1:length(idVec)
    id = idVec(i);
    for q=1:length(iqVec)
        fq(q,i) = interpn(idVec,iqVec',angleVec,fluxQ,id,iqVec(q),angleVec(x),'linear');
    end
    plot(iqVec,fq(:,i))
end
xlabel('i_q')
ylabel('\phi_q')
grid on

figure
hold on
for i=1:length(iqVec)
    iq = iqVec(i);
    for d=1:length(idVec)
        fd(d,i) = interpn(idVec,iqVec',angleVec,fluxD,idVec(d),iq,angleVec(x),'linear');
    end
    plot(idVec,fd(:,i))
end
xlabel('i_d')
ylabel('\phi_d')
grid on

% Representamos la corriente respecto al flujo
x = 1; % la dependencia con el ángulo no es importante para definir los punto de la LUT
figure
hold on
for i=1:length(Phid_Vec)
    Phid = Phid_Vec(i);
    for q=1:length(Phiq_Vec)
        iq(q,i) = interpn(Phid_Vec,Phiq_Vec',angleVec,Iq_orig,Phid,Phiq_Vec(q),angleVec(x),'linear');
    end
    plot(Phiq_Vec,iq(:,i))
end
ylabel('i_q')
xlabel('\phi_q')
grid on
% PUNTOS PARA phi_q --> [-0.3 .. -0.1, 0, 0.1 .. 0.3]

figure
hold on
for i=1:length(Phiq_Vec)
    Phiq = Phiq_Vec(i);
    for d=1:length(Phid_Vec)
        id(d,i) = interpn(Phid_Vec,Phiq_Vec',angleVec,Id_orig,Phid_Vec(d),Phiq,angleVec(x),'linear');
    end
    plot(Phid_Vec,id(:,i))
end
ylabel('i_d')
xlabel('\phi_d')
grid on
% PUNTOS PARA phi_d --> [-0.2 .. -0.1, 0, 0.1 .. 0.3]

figure
hold on
for d=1:length(idVec)
    id = idVec(d);
    for q=1:length(iqVec)
        iq = iqVec(q);
        Te(q,d) = interpn(idVec,iqVec',angleVec,torque,id,iq,angleVec(x),'linear');
    end
    plot(iqVec,Te(:,d))
end
xlabel('i_q')
ylabel('T_e')
grid on

figure
hold on
for q=1:length(iqVec)
    iq = iqVec(q);
    for d=1:length(idVec)
        id = idVec(d);
        Te(d,q) = interpn(idVec,iqVec',angleVec,torque,id,iq,angleVec(x),'linear');
    end
    plot(idVec,Te(:,q))
end
xlabel('i_d')
ylabel('T_e')
grid on

% EL PAR ES COMPLICADO