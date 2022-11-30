clear all
clc

%% Lectura del archivo excel
data_table = readtable('Parametrico_100A.csv');
data = table2array(data_table);

%% transformación abc --> dq
N = 5;

for i=1:length(data)
    ia = data(i,2);
    ib = data(i,3);
    ic = -ia-ib;
    thetar_deg = data(i,1);
    thetar = thetar_deg*pi/180;
    phi_a = data(i,4);
    phi_b = data(i,5);
    phi_c = data(i,6);

    T = 2/3*[ cos(N*thetar)    cos(N*thetar - 2*pi/3)    cos(N*thetar + 2*pi/3);
              sin(N*thetar)    sin(N*thetar - 2*pi/3)    sin(N*thetar + 2*pi/3);
              1/2              1/2                       1/2];

    idq0(i,:) = (T*[ia;ib;ic])';
    phi_dq0(i,:) = (T*[phi_a;phi_b;phi_c]);
end

% Montamos una nueva matriz con theta-id-iq-phid-phiq
data_dq = [data(:,1), idq0, phi_dq0];

% Creamos una tabla con los datos
data_dq = round(data_dq,4);
thR = data_dq(:,1);
iD = data_dq(:,2);
iQ = data_dq(:,3);
i0 = data_dq(:,4);
fD = data_dq(:,5);
fQ = data_dq(:,6);
f0 = data_dq(:,7);

data_dq_table = table(thR,iD,iQ,i0,fD,fQ,f0);

% Lo ordenamos por theta-id-iq
data_dq_table_sorted = sortrows(data_dq_table,{'thR','iD','iQ'},{'ascend','ascend','ascend'});
data_dq_sorted = table2array(data_dq_table_sorted);

tamletra = 10;
font = 'Times';

% phi_d para theta=0
fig1 = figure;
set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
scatter3(data_dq_sorted(1:441, 3), data_dq_sorted(1:441, 2), data_dq_sorted(1:441, 5), 'fill');
ylabel('{i_d}','FontSize',tamletra, 'FontName', font)
xlabel('{i_q}','FontSize',tamletra, 'FontName', font)
zlabel('{\it\phi_d}','FontSize',tamletra, 'FontName', font)
view([45 25])
%     zlim([0.05 0.25])
%     xlim([0 800])
grid on
colormap('parula')
set(fig1,'Position',[20 50 360 270])

% phi_q para theta=0
fig1 = figure;
set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
scatter3(data_dq_sorted(1:441, 3), data_dq_sorted(1:441, 2), data_dq_sorted(1:441, 6), 'fill');
ylabel('{i_d}','FontSize',tamletra, 'FontName', font)
xlabel('{i_q}','FontSize',tamletra, 'FontName', font)
zlabel('{\it\phi_q}','FontSize',tamletra, 'FontName', font)
view([45 25])
%     zlim([0.05 0.25])
%     xlim([0 800])
grid on
colormap('parula')
set(fig1,'Position',[20 50 360 270])


%% Creación LUTs phiabc
iaVec = [-100:10:100];
ibVec = iaVec;
thVec = [0:1:72];

nia = length(iaVec);
nib = length(ibVec);
nth = length(thVec);

% pasamos de theta-ia-ib --> ia-ib-theta
fluxA = permute(reshape(data(:,4),[nib,nia,nth]),[2,1,3]);
fluxB = permute(reshape(data(:,5),[nib,nia,nth]),[2,1,3]);
fluxC = permute(reshape(data(:,6),[nib,nia,nth]),[2,1,3]);

%% Representación de los flujos para cada valor de corrientes y de posición
tamletra = 10;
tamlinea = 1;
font = 'Times';
tammarker = 4;

for x=1:10:nth
    % phiA
    fig1 = figure;
    set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
    hold on
    surf(iaVec,ibVec,fluxA(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
    ylabel('{i_a}','FontSize',tamletra, 'FontName', font)
    xlabel('{i_b}','FontSize',tamletra, 'FontName', font)
    zlabel('{\phi_a}','FontSize',tamletra, 'FontName', font)
    view([45 25])
%     zlim([0.05 0.25])
%     xlim([0 800])
    grid
    colormap('parula')
    set(fig1,'Position',[20 50 360 270])

    % phiB
    fig1 = figure;
    set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
    hold on
    surf(iaVec,ibVec,fluxB(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
    ylabel('{i_a}','FontSize',tamletra, 'FontName', font)
    xlabel('{i_b}','FontSize',tamletra, 'FontName', font)
    zlabel('{\phi_b}','FontSize',tamletra, 'FontName', font)
    view([45 25])
%     zlim([0.05 0.25])
%     xlim([0 800])
    grid
    colormap('parula')
    set(fig1,'Position',[20 50 360 270])

    % phiC
    fig1 = figure;
    set(axes, 'Position',[0.16 0.15 0.78 0.815], 'FontSize', tamletra, 'FontName', font);
    hold on
    surf(iaVec,ibVec,fluxC(:,:,x),'FaceColor','interp','LineStyle','none','LineWidth',tamlinea,'FaceAlpha',1)
    ylabel('{i_a}','FontSize',tamletra, 'FontName', font)
    xlabel('{i_b}','FontSize',tamletra, 'FontName', font)
    zlabel('{\phi_c}','FontSize',tamletra, 'FontName', font)
    view([45 25])
%     zlim([0.05 0.25])
%     xlim([0 800])
    grid
    colormap('parula')
    set(fig1,'Position',[20 50 360 270])
end

%% Representación del flujo para un valor de corriente dado
close all;

omegar = 600*2*pi/60; %[rad/s]
N = 5;
omegae = N*omegar; %[rad/s]
fe = omegae/2/pi;
Te = 1/fe;
n = 1000;
thetar_ant = 0;
Ts = Te/n;

time = [0:Ts:Te];
thetar = zeros(1,length(time));

for i=1:length(time)
    ia2(i) = 10*cos(2*pi*fe*time(i));
    ib2(i) = 10*cos(2*pi*fe*time(i)-2*pi/3);
    ic2(i) = 10*cos(2*pi*fe*time(i)+2*pi/3);

    ia(i) = 10*sin(2*pi*fe*time(i));
    ib(i) = 10*sin(2*pi*fe*time(i)-2*pi/3);
    ic(i) = 10*sin(2*pi*fe*time(i)+2*pi/3);

    fa(i) = interpn(iaVec,ibVec',thVec,fluxA,ia(i),ib(i),thetar(i));
    fb(i) = interpn(iaVec,ibVec',thVec,fluxB,ia(i),ib(i),thetar(i));
    fc(i) = interpn(iaVec,ibVec',thVec,fluxC,ia(i),ib(i),thetar(i));
    
    thetae = N*thetar(i)*pi/180;

    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
             -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
              1/2            1/2                      1/2];

    phi_dq(:,i) = T*[fa(i);fb(i);fc(i)];

    thetar(i+1) = thetar_ant + omegar*Ts*180/pi;
    if thetar(i+1) > 360/N
        thetar(i+1) = thetar(i+1) - 360/N;
    end
    thetar_ant = thetar(i+1);
end

figure
hold on
plot(time,ia)
plot(time,ib)
plot(time,ic)
plot(time,ia2,'--')
plot(time,ib2,'--')
plot(time,ic2,'--')
grid on
legend('i_a','i_b','i_c','i_a','i_b','i_c')
xlabel('time [s]')
ylabel('i_{abc} [A]')

figure
hold on
plot(time,fa)
plot(time,fb)
plot(time,fc)
grid on
legend('\phi_a','\phi_b','\phi_c')
xlabel('time [s]')
ylabel('\phi_{abc} [Wb]')

figure
subplot(3,1,1)
plot(time,phi_dq(1,:))
legend('\phi_d')
xlabel('time [s]')
ylabel('\phi [Wb]')
subplot(3,1,2)
plot(time,phi_dq(2,:))
legend('\phi_q')
xlabel('time [s]')
ylabel('\phi [Wb]')
subplot(3,1,3)
plot(time,phi_dq(3,:))
legend('\phi_0')
xlabel('time [s]')
ylabel('\phi [Wb]')

figure
plot(time,thetar(1:end-1))
grid on

%% Flujos para corriente 0
% ia = 0;
ia = iaVec(11);
% ib = 0;
ib = ibVec(11);
ic = -ia-ib;

fphi_dq = zeros(3,nth);

for i = 1:nth
    fa(i) = interpn(iaVec,ibVec',thVec,fluxA,ia,ib,thVec(i));
    fb(i) = interpn(iaVec,ibVec',thVec,fluxB,ia,ib,thVec(i));
    fc(i) = interpn(iaVec,ibVec',thVec,fluxC,ia,ib,thVec(i));
    
    thetae = 5*thVec(i)*pi/180;
    T = 2/3*[ cos(thetae)    cos(thetae - 2*pi/3)    cos(thetae + 2*pi/3);
             -sin(thetae)   -sin(thetae - 2*pi/3)   -sin(thetae + 2*pi/3);
              1/2            1/2                     1/2];

    phi_dq(:,i) = T*[fa(i);fb(i);fc(i)];
end

figure
hold on
plot(thVec,fa)
plot(thVec,fb)
plot(thVec,fc)
grid on
legend('\phi_a','\phi_b','\phi_c')
xlabel('\theta [º]')
ylabel('\phi_{abc} [Wb]')

figure
subplot(3,1,1)
plot(thVec,phi_dq(1,:))
legend('\phi_d')
xlabel('\theta [º]')
ylabel('\phi [Wb]')
subplot(3,1,2)
plot(thVec,phi_dq(2,:))
legend('\phi_q')
xlabel('\theta [º]')
ylabel('\phi [Wb]')
subplot(3,1,3)
plot(thVec,phi_dq(3,:))
legend('\phi_0')
xlabel('\theta [º]')
ylabel('\phi [Wb]')