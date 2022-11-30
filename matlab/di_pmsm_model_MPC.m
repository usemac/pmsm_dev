%% FCS-MPC NON-SINUSOIDAL PMSM

clear all;
close all;
clc;

%% Parámetros del accionameinto PMSM+convertidor

% LUT from ANSYS Maxwell:
% id, iq, thetar --> fluxD, fluxQ, flux0, torque 
ee_ece_table;

% LUTs inversas de corrientes de Juan:
load('pmsm_phi2i_21p');

% Parámetros PMSM
R_s = 0.07;     % Stator resistance [Ohm]
PM = 0.1571;    % Permanent magnet flux
Ld = 0.0013;    % D-axis inductance
Lq = 0.0039;    % Q-axis inductance
L0 = 0.0013;    % Zero-sequence inductance
J = 0.02;       % Inertia
Tmax = 460;     % Maximum torque 
In = 400;       % Nominal current
ang = 2*pi/3;

% Parámetros Convertidor
param_vdc = 300;

param_T = [ 2 -1 -1;
           -1  2 -1;
           -1 -1  2]/3;

% Posibles vectores de disparo considerados por el optimizador
XI8 = [0 0 0;     % 0
       0 0 1;     % 1
       0 1 0;     % 2
       0 1 1;     % 3
       1 0 0;     % 4
       1 0 1;     % 5
       1 1 0;     % 6
       1 1 1];    % 7

param_NV = 8;

%% Definición de parámetros del experimento

% Parámetros controlador PI para el control de velocidad
Kp = 5;
Ki = 20;

% tiempo total del experimento
param_exp_Ttotal = 0.15;  % (s)
% paso de integración
param_sim_h = 1e-6; % (s)
% número de pasos a realizar para simular el experimento
param_sim_Nh = ceil(param_exp_Ttotal/param_sim_h)+1;
% Tiempo de muestreo
cte_contr_Tm_FSMPC = 100*param_sim_h;
% Vector de tiempo
tray_sim_tiempo = [0:param_sim_h:param_exp_Ttotal]';

%% Espacios para históricos

% en pasos de simulación h
tray_sim_idq_real    = zeros( param_sim_Nh+1, 2); % corriente estátor ejes d-q
tray_sim_idq_refe    = zeros( param_sim_Nh+1, 2); % corriente estátor ejes d-q de referencia
tray_sim_ifase_real  = zeros( param_sim_Nh+1, 3); % corriente estátor de fase
tray_sim_ifase_refe  = zeros( param_sim_Nh+1, 3); % corriente estátor de fase de referencia
tray_sim_indopti     = zeros( param_sim_Nh+1, 1); % índice elegido en el último muestreo
tray_sim_velocidad   = zeros( param_sim_Nh+1, 1); % velocidad mecánica wr
tray_sim_thetar      = zeros( param_sim_Nh+1, 1);
tray_sim_thetar_deg  = zeros( param_sim_Nh+1, 1);

% en periodos de muestreo del controlador
tamaprox = 1+ceil( param_sim_Nh*param_sim_h/cte_contr_Tm_FSMPC);
tray_contr_idq_real    = zeros( tamaprox, 2 );
tray_contr_idq_refe    = zeros( tamaprox, 2 );
tray_contr_ifase_real  = zeros( tamaprox, 3 );
tray_contr_ifase_refe  = zeros( tamaprox, 3 );
tray_contr_tiempoc     = zeros( tamaprox, 1 );
tray_contr_indopti     = zeros( tamaprox, 1 ); 
tray_contr_velocidad   = zeros( tamaprox, 1 ); 

%% Estado de partida de la MI para el experimento

% Referencias
% Par de Carga/corriente iq de referencia
% rt = 0.45;  % El valor del par nominal está en el archivo de parámetros m_par_m5f
% TL = rt*Tmax;
isq_ref = 40;

% Corriente id de referencia
isd_ref = 0;  % (A)

% Velocidades
omegar_ref = 1000;  % (rpm)
omegar = omegar_ref*pi/30; % (rad/s)
omegae_ref = N*omegar;  % (rad/s)
omegae = omegae_ref;    % (rad/s)

% Inicializaciones
% valor inicial corrientes estátor d-q
param_exp_idq_inicial = [0 0];
var_sim_x = param_exp_idq_inicial';
var_sim_idq_real = param_exp_idq_inicial;

% Variables control e integrador de ángulo
var_contr_thetar_ant = 0;
var_contr_thetar = 0;
var_contr_thetar_deg_ant = 0;
var_contr_thetar_deg = 0;
var_sim_thetar_ant = 0;
var_sim_thetar = 0;
var_sim_thetar_deg_ant = 0;
var_sim_thetar_deg = 0;
e_km1 = 0;
inte_km1 = 0;
var_contr_fdq = [0;0];
var_sim_fdq = [0;0];

% duración del tiempo de muestreo actual (s)
var_contr_Tm = cte_contr_Tm_FSMPC;
var_contr_Tm_norm = ceil(var_contr_Tm/param_sim_h); % valor normalizado

% cronómetro entre muestreos (s)
var_contr_cronometro = var_contr_Tm;
var_contr_cronometro_norm = var_contr_Tm_norm; % valor normalizado

% contador de periodos de muestreo
var_contr_contadorPM = 0;

% valor inicial acción de control
var_contr_vdq_opt = [0;0];
var_contr_vdq_opt_sig = var_contr_vdq_opt;
tray_contr_indopti(1) = 1;

method = 'linear';

var_contr_vabc = [100*cos(omegae*tray_sim_tiempo), 100*cos(omegae*tray_sim_tiempo-2*pi/3), 100*cos(omegae*tray_sim_tiempo-4*pi/3)];

%% Control de la PMSM con control predictivo
for h=1:param_sim_Nh-1
    
%     % (si se cumple el Tm fijado con anterioridad)
%     if var_contr_cronometro_norm >= var_contr_Tm_norm
%         tiempo_transcurrido = var_contr_cronometro_norm;
%         % puesta a cero del cronómetro
%         var_contr_cronometro_norm = 0;
%             
% %         % CONTROL DE VELOCIDAD (PI)
% %         % corriente de refererencia, salida del PI
% %         werror = omegae_ref - omegae;
% %         inte_k = 0.5*cte_contr_Tm_FSMPC*(werror + e_km1) + inte_km1;
% %         isq_ref = Kp*(werror + Ki*inte_k);
% %         if (isq_ref >= In)
% %             isq_ref = In;
% %         elseif (isq_ref <= -In)
% %             isq_ref = -In;
% %         else
% %             e_km1 = werror;
% %             inte_km1 = inte_k;
% %         end
%         
%         % referencia actual
%         var_contr_idq_refe = [isd_ref isq_ref];
%         % lectura de corrientes dq reales 
%         var_contr_idq_real = var_sim_x;
%         % actuación con la señal de control calculada en el Tm anterior
%         var_contr_vdq = var_contr_vdq_opt_sig;
% 
%         % matriz de transformación de dq a fase
%         param_Mtransf = 2/3*[ cos(N*var_contr_thetar),   cos(N*var_contr_thetar - ang),     cos(N*var_contr_thetar - 2*ang); 
%                               sin(N*var_contr_thetar),   sin(N*var_contr_thetar - ang),     sin(N*var_contr_thetar - 2*ang);
%                                                   1/2,                             1/2,                                1/2];
%         param_Mtransf_inv = inv(param_Mtransf); 
% 
%         % corrientes de fase reales
%         var_contr_ifase_real = param_Mtransf_inv(:,1:2)*var_contr_idq_real;
%         var_contr_ifase_refe = param_Mtransf_inv(:,1:2)*var_contr_idq_refe';
%         
%         % Se le añade RUIDO
% %         var_contr_idq_real(1) = var_contr_idq_real(1)+0.07*randn(1);
% %         var_contr_idq_real(2) = var_contr_idq_real(2)+0.07*randn(1);
%                 
%         % CONTROL DE CORRIENTE FCS-MPC
%         % Predicción para k+1
%         % Flujos en k+1
%         var_contr_fdq_p1p = var_contr_fdq + cte_contr_Tm_FSMPC*(var_contr_vdq - R_s*var_contr_idq_real + omegae*[0 -1;1 0]*var_contr_fdq);
%         % Estimación del ángulo del rotor para k+1 y k+2
%         dthetar_dt_contr = omegar*cte_contr_Tm_FSMPC;
%         var_contr_thetar = var_contr_thetar_ant + dthetar_dt_contr;
%         var_contr_thetar_ant = var_contr_thetar;
%         var_contr_thetar_deg = var_contr_thetar_deg_ant + dthetar_dt_contr*180/pi;
%         if var_contr_thetar_deg > 360/3/N
%             var_contr_thetar_deg = var_contr_thetar_deg - 360/3/N;
%         end
%         var_contr_thetar_deg_ant = var_contr_thetar_deg;
%         var_contr_thetar_degk2 = var_contr_thetar_deg_ant + dthetar_dt_contr*180/pi;
%         if var_contr_thetar_degk2 > 360/3/N
%             var_contr_thetar_degk2 = var_contr_thetar_degk2 - 360/3/N;
%         end
%         % Corrientes en k+1
%         var_contr_idq_p1p(1) = interpn(data.phi_d,data.phi_q',data.theta,data.i_d,var_contr_fdq_p1p(1),var_contr_fdq_p1p(2),var_contr_thetar_deg,method);
%         var_contr_idq_p1p(2) = interpn(data.phi_d,data.phi_q',data.theta,data.i_q,var_contr_fdq_p1p(2),var_contr_fdq_p1p(2),var_contr_thetar_deg,method);
% 
%         % optimización por búsqueda exhaustiva
%         var_contr_Joptimo = inf;
%         var_contr_ioptimo = 0;
%         for vv=1:param_NV
%             % Tensión dq
%             var_contr_vdq = param_vdc*XI8(vv,:)*param_T;
%             % Flujos en k+2
%             var_contr_fdq_p2p = var_contr_fdq_p1p + cte_contr_Tm_FSMPC*(var_contr_vdq(1:2)' - R_s*var_contr_idq_p1p + omegae*[0 -1;1 0]*var_contr_fdq_p1p);
%             % Corrientes en k+1
%             var_contr_idq_p2p(1) = interpn(data.phi_d,data.phi_q',data.theta,data.i_d,var_contr_fdq_p2p(1),var_contr_fdq_p2p(2),var_contr_thetar_degk2,method);
%             var_contr_idq_p2p(2) = interpn(data.phi_d,data.phi_q',data.theta,data.i_q,var_contr_fdq_p2p(2),var_contr_fdq_p2p(2),var_contr_thetar_degk2,method);
%             % discrepancia con la referencia = error de control a 2 pasos
%             var_contr_idq_errc2p = var_contr_idq_refe'-var_contr_idq_p2p;
%             % optimización de la función de coste J
%             var_contr_J = sum(var_contr_idq_errc2p.^2);
%             if var_contr_J < var_contr_Joptimo
%                 var_contr_Joptimo = var_contr_J;
%                 var_contr_ioptimo = vv;
%                 var_contr_vdq_opt = var_contr_vdq(1:2)';
%             end
%         end
%         
%         % señal de control para el siguiente Tm
%         var_contr_vdq_opt_sig = var_contr_vdq_opt;
%         
%         % actualización del contador
%         var_contr_contadorPM = var_contr_contadorPM + 1;
%                         
%         % históricos en periodos de muestreo
%         tray_contr_idq_real(var_contr_contadorPM,:)    = var_contr_idq_real;
%         tray_contr_idq_refe(var_contr_contadorPM,:)    = var_contr_idq_refe;
%         tray_contr_ifase_real(var_contr_contadorPM,:)  = var_contr_ifase_real;
%         tray_contr_ifase_refe(var_contr_contadorPM,:)  = var_contr_ifase_refe;
%         tray_contr_tiempoc(var_contr_contadorPM)       = (h-1)*param_sim_h;
%         tray_contr_indopti(var_contr_contadorPM+1,:)   = var_contr_ioptimo;
%         tray_contr_velocidad(var_contr_contadorPM,:)   = omegar;         
%     end
%    

    % ACTUALIZACIÓN DEL SISTEMA (método Euler)
    % Par en h
    Te = interpn(idVec,iqVec',angleVec,torque,var_sim_x(1),var_sim_x(2),var_sim_thetar_deg,method);
    
    % matriz de transformación de dq a fase
    param_Mtransf = 2/3*[ cos(N*var_sim_thetar),   cos(N*var_sim_thetar - ang),     cos(N*var_sim_thetar - 2*ang); 
                          sin(N*var_sim_thetar),   sin(N*var_sim_thetar - ang),     sin(N*var_sim_thetar - 2*ang);
                                            1/2,                           1/2,                              1/2];
    var_contr_vdq_real = param_Mtransf*var_contr_vabc(h,:)';
    var_contr_vdq_opt(:,h) = var_contr_vdq_real(1:2);

    var_sim_vdq_real = var_contr_vdq_opt(:,h) ;
    
    % Flujos en h+1
    var_sim_fdq = var_sim_fdq + param_sim_h*(var_sim_vdq_real - R_s*var_sim_x + omegae*[0 -1;1 0]*var_sim_fdq);

    % Estimación del ángulo del rotor para h+1
    dthetar_dt = omegar*param_sim_h;
    var_sim_thetar = var_sim_thetar_ant + dthetar_dt;
    var_sim_thetar_ant = var_sim_thetar;
    var_sim_thetar_deg = var_sim_thetar_deg_ant + dthetar_dt*180/pi;
    if var_sim_thetar_deg > 360/3/N
        var_sim_thetar_deg = var_sim_thetar_deg - 360/3/N;
    end
    var_sim_thetar_deg_ant = var_sim_thetar_deg;
    
    % Corrientes en h+1
    var_sim_x(1) = interpn(data.phi_d,data.phi_q',data.theta,data.i_d,var_sim_fdq(1),var_sim_fdq(2),var_sim_thetar_deg,method);
    var_sim_x(2) = interpn(data.phi_d,data.phi_q',data.theta,data.i_q,var_sim_fdq(1),var_sim_fdq(2),var_sim_thetar_deg,method);
    var_sim_idq_real = var_sim_x';

    % Velocidad en h+1
%     omegar = omegar + (param_sim_h/param_J)*(Te-TL);
%     omegae = N*omegar;
    
    % matriz de transformación de dq a fase
    param_Mtransf = 2/3*[ cos(N*var_sim_thetar),   cos(N*var_sim_thetar - ang),     cos(N*var_sim_thetar - 2*ang); 
                          sin(N*var_sim_thetar),   sin(N*var_sim_thetar - ang),     sin(N*var_sim_thetar - 2*ang);
                                            1/2,                           1/2,                              1/2];
    param_Mtransf_inv = inv(param_Mtransf); 

    % corrientes de fase en h+1
    var_sim_ifase_real = param_Mtransf_inv(:,1:2)*var_sim_x;
    
    % actualización variables
    var_contr_cronometro_norm = var_contr_cronometro_norm + 1;
    
    % históricos en pasos de simulación
    tray_sim_idq_real(h+1,:)    = var_sim_idq_real;
%     tray_sim_idq_refe(h ,:)     = var_contr_idq_refe;
    tray_sim_ifase_real(h+1,:)  = var_sim_ifase_real;
%     tray_sim_ifase_refe(h,:)    = var_contr_ifase_refe;
%     tray_sim_indopti(h ,:)      = var_contr_ioptimo;
    tray_sim_velocidad(h+1,:)   = omegar;
    tray_sim_thetar(h+1,:)      = var_sim_thetar;
    tray_sim_thetar_deg(h+1,:)  = var_sim_thetar_deg;
    tray_sim_torque(h)          = Te;
end

%% se eliminan muestras no usadas de históricos cada Tm
if var_contr_contadorPM < tamaprox
    tray_contr_idq_real( var_contr_contadorPM+1:end,:)    = [];
    tray_contr_idq_refe( var_contr_contadorPM+1:end,:)    = [];
    tray_contr_ifase_real( var_contr_contadorPM+1:end,:)  = [];
    tray_contr_ifase_refe( var_contr_contadorPM+1:end,:)  = [];
    tray_contr_indopti( var_contr_contadorPM+1:end,:)     = [];
    tray_contr_tiempoc( var_contr_contadorPM+1:end,:)     = [];
    tray_contr_velocidad( var_contr_contadorPM+1:end,:)   = []; 
end

%% Representación de resultados
figure
hold on
plot(tray_sim_tiempo,tray_sim_idq_real(1:end-1,:))
grid on

figure
hold on
plot(tray_contr_tiempoc,tray_contr_idq_real)
plot(tray_contr_tiempoc,tray_contr_idq_refe,'k--')
grid on

figure
hold on
plot(tray_contr_tiempoc,tray_contr_ifase_real)
plot(tray_contr_tiempoc,tray_contr_ifase_refe,'--')
grid on