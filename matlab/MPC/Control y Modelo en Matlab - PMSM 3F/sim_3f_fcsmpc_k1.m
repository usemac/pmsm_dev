% C0NTROL FCSMPC PARA PMSM DE 3 FASES
% CREADO POR CMT, 31/03/2022

clear all
close all
clc

%% definición de parámetros del experimento

% global param_Ld1 param_flux_max_1 param_Lq1 param_Ld3 param_flux_max_3 param_Lq3
% global param_p Vmax_ph_ph Imax_ph
% global omegam_ref Tem_ref_real

% Load torque
TL = 10;
% Velocidad mecánica de referencia
omegam_ref = 1000*2*pi/60;  % rad/s
% % Par de referencia
% Tem_ref = 8;     %Nm
% Load of the curve speed-torque for the saturation
% % caso 3) h1+h3
% w_Tem_val_curve_h1h3        % Load of the curve speed-torque (saturation)
% w_curve = w_Tem_curve(:,1); % Speed values of the curve
% Tem_max = w_Tem_curve(:,2); % Torque values of the curve
% % caso 4) h1+h3+h5
% w_Tem_val_curve_h1h3h5        % Load of the curve speed-torque (saturation)
% w_curve = w_Tem_curve(:,1); % Speed values of the curve
% Tem_max = w_Tem_curve(:,2); % Torque values of the curve
% Load of the Look-up tables for current references
% % caso 3) h1+h3
% LUT_id1_h1h3                  % Load of the id1_ref values of the LUT
% LUT_iq1_h1h3                  % Load of the iq1_ref values of the LUT
% LUT_id3_h1h3                  % Load of the id3_ref values of the LUT
% LUT_iq3_h1h3                  % Load of the iq3_ref values of the LUT
% LUT_w_val_h1h3                % Load of speed values of the LUT (rows)
% LUT_Tem_val_h1h3              % Load of torque values of the LUT (columns)
% % caso 4) h1+h3+h5
% LUT_id1_h1h3h5                  % Load of the id1_ref values of the LUT
% LUT_iq1_h1h3h5                  % Load of the iq1_ref values of the LUT
% LUT_id3_h1h3h5                  % Load of the id3_ref values of the LUT
% LUT_iq3_h1h3h5                  % Load of the iq3_ref values of the LUT
% LUT_w_val_h1h3h5                % Load of speed values of the LUT (rows)
% LUT_Tem_val_h1h3h5              % Load of torque values of the LUT (columns)
% mat_w_val = ones(length(w_val),length(Tem_val));
% for i=1:length(Tem_val)
%     mat_w_val(:,i) = w_val';
% end
% mat_Tem_val = ones(length(w_val),length(Tem_val));
% for i=1:length(w_val)
%     mat_Tem_val(i,:) = Tem_val;
% end

% tiempo total del experimento
param_exp_Ttotal = 2;


%% definición de parámetros del simulador (integración ecs. dif.)

% paso de integración
param_sim_h = 1e-6; % (s)
% número de pasos a realizar para simular el experimento
param_sim_Nh = ceil(param_exp_Ttotal/param_sim_h)+1;

%% Definición de parámetros del controlador

% Tiempo de muestreo del controlador
cte_contr_Tm_FSMPC = 50*param_sim_h;
% % Tiempo de muestreo del optimizador
% cte_opt_Tm = 3000*param_sim_h;
% cte_opt_Tm = cte_contr_Tm_FSMPC;

%% definición de la PMSM 3F

% definición de la PMSM de 3F
param_fcsmpc;

%% definiciones y inicios de variables

tray_sim_tiempo = [0:param_sim_h:param_exp_Ttotal]';

%% espacios para históricos

% en pasos de simulación h
tray_sim_ifase_real      = zeros( param_sim_Nh, 3); 
tray_sim_ifase_refe      = zeros( param_sim_Nh, 3);
tray_sim_idq_real        = zeros( param_sim_Nh+1, 2);
tray_sim_idq_refe        = zeros( param_sim_Nh+1, 2);
tray_sim_bemf_fase       = zeros( param_sim_Nh, 3);
tray_sim_bemf_dq         = zeros( param_sim_Nh, 2);
tray_sim_indopti         = zeros( param_sim_Nh, 1); % índice elegido en el último muestreo
tray_sim_velocidad       = zeros( param_sim_Nh+1, 1); % índice elegido en el último muestreo
tray_sim_velocidad_refe  = zeros( param_sim_Nh, 1);
tray_sim_par_real        = zeros( param_sim_Nh+1, 1);
tray_sim_par_refe        = zeros( param_sim_Nh, 1);
tray_sim_vfase_real      = zeros( param_sim_Nh, 3);
tray_sim_vdq_real        = zeros( param_sim_Nh, 2);
tray_sim_theta           = zeros( param_sim_Nh+1, 1);
tray_sim_vfase_vsi       = zeros( param_sim_Nh, 3);

% % en periodos de muestreo del controlador
% tamaprox = 1+ceil( param_sim_Nh*param_sim_h/cte_contr_Tm_FSMPC);
% tray_contr_is_real    = zeros( tamaprox, 4 );
% tray_contr_is_refe    = zeros( tamaprox, 4 );
% tray_contr_idq_real   = zeros( tamaprox, 2 );
% tray_contr_idq_refe   = zeros( tamaprox, 2 );
% tray_contr_tiempoc    = zeros( tamaprox, 1 );
% tray_contr_indopti    = zeros( tamaprox, 1 );  % se almacena el índice nada más
% tray_contr_ir_obsv    = zeros( tamaprox, 2 );
% tray_contr_ir_real    = zeros( tamaprox, 2 );
% tray_contr_is_obsv    = zeros( tamaprox, 4 );
% tray_contr_ifase_real = zeros( tamaprox, 5 );
% tray_contr_ifase_refe = zeros( tamaprox, 5 );
% tray_contr_velocidad  = zeros( tamaprox, 1 ); 
% tray_contr_welect     = zeros( tamaprox, 1 ); 

%% estado de partida de la MI para el experimento

% Velocidades iniciales
omegam = omegam_ref;
var_sim_x = [0 0]';
% var_sim_idq_real = param_exp_idq_inicial';
theta_ant = 0;
theta = 0;
% theta_contr_ant = 0;
% theta_contr = 0;
% theta_contr_p = 0;
e_km1 = 0;
inte_km1 = 0;
isd_ref = 0;

% isd_ref = -120;
% isq_ref = 12;

% duración del tiempo de muestreo del controlador actual 
var_contr_Tm = cte_contr_Tm_FSMPC;
var_contr_Tm_norm = ceil(var_contr_Tm/param_sim_h);
% % duración del tiempo de muestreo del optimizador actual 
% var_opt_Tm = cte_opt_Tm;
% var_opt_Tm_norm = ceil(var_opt_Tm/cte_contr_Tm_FSMPC);

% cronómetro entre muestreos de controlador
var_contr_cronometro = var_contr_Tm;
var_contr_cronometro_norm = var_contr_Tm_norm;
% % cronómetro entre muestreos de optimización
% var_opt_cronometro = var_opt_Tm;
% var_opt_cronometro_norm = var_opt_Tm_norm;

% contador de periodos de muestreo y tiempo transcurrido
var_contr_contadorPM = 0;
var_opt_contadorPM = 0;
tiempo_real = 0;

% valor inicial acción de control
var_contr_sopt = XI8(1,:);
var_contr_sopt_sig = var_contr_sopt;

% % Seed and options for 'fmincon'
% i0 = [-10 ; 10 ; 10 ; -10];
% options = optimset('Algorithm','sqp','MaxFunEvals',1000,'MaxIter',1000);
% % Declaration of linear constraints to use in 'fmincon':
% % A * i <= b
% % Aeq * i = beq
% % lb <= i <= ub
% A   = [];
% b   = [];
% Aeq = [];
% beq = [];
% lb  = [];
% ub  = [];inicializaciones para el proceso de optimización
% % 

% referencias optimas analíticas sin límites
% id1_ref_real = 0;
% iq1_ref_real = Tem_ref*epsilon_q1/(epsilon_q1^2+epsilon_q3^2);
% id3_ref_real = 0;
% iq3_ref_real = Tem_ref*epsilon_q3/(epsilon_q1^2+epsilon_q3^2);    


%% Experimento de control de la M6F con el controlador por persecución
for h=1:param_sim_Nh
    
    % (si se cumple el Tm fijado con anterioridad)
    if var_contr_cronometro_norm >= var_contr_Tm_norm
        % puesta a cero del cronómetro
        var_contr_cronometro_norm = 0;
        
%         % Actualización de la velocidad de referencia
%         if (tiempo_real>=0 && tiempo_real<0.6)
%             omegam_ref = 50;
%         elseif (tiempo_real>=0.6 && tiempo_real<0.7)
%             omegam_ref = 1000*tiempo_real - 550;   
%         elseif (tiempo_real>=0.7 && tiempo_real<1.1)
%             omegam_ref = 150;
%         elseif (tiempo_real>=1.1 && tiempo_real<1.2)
%             omegam_ref = 500*tiempo_real - 400;
%         elseif (tiempo_real>=1.2)
%             omegam_ref = 200;
%         end

%         omegam_ref = (240/1.5)*tiempo_real + 1;
                
        % CONTROL DE VELOCIDAD (PI)
        % corriente de refererencia, salida del PI
        werror = omegam_ref - omegam;
        inte_k = 0.5*cte_contr_Tm_FSMPC*(werror + e_km1) + inte_km1;
        isq_ref = Kp_w*(werror + Ki_w*inte_k);
        if (isq_ref >= In) 
            isq_ref = In;
        elseif (isq_ref <= -In)
            isq_ref = -In;
        else
            e_km1 = werror;
            inte_km1 = inte_k;
        end

%         omegam = omegam_ref;

%         % obtención del par de referencia a partir de la curva Te-w
%         Tem_ref_max = interp1(w_curve,Tem_max,omegam_ref,'linear');
%         if Tem_ref_max >= Tem_ref
%             Tem_ref_real = Tem_ref;
%         else
%             Tem_ref_real = Tem_ref_max;
%         end
%         
%         % actualizamos la referencia de corrientes si ha pasado el tiempo
%         % Tm del optimizador
%         if var_opt_cronometro_norm >= var_opt_Tm_norm
%             % puesta a cero del cronómetro
%             var_opt_cronometro_norm = 0;
%             
%             % Aplicación de los límites
%             % usando la función de optimización
%             [sol,loss,flag] = fmincon(@id_iq_function_5ph_phph,i0,A,b,Aeq,beq,lb,ub,@id_iq_constraint_5ph_phph,options);
%             id1_ref_real = sol(1);
%             iq1_ref_real = sol(2);
%             id3_ref_real = sol(3);
%             iq3_ref_real = sol(4);   
%             var_contr_idq_refe = [id1_ref_real, iq1_ref_real, id3_ref_real, iq3_ref_real];
% %             % usando las lookup tables precalculadas
% %             id1_ref_real = interp2(mat_Tem_val,mat_w_val,id1_ref,Tem_ref_real,omegam_ref);
% %             iq1_ref_real = interp2(mat_Tem_val,mat_w_val,iq1_ref,Tem_ref_real,omegam_ref);
% %             id3_ref_real = interp2(mat_Tem_val,mat_w_val,id3_ref,Tem_ref_real,omegam_ref);
% %             iq3_ref_real = interp2(mat_Tem_val,mat_w_val,iq3_ref,Tem_ref_real,omegam_ref);      
% %             var_contr_idq_refe = [id1_ref_real, iq1_ref_real, id3_ref_real, iq3_ref_real];
%             
%             % actualización del contador
%             var_opt_contadorPM = var_opt_contadorPM + 1;
%             
%         end
%         var_opt_cronometro_norm = var_opt_cronometro_norm + 1;
        
        % Referencia actual
        var_contr_idq_refe = [isd_ref, -isq_ref];
        % lectura de corrientes dq reales 
        var_contr_idq_real = var_sim_x;
        
        % matriz de transformación de dq a fase
%         param_contr_Mtransf = 2/3*[ cos(param_p*theta),   cos(param_p*theta - ang),     cos(param_p*theta - 2*ang); 
%                                    -sin(param_p*theta),  -sin(param_p*theta - ang),    -sin(param_p*theta - 2*ang);
%                                                    1/2,                        1/2,                            1/2];   
        param_contr_Mtransf = 2/3*[ cos(param_p*theta),   cos(param_p*theta - ang),     cos(param_p*theta - 2*ang); 
                                    sin(param_p*theta),   sin(param_p*theta - ang),     sin(param_p*theta - 2*ang);
                                                   1/2,                        1/2,                            1/2];  
        
        % CONTROL DE CORRIENTE FCS-MPC a 1 paso
        % optimización por búsqueda exhaustiva
        var_contr_Joptimo = inf;
        var_contr_ioptimo = 0;
        for vv=1:param_NV
            % Predicción para k+1
            % vector de voltaje en fase real
            var_opt_vfase = Vbus*param_T*XI8(vv,:)';
            % vector de voltaje d-q
            var_opt_vdq = param_contr_Mtransf(1:2,:)*var_opt_vfase;
            % cálculo de la bemf
            var_contr_e = omegam*param_PMSM_D*var_contr_idq_real + omegam*param_PMSM_E;
            % cálculo de las corrientes d-q para h+1
            var_contr_idq_p1p = var_contr_idq_real + cte_contr_Tm_FSMPC*(param_PMSM_A*var_opt_vdq + param_PMSM_B*var_contr_idq_real + param_PMSM_C*var_contr_e);
            % discrepancia con la referencia = error de control a 2 pasos
            var_contr_idq_err = var_contr_idq_refe'-var_contr_idq_p1p;
            % optimización de la función de coste J
            var_contr_J = sum(var_contr_idq_err.^2);
            if var_contr_J < var_contr_Joptimo
                var_contr_Joptimo = var_contr_J;
                var_contr_ioptimo = vv;
                var_contr_vdq_opt = var_opt_vdq;
            end
        end
               
        % actualización del contador
        var_contr_contadorPM = var_contr_contadorPM + 1;
        
%         % Estimación del ángulo del rotor
%         theta_contr = theta_contr_ant + omegam_ref*cte_contr_Tm_FSMPC;
%         theta_contr_p = theta_contr_ant + 3*omegam_ref*cte_contr_Tm_FSMPC;
%         theta_contr_ant = theta_contr;
%                 
%         % históricos en periodos de muestreo (ojo Tm no constante)
%         tray_contr_is_real( var_contr_contadorPM,:)   = var_contr_is_real;
%         tray_contr_is_refe( var_contr_contadorPM,:)   = var_contr_is_refe;
%         tray_contr_idq_real(var_contr_contadorPM,:)   = var_contr_idq_real;
%         tray_contr_idq_refe(var_contr_contadorPM,:)   = var_contr_idq_refe;
%         tray_contr_ifase_real(var_contr_contadorPM,:) = var_contr_is_real'*param_Minv;
%         tray_contr_ifase_refe(var_contr_contadorPM,:) = var_contr_is_refe*param_Minv;
%         tray_contr_tiempoc( var_contr_contadorPM)     = (h-1)*param_sim_h;
%         tray_contr_indopti( var_contr_contadorPM+1,:) = var_contr_ioptimo;
%         tray_contr_is_obsv( var_contr_contadorPM,:)   = var_obs_z(1:4);
%         tray_contr_ir_obsv( var_contr_contadorPM,:)   = var_obs_z(5:6);
%         tray_contr_ir_real( var_contr_contadorPM,:)   = var_sim_ir_real;
%         tray_contr_velocidad(var_contr_contadorPM,:)  = omegam; 
%         tray_contr_welect(var_contr_contadorPM,:)     = we; 
        
    end

    % matriz de transformación de dq a fase
%     param_Mtransf = 2/3*[ cos(param_p*theta),   cos(param_p*theta - ang),     cos(param_p*theta - 2*ang); 
%                          -sin(param_p*theta),  -sin(param_p*theta - ang),    -sin(param_p*theta - 2*ang);
%                                          1/2,                        1/2,                           1/2];
    param_Mtransf = 2/3*[ cos(param_p*theta),   cos(param_p*theta - ang),     cos(param_p*theta - 2*ang); 
                          sin(param_p*theta),   sin(param_p*theta - ang),     sin(param_p*theta - 2*ang);
                                         1/2,                        1/2,                           1/2];
    param_Mtransf_inv = inv(param_Mtransf); 
    % corrientes de fase en h
    var_sim_ifase_real = param_Mtransf_inv(:,1:2)*var_sim_x;
    % corrientes de referencia de fase en h
    var_sim_ifase_refe = param_Mtransf_inv(:,1:2)*var_contr_idq_refe';
    
    % voltaje optimo a aplicar  
    % en d-q
    var_sim_vdq_real = var_contr_vdq_opt;
    % en fase
    var_sim_vfase_real = param_Mtransf_inv(:,1:2)*var_sim_vdq_real;
       
    % actualización del estado de la máquina (método Euler)
    % cálculo bemf d-q
    var_sim_e = omegam*param_PMSM_D*var_sim_x + omegam*param_PMSM_E;
    % cálculo bemf fase
    var_sim_efase = param_Mtransf_inv(:,1:2)*var_sim_e;
    % cálculo de las corrientes d-q para h+1
    var_sim_x = var_sim_x + param_sim_h*(param_PMSM_A*var_sim_vdq_real + param_PMSM_B*var_sim_x + param_PMSM_C*var_sim_e);
    
    % Ecuaciones mecánicas
%     Tem_real =   3/2*param_p*(param_Ld-param_Lq)*var_sim_x(1)*var_sim_x(2) + ...
%                  3/2*param_p*param_flux_m*var_sim_x(2); 
    Tem_real =   3/2*param_p*(-param_Ld+param_Lq)*var_sim_x(1)*var_sim_x(2) - ...
                 3/2*param_p*param_flux_m*var_sim_x(2); 
    omegam = omegam + (param_sim_h/param_J)*(Tem_real-TL-param_f*omegam);
    
    % Estimación del ángulo del rotor
    theta = theta_ant + omegam*param_sim_h;
    theta_ant = theta;
    
    % actualización variables
    var_contr_cronometro_norm = var_contr_cronometro_norm + 1;
    tiempo_real = tiempo_real + param_sim_h;    
     
    % históricos en pasos de simulación
    tray_sim_ifase_real(h,:)     = var_sim_ifase_real';
    tray_sim_ifase_refe(h,:)     = var_sim_ifase_refe';
    tray_sim_idq_real(h+1,:)     = var_sim_x';
    tray_sim_idq_refe(h,:)       = var_contr_idq_refe;
    tray_sim_bemf_fase(h,:)      = var_sim_efase';
    tray_sim_bemf_dq(h,:)        = var_sim_e';
    tray_sim_indopti(h,:)        = var_contr_ioptimo;
    tray_sim_velocidad(h+1,:)    = omegam;
    tray_sim_velocidad_refe(h,:) = omegam_ref;
    tray_sim_par_real(h+1,:)     = Tem_real;
    tray_sim_par_refe(h,:)       = TL;
    tray_sim_vfase_real(h,:)     = var_sim_vfase_real';
    tray_sim_vdq_real(h,:)       = var_sim_vdq_real';
    tray_sim_theta(h+1,:)        = theta;
%     tray_sim_vfase_vsi(h,:)    = var_sim_vvsi';
%     tray_sim_theta_contr(h+1,:)= theta_contr;
end

% % se eliminan muestras no usadas de históricos cada Tm
% if var_contr_contadorPM < tamaprox
%     tray_contr_is_real( var_contr_contadorPM+1:end,:)    = [];
%     tray_contr_is_refe( var_contr_contadorPM+1:end,:)    = [];
%     tray_contr_idq_refe( var_contr_contadorPM+1:end,:)   = [];
%     tray_contr_idq_real( var_contr_contadorPM+1:end,:)   = [];
%     tray_contr_ifase_real( var_contr_contadorPM+1:end,:) = [];
%     tray_contr_ifase_refe( var_contr_contadorPM+1:end,:) = [];
%     tray_contr_indopti( var_contr_contadorPM+1:end,:)    = [];
%     tray_contr_tiempoc( var_contr_contadorPM+1:end,:)    = [];
%     tray_contr_ir_obsv( var_contr_contadorPM+1:end,:)    = [];
%     tray_contr_ir_real( var_contr_contadorPM+1:end,:)    = [];
%     tray_contr_is_obsv( var_contr_contadorPM+1:end,:)    = [];
%     tray_contr_velocidad( var_contr_contadorPM+1:end,:)  = []; 
%     tray_contr_welect( var_contr_contadorPM+1:end,:)     = [];
% end

%% gráficas
figure
plot(tray_sim_tiempo,tray_sim_idq_real(1:end-1,:))
hold on
plot(tray_sim_tiempo,tray_sim_idq_refe(1:end-1,:),'k--')

figure
plot(tray_sim_tiempo,tray_sim_ifase_real)
hold on
plot(tray_sim_tiempo,tray_sim_ifase_refe,'k--')
% plot([0,tray_sim_tiempo(end)],[Imax_ph,Imax_ph],'--k')
% plot([0,tray_sim_tiempo(end)],[-Imax_ph,-Imax_ph],'--k')

% figure
% plot(tray_sim_tiempo,tray_sim_vdq_real)
% 
% figure
% plot(tray_sim_tiempo,tray_sim_vfase_real)

% % Tensiones de fase filtradas
% Va_f = zeros(length(tray_sim_tiempo),1);
% Vb_f = zeros(length(tray_sim_tiempo),1);
% Vc_f = zeros(length(tray_sim_tiempo),1);
% Vd_f = zeros(length(tray_sim_tiempo),1);
% Ve_f = zeros(length(tray_sim_tiempo),1);
% for i=2:length(tray_sim_tiempo)
%     Va_f(i) = Va_f(i-1) + 10000*param_sim_h*(tray_sim_vfase_real(i-1,1)-Va_f(i-1));
%     Vb_f(i) = Vb_f(i-1) + 10000*param_sim_h*(tray_sim_vfase_real(i-1,2)-Vb_f(i-1));
%     Vc_f(i) = Vc_f(i-1) + 10000*param_sim_h*(tray_sim_vfase_real(i-1,3)-Vc_f(i-1));
%     Vd_f(i) = Vd_f(i-1) + 10000*param_sim_h*(tray_sim_vfase_real(i-1,4)-Vd_f(i-1));
%     Ve_f(i) = Ve_f(i-1) + 10000*param_sim_h*(tray_sim_vfase_real(i-1,5)-Ve_f(i-1));
% end
% tray_sim_vfase_real_filt = [Va_f, Vb_f, Vc_f, Vd_f, Ve_f];
% figure
% % plot(tray_sim_tiempo,tray_sim_vfase_real)
% hold on
% plot(tray_sim_tiempo,tray_sim_vfase_real_filt)
% plot([0,tray_sim_tiempo(end)],[Vbus/2,Vbus/2],'--k')
% plot([0,tray_sim_tiempo(end)],[-Vbus/2,-Vbus/2],'--k')

% Prueba de filtrado de la tensión a partir del modelo de la máquina
Vfase_f = zeros(3,length(tray_sim_tiempo));
for i=2:length(tray_sim_tiempo)
    bemf = tray_sim_velocidad(i)*param_PMSM_D*tray_sim_idq_real(i,:)' + tray_sim_velocidad(i)*param_PMSM_E;
    Vdq_f = inv(param_PMSM_A)*(tray_sim_idq_real(i,:)' - tray_sim_idq_real(i-1,:)' - param_PMSM_B*tray_sim_idq_real(i,:)' - param_PMSM_C*bemf);
    theta = tray_sim_theta(i);
    param_Mtransf = 2/3*[ cos(param_p*theta),   cos(param_p*theta - ang),     cos(param_p*theta - 2*ang); 
                         -sin(param_p*theta),  -sin(param_p*theta - ang),    -sin(param_p*theta - 2*ang);
                         1/2,                        1/2,                            1/2];  
    param_Mtransf_inv = inv(param_Mtransf); 
    Vfase_f(:,i) = param_Mtransf_inv(:,1:2)*Vdq_f;
end
% figure
% % plot(tray_sim_tiempo,tray_sim_vfase_real)
% hold on
% plot(tray_sim_tiempo,Vfase_f)
% plot([0,tray_sim_tiempo(end)],[Vbus/2,Vbus/2],'--k')
% plot([0,tray_sim_tiempo(end)],[-Vbus/2,-Vbus/2],'--k')

% tensiones ph-ph ya filtradas
Vab = Vfase_f(1,:)-Vfase_f(2,:);
Vac = Vfase_f(1,:)-Vfase_f(3,:);
Vfase_fase_f = [Vab; Vac];
figure
hold on
plot(tray_sim_tiempo,Vfase_fase_f)
% plot([0,tray_sim_tiempo(end)],[Vmax_ph_ph,Vmax_ph_ph],'--k')
% plot([0,tray_sim_tiempo(end)],[-Vmax_ph_ph,-Vmax_ph_ph],'--k')

% % tensiones en el VSI
% Vvsi = zeros(length(tray_sim_tiempo),5);
% for i=1:length(tray_sim_tiempo)
%     Vvsi(i,:) = Vfase_f(:,i)' - 0.5*(max(Vfase_f(:,i)) + min(Vfase_f(:,i)));
% end
% figure
% hold on
% plot(tray_sim_tiempo,Vvsi)
% plot([0,tray_sim_tiempo(end)],[Vmax_ph_VSI,Vmax_ph_VSI],'--k')
% plot([0,tray_sim_tiempo(end)],[-Vmax_ph_VSI,-Vmax_ph_VSI],'--k')

% vn = sum(XI32(tray_sim_indopti,:),2)*Vbus/5;
% Vvsi = zeros(length(vn),5);
% for i=1:length(vn)
%     Vvsi(i,:) = Vfase_f(:,i)' - Vbus/2 + vn(i);
% end
% figure
% hold on
% plot(tray_sim_tiempo,Vvsi)
% plot([0,tray_sim_tiempo(end)],[Vmax_ph_VSI,Vmax_ph_VSI],'--k')
% plot([0,tray_sim_tiempo(end)],[-Vmax_ph_VSI,-Vmax_ph_VSI],'--k')

figure
plot(tray_sim_tiempo,tray_sim_par_real(1:end-1))
hold on 
plot(tray_sim_tiempo,tray_sim_par_refe,'r')

figure
plot(tray_sim_tiempo,tray_sim_velocidad(1:end-1))
hold on
plot(tray_sim_tiempo,tray_sim_velocidad_refe,'r')

% figure
% plot(tray_sim_tiempo,tray_sim_bemf_dq)
% 
% figure
% plot(tray_sim_tiempo,tray_sim_bemf_fase)

%% Gráficas de resultados
% Se guardan las variables a representar
tiempo = tray_sim_tiempo;
Ia = tray_sim_ifase_real(:,1);
Ib = tray_sim_ifase_real(:,2);
Ic = tray_sim_ifase_real(:,3);
Ia_refe = tray_sim_ifase_refe(:,1);
Ib_refe = tray_sim_ifase_refe(:,2);
Ic_refe = tray_sim_ifase_refe(:,3);
Id_refe = tray_sim_idq_refe(1:param_sim_Nh,1);
Iq_refe = tray_sim_idq_refe(1:param_sim_Nh,2);
Id = tray_sim_idq_real(1:param_sim_Nh,1);
Iq = tray_sim_idq_real(1:param_sim_Nh,2);
wreal = tray_sim_velocidad(1:param_sim_Nh);
wref  = tray_sim_velocidad_refe(1:param_sim_Nh);
torque_ref = tray_sim_par_refe;
torque = tray_sim_par_real(1:param_sim_Nh);
Va = tray_sim_vfase_real(1:param_sim_Nh,1);
Vb = tray_sim_vfase_real(1:param_sim_Nh,2);
Vc = tray_sim_vfase_real(1:param_sim_Nh,3);
% Va_VSI = tray_sim_vfase_vsi(1:param_sim_Nh,1);
% Vb_VSI = tray_sim_vfase_vsi(1:param_sim_Nh,2);
% Vc_VSI = tray_sim_vfase_vsi(1:param_sim_Nh,3);
% Vd_VSI = tray_sim_vfase_vsi(1:param_sim_Nh,4);
% Ve_VSI = tray_sim_vfase_vsi(1:param_sim_Nh,5);
% Vd1 = tray_sim_vdq_real(1:param_sim_Nh,1);
% Vq1 = tray_sim_vdq_real(1:param_sim_Nh,2);
% Vd3 = tray_sim_vdq_real(1:param_sim_Nh,3);
% Vq3 = tray_sim_vdq_real(1:param_sim_Nh,4);

% Establecemos las características de las figuras
tam_lin=1.5;      % Tamaño de las líneas

% Creamos los títulos de las figuras
p11  = figure('Name','Iphase','NumberTitle','off');
p12  = figure('Name','Iphase','NumberTitle','off');
p13  = figure('Name','Iphase','NumberTitle','off');
p41  = figure('Name','Idq','NumberTitle','off');
p5  = figure('Name','Speed','NumberTitle','off');
p6  = figure('Name','Torque','NumberTitle','off');
p151 = figure('Name','Vphase-to-phase filtered','NumberTitle','off');
p152 = figure('Name','Vphase-to-phase filtered','NumberTitle','off');
p153 = figure('Name','Vphase-to-phase filtered','NumberTitle','off');

% Figura de velocidad
figure(p5)
set(p5,'OuterPosition',[1 210 430 300]);
plot(tiempo,wreal,'Color',[0, 0, 153]./256,'LineWidth',tam_lin)
hold on
plot(tiempo,wref,'Color','r','LineWidth',tam_lin)
vel=legend('\it{\omega}','\it{\omega}^{*}','Location','NorthEast','Orientation','Horizontal');
set(vel,'FontSize',10,'FontName','Times New Roman')
xlabel('Time (s)','FontSize',10)
ylabel('Speed (rad/s)','FontSize',10)
grid on
% xlim([0 tsim-delay])
% ylim([omegam_ref-(1*2*pi/60) omegam_ref+(1*2*pi/60)])
box on

% Figura de par
figure(p6)
set(p6,'OuterPosition',[1 210 430 300]);
hold on
h(2)=plot(tiempo,torque,'Color',[0, 0, 153]./256,'LineWidth',tam_lin);
h(1)=plot(tiempo,torque_ref,'Color','r','LineWidth',tam_lin);
c={'\it{T^*_e_m}','\it{T_e_m}'};
order=[2 1];
tor=legend(h(order),c(order));
set(tor,'FontSize',10,'FontName','Times New Roman','Location','NorthEast','Orientation','Horizontal')
xlabel('Time (s)','FontSize',10)
ylabel('Torque (N·m)','FontSize',10)
grid on
% xlim([0 tsim-delay])
% ylim([-10 10])
box on


% Figura de Idq
figure(p41)
set(p41,'OuterPosition',[1 210 430 300]);
plot(tiempo,Id,'Color',[255, 140, 0]./256,'LineWidth',tam_lin);
hold on
plot(tiempo,Iq,'Color',[0, 128, 0]./256,'LineWidth',tam_lin)
hold on
% plot(tiempo,Id3,'Color',[204, 0, 51]./256,'LineWidth',tam_lin)
% hold on
% plot(tiempo,Iq3,'Color',[0, 0, 153]./256,'LineWidth',tam_lin)
% hold on
plot(tiempo,[Id_refe, Iq_refe],'k--','LineWidth',tam_lin)
hold on
dq=legend('\it{i_{d}}','\it{i_{q}}');
set(dq,'FontSize',10,'FontName','Times New Roman','Location','NorthEast','Orientation','Horizontal')
xlabel('Time (s)','FontSize',10)
ylabel('Current (A)','FontSize',10)
grid on
% xlim([0 tsim-delay])
% ylim([-10 10])
box on


% Figura de corrientes de fase
figure(p11)
set(p11,'OuterPosition',[1 210 430 300]);
plot(tiempo,Ia,'Color',[255, 140, 0]./256,'LineWidth',tam_lin)
hold on
plot(tiempo,Ib,'Color',[0, 128, 0]./256,'LineWidth',tam_lin)
hold on
plot(tiempo,Ic,'Color',[204, 0, 51]./256,'LineWidth',tam_lin)
hold on
plot(tiempo,[Ia_refe, Ib_refe, Ic_refe],'k--','LineWidth',tam_lin)
hold on
% plot([0,tiempo(end)],[Imax_ph,Imax_ph],'--k','LineWidth',tam_lin)
% plot([0,tiempo(end)],[-Imax_ph,-Imax_ph],'--k','LineWidth',tam_lin)
fase=legend('\it{i_{a}}','\it{i_{b}}','\it{i_{c}}','Location','NorthEast','Orientation','Horizontal');
set(fase,'FontSize',10,'FontName','Times New Roman')
xlabel('Time (s)','FontSize',10)
ylabel('Current (A)','FontSize',10)
grid on
xlim([0.44 0.48])
% ylim([-5 5])
box on



% figure(p12)
% set(p12,'OuterPosition',[1 210 430 300]);
% plot(tiempo,Ia,'Color',[255, 140, 0]./256,'LineWidth',tam_lin)
% hold on
% plot(tiempo,Ib,'Color',[0, 128, 0]./256,'LineWidth',tam_lin)
% hold on
% plot(tiempo,Ic,'Color',[204, 0, 51]./256,'LineWidth',tam_lin)
% hold on
% plot(tiempo,[Ia_refe, Ib_refe, Ic_refe],'k--','LineWidth',tam_lin)
% hold on
% % plot([0,tiempo(end)],[Imax_ph,Imax_ph],'--k','LineWidth',tam_lin)
% % plot([0,tiempo(end)],[-Imax_ph,-Imax_ph],'--k','LineWidth',tam_lin)
% fase=legend('\it{i_{a}}','\it{i_{b}}','\it{i_{c}}','Location','NorthEast','Orientation','Horizontal');
% set(fase,'FontSize',10,'FontName','Times New Roman')
% xlabel('Time (s)','FontSize',10)
% ylabel('Current (A)','FontSize',10)
% grid on
% % xlim([0.92 0.94])
% % ylim([-75 125])
% box on


% figure(p13)
% set(p13,'OuterPosition',[1 210 430 300]);
% plot(tiempo,Ia,'Color',[255, 140, 0]./256,'LineWidth',tam_lin)
% hold on
% plot(tiempo,Ib,'Color',[0, 128, 0]./256,'LineWidth',tam_lin)
% hold on
% plot(tiempo,Ic,'Color',[204, 0, 51]./256,'LineWidth',tam_lin)
% hold on
% plot(tiempo,[Ia_refe, Ib_refe, Ic_refe],'k--','LineWidth',tam_lin)
% hold on
% % plot([0,tiempo(end)],[Imax_ph,Imax_ph],'--k','LineWidth',tam_lin)
% % plot([0,tiempo(end)],[-Imax_ph,-Imax_ph],'--k','LineWidth',tam_lin)
% fase=legend('\it{i_{a}}','\it{i_{b}}','\it{i_{c}}','Location','NorthEast','Orientation','Horizontal');
% set(fase,'FontSize',10,'FontName','Times New Roman')
% xlabel('Time (s)','FontSize',10)
% ylabel('Current (A)','FontSize',10)
% grid on
% % xlim([1.32 1.33])
% % ylim([-75 125])
% box on


% Figura de tensiones de fase a fase filtradas
figure(p151)
set(p151,'OuterPosition',[1 210 430 300]);
plot(tiempo,Vab,'Color',[255, 140, 0]./256,'LineWidth',tam_lin)
hold on
plot(tiempo,Vac,'Color',[0, 128, 0]./256,'LineWidth',tam_lin)
hold on
plot([0,tiempo(end)],[Vbus,Vbus],'--k','LineWidth',tam_lin)
plot([0,tiempo(end)],[-Vbus,-Vbus],'--k','LineWidth',tam_lin)
fase=legend('\it{v_{ab}}','\it{v_{ac}}','Location','NorthEast','Orientation','Horizontal');
set(fase,'FontSize',10,'FontName','Times New Roman')
xlabel('Time (s)','FontSize',10)
ylabel('Voltage (V)','FontSize',10)
grid on
xlim([0.42 0.48])
% ylim([-50 100])
box on

 
% figure(p152)
% set(p152,'OuterPosition',[1 210 430 300]);
% plot(tiempo,Vab,'Color',[255, 140, 0]./256,'LineWidth',tam_lin)
% hold on
% plot(tiempo,Vac,'Color',[0, 128, 0]./256,'LineWidth',tam_lin)
% hold on
% plot([0,tiempo(end)],[Vbus,Vbus],'--k','LineWidth',tam_lin)
% plot([0,tiempo(end)],[-Vbus,-Vbus],'--k','LineWidth',tam_lin)
% fase=legend('\it{v_{ab}}','\it{v_{ac}}','Location','NorthEast','Orientation','Horizontal');
% set(fase,'FontSize',10,'FontName','Times New Roman')
% xlabel('Time (s)','FontSize',10)
% ylabel('Voltage (V)','FontSize',10)
% grid on
% % xlim([0.92 0.94])
% % ylim([-50 100])
% box on
% 
% 
% figure(p153)
% set(p153,'OuterPosition',[1 210 430 300]);
% plot(tiempo,Vab,'Color',[255, 140, 0]./256,'LineWidth',tam_lin)
% hold on
% plot(tiempo,Vac,'Color',[0, 128, 0]./256,'LineWidth',tam_lin)
% hold on
% % plot([0,tiempo(end)],[Vbus,Vbus],'--k','LineWidth',tam_lin)
% % plot([0,tiempo(end)],[-Vbus,-Vbus],'--k','LineWidth',tam_lin)
% fase=legend('\it{v_{ab}}','\it{v_{ac}}','Location','NorthEast','Orientation','Horizontal');
% set(fase,'FontSize',10,'FontName','Times New Roman')
% xlabel('Time (s)','FontSize',10)
% ylabel('Voltage (V)','FontSize',10)
% grid on
% % xlim([1.32 1.33])
% % ylim([-50 100])
% box on


