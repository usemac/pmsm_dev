syms eq_i_a l_aa di_a l_ab di_b l_ac di_c v_a R_s i_a phi_ar omega_r  ... 
eq_i_b l_ba di_a l_bb di_b l_bc di_c v_b R_s i_b phi_br omega_r ...
eq_i_c l_ca di_a l_cb di_b l_cc di_c v_c R_s i_c phi_cr omega_r real

eq_i_a =  l_aa*di_a + l_ab*di_b + l_ac*di_c - v_a + R_s * i_a + phi_ar * omega_r;    
eq_i_b =  l_ba*di_a + l_bb*di_b + l_bc*di_c - v_b + R_s * i_b + phi_br * omega_r;
eq_i_c =  l_ca*di_a + l_cb*di_b + l_cc*di_c - v_c + R_s * i_c + phi_cr * omega_r;

%% Solve
[eq_dia, eq_dib, eq_dic] = solve([eq_i_a, eq_i_b, eq_i_c], di_a, di_b, di_c)

%% Simplificación

eq_dia_symp = simplify(eq_dia)
eq_dib_symp = simplify(eq_dib)
eq_dic_symp = simplify(eq_dic)

%% Derivación de las ecuaciones en dq con flujos
syms theta

a11 = simplify(2/3*(-cos(theta)*sin(theta)-cos(theta-2*pi/3)*sin(theta-2*pi/3)-cos(theta+2*pi/3)*sin(theta+2*pi/3)))
a12 = simplify(2/3*(cos(theta)*cos(theta)+cos(theta-2*pi/3)*cos(theta-2*pi/3)+cos(theta+2*pi/3)*cos(theta+2*pi/3)))
a21 = simplify(2/3*(-sin(theta)*sin(theta)-sin(theta-2*pi/3)*sin(theta-2*pi/3)-sin(theta+2*pi/3)*sin(theta+2*pi/3)))
a22 = simplify(2/3*(sin(theta)*cos(theta)+sin(theta-2*pi/3)*cos(theta-2*pi/3)+sin(theta+2*pi/3)*cos(theta+2*pi/3)))
a31 = simplify(2/3*(-sin(theta)*1/2-sin(theta-2*pi/3)*1/2-sin(theta+2*pi/3)*1/2))
a32 = simplify(2/3*(1/2*cos(theta)+1/2*cos(theta-2*pi/3)+1/2*cos(theta+2*pi/3)))
