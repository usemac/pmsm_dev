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