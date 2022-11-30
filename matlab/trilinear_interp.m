function value = trilinear_interp(x,y,z,xVec,yVec,zVec,LUT) 

% Buscamos los vértices del cubo
x_aux = find(xVec <= x);
index_x0 = x_aux(end);
index_x1 = x_aux(end)+1;
x_0 = xVec(index_x0);
x_1 = xVec(index_x1);

y_aux = find(yVec <= y);
index_y0 = y_aux(end);
index_y1 = y_aux(end)+1;
y_0 = yVec(index_y0);
y_1 = yVec(index_y1);

z_aux = find(zVec <= z);
index_z0 = z_aux(end);
index_z1 = z_aux(end)+1;
z_0 = zVec(index_z0);
z_1 = zVec(index_z1);

C_000 = LUT(index_x0, index_y0, index_z0);
C_001 = LUT(index_x0, index_y0, index_z1);
C_010 = LUT(index_x0, index_y1, index_z0);
C_011 = LUT(index_x0, index_y1, index_z1);
C_100 = LUT(index_x1, index_y0, index_z0);
C_101 = LUT(index_x1, index_y0, index_z1);
C_110 = LUT(index_x1, index_y1, index_z0);
C_111 = LUT(index_x1, index_y1, index_z1);

a_0 = (-C_000*x_1*y_1*z_1 + C_001*x_1*y_1*z_0 + C_010*x_1*y_0*z_1 - C_011*x_1*y_0*z_0 + C_100*x_0*y_1*z_1 - C_101*x_0*y_1*z_0 - C_110*x_0*y_0*z_1 + C_111*x_0*y_0*z_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1);
a_1 = (C_000*y_1*z_1 - C_001*y_1*z_0 - C_010*y_0*z_1 + C_011*y_0*z_0 - C_100*y_1*z_1 + C_101*y_1*z_0 + C_110*y_0*z_1 - C_111*y_0*z_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1);
a_2 = (C_000*x_1*z_1 - C_001*x_1*z_0 - C_010*x_1*z_1 + C_011*x_1*z_0 - C_100*x_0*z_1 + C_101*x_0*z_0 + C_110*x_0*z_1 - C_111*x_0*z_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1);
a_3 = (C_000*x_1*y_1 - C_001*x_1*y_1 - C_010*x_1*y_0 + C_011*x_1*y_0 - C_100*x_0*y_1 + C_101*x_0*y_1 + C_110*x_0*y_0 - C_111*x_0*y_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1);
a_4 = (-C_000*z_1 + C_001*z_0 + C_010*z_1 - C_011*z_0 + C_100*z_1 - C_101*z_0 - C_110*z_1 + C_111*z_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1);
a_5 = (-C_000*y_1 + C_001*y_1 + C_010*y_0 - C_011*y_0 + C_100*y_1 - C_101*y_1 - C_110*y_0 + C_111*y_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1);
a_6 = (-C_000*x_1 + C_001*x_1 + C_010*x_1 - C_011*x_1 + C_100*x_0 - C_101*x_0 - C_110*x_0 + C_111*x_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1);
a_7 = (C_000 - C_001 - C_010 + C_011 - C_100 + C_101 + C_110 - C_111)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1);
    
value = (a_0 + a_1*x + a_2*y + a_3*z + a_4*x*y + a_5*x*z + a_6*y*z + a_7*x*y*z);