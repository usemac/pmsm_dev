import numpy as np
import numba


def val2idx_i2p(i_d,i_q,theta,I_d_i2p,I_q_i2p,Theta_i2p):
    '''
    Returns indices in vectors I_d, I_q and Theta for given i_d, i_q and theta values.
    '''
    
    i_d_min,i_d_max,N_i_d = I_d_i2p[0],I_d_i2p[-1],len(I_d_i2p)
    i_q_min,i_q_max,N_i_q = I_q_i2p[0],I_q_i2p[-1],len(I_q_i2p)
    theta_min,theta_max,N_th = Theta_i2p[0],Theta_i2p[-1],len(Theta_i2p)
   
    i_d_idx = int((i_d - i_d_min)/(i_d_max - i_d_min)*(N_i_d))
    i_q_idx = int((i_q - i_q_min)/(i_q_max - i_q_min)*(N_i_q))
    theta_idx = int((theta - theta_min)/(theta_max - theta_min)*(N_th))
    
    return i_d_idx,i_q_idx,theta_idx

def phi_eval(theta,i_d,i_q):
    
    i_d_idx = int((i_d - I_d_min)/(I_d_max - I_d_min)*(N_i_d))
    i_q_idx = int((i_q - I_q_min)/(I_q_max - I_q_min)*(N_i_q))
    theta_idx = int(((theta) - theta_min)/(theta_max - theta_min)*(N_th))

    phi_d = data['phi_d'][i_d_idx,i_q_idx,theta_idx]
    phi_q = data['phi_q'][i_d_idx,i_q_idx,theta_idx]
    
    #dphi_d_i_d = -(phi_d - data['phi_d'][i_d_idx+1,  i_q_idx,theta_idx])/30
    #dphi_q_i_q = -(phi_q - data['phi_q'][  i_d_idx,i_q_idx+1,theta_idx])/30
    
    L_dd = (data['phi_d'][i_d_idx+1,  i_q_idx,0]-data['phi_d'][i_d_idx-1,  i_q_idx,0])/(data['i_d'][i_d_idx+1]-data['i_d'][i_d_idx-1])
    L_qq = (data['phi_q'][  i_d_idx,i_q_idx+1,0]-data['phi_q'][  i_d_idx,i_q_idx-1,0])/(data['i_q'][i_q_idx+1]-data['i_q'][i_q_idx-1])
    L_dq = (data['phi_d'][  i_d_idx,i_q_idx+1,0]-data['phi_d'][  i_d_idx,i_q_idx-1,0])/(data['i_q'][i_q_idx+1]-data['i_q'][i_q_idx-1])
    L_qd = (data['phi_q'][i_d_idx+1,  i_q_idx,0]-data['phi_q'][i_d_idx-1,  i_q_idx,0])/(data['i_d'][i_d_idx+1]-data['i_d'][i_d_idx-1])

    L_s = np.array([[L_dd,L_dq],
                    [L_qd,L_qq]])
    
    return phi_d,phi_q,L_s

@numba.njit()
def phi_eval_2(i_d,i_q,theta,Phi_d,Phi_q,I_d,I_q,Theta):
    
    I_d_min = I_d[0]
    I_d_max = I_d[-1]
    N_i_d = len(I_d)
    i_d_idx_0 = np.int32(np.floor((i_d - I_d_min)/(I_d_max - I_d_min)*(N_i_d-1)))
    i_d_idx_1 = i_d_idx_0
    x_0 = I_d[i_d_idx_0]
    x_1 = x_0 + 0.01
    if i_d_idx_0 < (N_i_d-1):
        i_d_idx_1 = i_d_idx_0 + 1
        x_1 = I_d[i_d_idx_1] 
        
    I_q_min = I_q[0]
    I_q_max = I_q[-1]
    N_i_q = len(I_q)
    i_q_idx_0 = np.int32(np.floor((i_q - I_q_min)/(I_q_max - I_q_min)*(N_i_q-1)))
    i_q_idx_1 = i_q_idx_0
    y_0 = I_q[i_q_idx_0]
    y_1 = y_0 + 0.01
    if i_q_idx_0 < (N_i_q-1):
        i_q_idx_1 = i_q_idx_0 + 1
        y_1 = I_q[i_q_idx_1] 

    Theta_min = Theta[0]
    Theta_max = Theta[-1]
    N_th = len(Theta)
    #if theta==Theta_max: theta=0.99*theta # avoids array overflow
    theta_idx_0 = np.int32(np.floor(((theta) - Theta_min)/(Theta_max - Theta_min)*(N_th-1)))
    theta_idx_1 = theta_idx_0
    z_0 = Theta[theta_idx_0]
    z_1 = 0.01
    if theta_idx_0 < (N_th-1):
        theta_idx_1 = theta_idx_0 + 1
        z_1 = Theta[theta_idx_1]
 
   
    x = i_d
    y = i_q
    z = theta
    
    M = np.array([[1, x_0, y_0, z_0, x_0*y_0, x_0*z_0, y_0*z_0, x_0*y_0*z_0],    
                  [1, x_1, y_0, z_0, x_1*y_0, x_1*z_0, y_0*z_0, x_1*y_0*z_0],    
                  [1, x_0, y_1, z_0, x_0*y_1, x_0*z_0, y_1*z_0, x_0*y_1*z_0],    
                  [1, x_1, y_1, z_0, x_1*y_1, x_1*z_0, y_1*z_0, x_1*y_1*z_0],    
                  [1, x_0, y_0, z_1, x_0*y_0, x_0*z_1, y_0*z_1, x_0*y_0*z_1],    
                  [1, x_1, y_0, z_1, x_1*y_0, x_1*z_1, y_0*z_1, x_1*y_0*z_1],    
                  [1, x_0, y_1, z_1, x_0*y_1, x_0*z_1, y_1*z_1, x_0*y_1*z_1],    
                  [1, x_1, y_1, z_1, x_1*y_1, x_1*z_1, y_1*z_1, x_1*y_1*z_1]])
    
    C_000 = Phi_d[i_d_idx_0,i_q_idx_0,theta_idx_0]
    C_100 = Phi_d[i_d_idx_1,i_q_idx_0,theta_idx_0]   
    C_010 = Phi_d[i_d_idx_0,i_q_idx_1,theta_idx_0]  
    C_110 = Phi_d[i_d_idx_1,i_q_idx_1,theta_idx_0]  
    C_001 = Phi_d[i_d_idx_0,i_q_idx_0,theta_idx_1]
    C_101 = Phi_d[i_d_idx_1,i_q_idx_0,theta_idx_1]   
    C_011 = Phi_d[i_d_idx_0,i_q_idx_1,theta_idx_1]  
    C_111 = Phi_d[i_d_idx_1,i_q_idx_1,theta_idx_1]   
    
    C_d = np.array([[C_000],
                    [C_100],
                    [C_010],
                    [C_110],
                    [C_001],
                    [C_101],
                    [C_011],
                    [C_111]])
    
    a = np.linalg.solve(M,C_d)
    
    phi_d = (a[0] + a[1]*x + a[2]*y + a[3]*z + a[4]*x*y + a[5]*x*z + a[6]*y*z + a[7]*x*y*z)[0]

    C_000 = Phi_q[i_d_idx_0,i_q_idx_0,theta_idx_0]
    C_100 = Phi_q[i_d_idx_1,i_q_idx_0,theta_idx_0]   
    C_010 = Phi_q[i_d_idx_0,i_q_idx_1,theta_idx_0]  
    C_110 = Phi_q[i_d_idx_1,i_q_idx_1,theta_idx_0]  
    C_001 = Phi_q[i_d_idx_0,i_q_idx_0,theta_idx_1]
    C_101 = Phi_q[i_d_idx_1,i_q_idx_0,theta_idx_1]   
    C_011 = Phi_q[i_d_idx_0,i_q_idx_1,theta_idx_1]  
    C_111 = Phi_q[i_d_idx_1,i_q_idx_1,theta_idx_1]   
    
    C_q = np.array([[C_000],
                    [C_100],
                    [C_010],
                    [C_110],
                    [C_001],
                    [C_101],
                    [C_011],
                    [C_111]])
    
    a = np.linalg.solve(M,C_q)

    phi_q = (a[0] + a[1]*x + a[2]*y + a[3]*z + a[4]*x*y + a[5]*x*z + a[6]*y*z + a[7]*x*y*z)[0]
    
    if i_d_idx_0<(N_i_d-1) and i_q_idx_0<(N_i_q-1):
        L_dd = (Phi_d[i_d_idx_0+1,   i_q_idx_0,0] - Phi_d[i_d_idx_0-1, i_q_idx_0,0])/(I_d[i_d_idx_0+1]-I_d[i_d_idx_0-1])
        L_qq = (Phi_q[  i_d_idx_0, i_q_idx_0+1,0] - Phi_q[i_d_idx_0, i_q_idx_0-1,0])/(I_q[i_q_idx_0+1]-I_q[i_q_idx_0-1])
        L_dq = (Phi_d[  i_d_idx_0, i_q_idx_0+1,0] - Phi_d[i_d_idx_0, i_q_idx_0-1,0])/(I_q[i_q_idx_0+1]-I_q[i_q_idx_0-1])
        L_qd = (Phi_q[i_d_idx_0+1, i_q_idx_0,0] - Phi_q[i_d_idx_0-1, i_q_idx_0,0])/(I_d[i_d_idx_0+1]-I_d[i_d_idx_0-1])
    else:
        L_dd = 1.0
        L_qq = 1.0
        L_dq = 1.0
        L_qd = 1.0      

    L_s = np.array([[L_dd,L_dq],
                    [L_qd,L_qq]])
    
    return phi_d,phi_q,L_s


@numba.njit()
def i_eval_2(phi_d,phi_q,theta,Phi_d,Phi_q,I_d,I_q,Theta):
    
    Phi_d_min = Phi_d[0]
    Phi_d_max = Phi_d[-1]
    N_phi_d = len(Phi_d)
    phi_d_idx_0 = np.int32(np.floor((phi_d - Phi_d_min)/(Phi_d_max - Phi_d_min)*(N_phi_d-1)))
    phi_d_idx_1 = phi_d_idx_0 + 1

    Phi_q_min = Phi_q[0]
    Phi_q_max = Phi_q[-1]
    N_phi_q = len(Phi_q)
    phi_q_idx_0 = np.int32(np.floor((phi_q - Phi_q_min)/(Phi_q_max - Phi_q_min)*(N_phi_q-1)))
    phi_q_idx_1 = phi_q_idx_0 + 1
 
    Theta_min = Theta[0]
    Theta_max = Theta[-1]
    N_th = len(Theta)
    theta_idx_0 = np.int32(np.floor(((theta) - Theta_min)/(Theta_max - Theta_min)*(N_th-1)))
    theta_idx_1 = theta_idx_0
    z_0 = Theta[theta_idx_0]
    z_1 = z_0 + 0.0001
    if theta_idx_0 < (N_th-1):
        theta_idx_1 = theta_idx_0 + 1
        z_1 = Theta[theta_idx_1]
        
 
    x_0 = Phi_d[phi_d_idx_0]
    x_1 = Phi_d[phi_d_idx_1]
    y_0 = Phi_q[phi_q_idx_0]
    y_1 = Phi_q[phi_q_idx_1]

    
    x = phi_d
    y = phi_q
    z = theta
    
    M = np.array([[1, x_0, y_0, z_0, x_0*y_0, x_0*z_0, y_0*z_0, x_0*y_0*z_0],    
                  [1, x_1, y_0, z_0, x_1*y_0, x_1*z_0, y_0*z_0, x_1*y_0*z_0],    
                  [1, x_0, y_1, z_0, x_0*y_1, x_0*z_0, y_1*z_0, x_0*y_1*z_0],    
                  [1, x_1, y_1, z_0, x_1*y_1, x_1*z_0, y_1*z_0, x_1*y_1*z_0],    
                  [1, x_0, y_0, z_1, x_0*y_0, x_0*z_1, y_0*z_1, x_0*y_0*z_1],    
                  [1, x_1, y_0, z_1, x_1*y_0, x_1*z_1, y_0*z_1, x_1*y_0*z_1],    
                  [1, x_0, y_1, z_1, x_0*y_1, x_0*z_1, y_1*z_1, x_0*y_1*z_1],    
                  [1, x_1, y_1, z_1, x_1*y_1, x_1*z_1, y_1*z_1, x_1*y_1*z_1]])
    
    C_000 = I_d[phi_d_idx_0,phi_q_idx_0,theta_idx_0]
    C_100 = I_d[phi_d_idx_1,phi_q_idx_0,theta_idx_0]   
    C_010 = I_d[phi_d_idx_0,phi_q_idx_1,theta_idx_0]  
    C_110 = I_d[phi_d_idx_1,phi_q_idx_1,theta_idx_0]  
    C_001 = I_d[phi_d_idx_0,phi_q_idx_0,theta_idx_1]
    C_101 = I_d[phi_d_idx_1,phi_q_idx_0,theta_idx_1]   
    C_011 = I_d[phi_d_idx_0,phi_q_idx_1,theta_idx_1]  
    C_111 = I_d[phi_d_idx_1,phi_q_idx_1,theta_idx_1]   
    
    C_d = np.array([[C_000],
                    [C_100],
                    [C_010],
                    [C_110],
                    [C_001],
                    [C_101],
                    [C_011],
                    [C_111]])
    
    a = np.linalg.solve(M,C_d)
    
    i_d = (a[0] + a[1]*x + a[2]*y + a[3]*z + a[4]*x*y + a[5]*x*z + a[6]*y*z + a[7]*x*y*z)[0]

    C_000 = I_q[phi_d_idx_0,phi_q_idx_0,theta_idx_0]
    C_100 = I_q[phi_d_idx_1,phi_q_idx_0,theta_idx_0]   
    C_010 = I_q[phi_d_idx_0,phi_q_idx_1,theta_idx_0]  
    C_110 = I_q[phi_d_idx_1,phi_q_idx_1,theta_idx_0]  
    C_001 = I_q[phi_d_idx_0,phi_q_idx_0,theta_idx_1]
    C_101 = I_q[phi_d_idx_1,phi_q_idx_0,theta_idx_1]   
    C_011 = I_q[phi_d_idx_0,phi_q_idx_1,theta_idx_1]  
    C_111 = I_q[phi_d_idx_1,phi_q_idx_1,theta_idx_1]   
    
    C_q = np.array([[C_000],
                    [C_100],
                    [C_010],
                    [C_110],
                    [C_001],
                    [C_101],
                    [C_011],
                    [C_111]])
 
    a = np.linalg.solve(M,C_q)

    i_q = (a[0] + a[1]*x + a[2]*y + a[3]*z + a[4]*x*y + a[5]*x*z + a[6]*y*z + a[7]*x*y*z)[0]

    return i_d,i_q

@numba.njit()
def trilinear(x,y,z,x_0,y_0,z_0,x_1,y_1,z_1,C_000,C_100,C_010,C_110,C_001,C_101,C_011,C_111):

    a_0 = (-C_000*x_1*y_1*z_1 + C_001*x_1*y_1*z_0 + C_010*x_1*y_0*z_1 - C_011*x_1*y_0*z_0 + C_100*x_0*y_1*z_1 - C_101*x_0*y_1*z_0 - C_110*x_0*y_0*z_1 + C_111*x_0*y_0*z_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1)
    a_1 = (C_000*y_1*z_1 - C_001*y_1*z_0 - C_010*y_0*z_1 + C_011*y_0*z_0 - C_100*y_1*z_1 + C_101*y_1*z_0 + C_110*y_0*z_1 - C_111*y_0*z_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1)
    a_2 = (C_000*x_1*z_1 - C_001*x_1*z_0 - C_010*x_1*z_1 + C_011*x_1*z_0 - C_100*x_0*z_1 + C_101*x_0*z_0 + C_110*x_0*z_1 - C_111*x_0*z_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1)
    a_3 = (C_000*x_1*y_1 - C_001*x_1*y_1 - C_010*x_1*y_0 + C_011*x_1*y_0 - C_100*x_0*y_1 + C_101*x_0*y_1 + C_110*x_0*y_0 - C_111*x_0*y_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1)
    a_4 = (-C_000*z_1 + C_001*z_0 + C_010*z_1 - C_011*z_0 + C_100*z_1 - C_101*z_0 - C_110*z_1 + C_111*z_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1)
    a_5 = (-C_000*y_1 + C_001*y_1 + C_010*y_0 - C_011*y_0 + C_100*y_1 - C_101*y_1 - C_110*y_0 + C_111*y_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1)
    a_6 = (-C_000*x_1 + C_001*x_1 + C_010*x_1 - C_011*x_1 + C_100*x_0 - C_101*x_0 - C_110*x_0 + C_111*x_0)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1)
    a_7 = (C_000 - C_001 - C_010 + C_011 - C_100 + C_101 + C_110 - C_111)/(x_0*y_0*z_0 - x_0*y_0*z_1 - x_0*y_1*z_0 + x_0*y_1*z_1 - x_1*y_0*z_0 + x_1*y_0*z_1 + x_1*y_1*z_0 - x_1*y_1*z_1)
    
    value = (a_0 + a_1*x + a_2*y + a_3*z + a_4*x*y + a_5*x*z + a_6*y*z + a_7*x*y*z)
    
    return value 


@numba.njit()
def phi_eval_3(i_d,i_q,theta,Phi_d,Phi_q,Torque,I_d,I_q,Theta):
    
    I_d_min = I_d[0]
    I_d_max = I_d[-1]
    N_i_d = len(I_d)
    i_d_idx_0 = np.int32(np.floor((i_d - I_d_min)/(I_d_max - I_d_min)*(N_i_d-1)))
    i_d_idx_1 = i_d_idx_0
    x_0 = I_d[i_d_idx_0]
    x_1 = x_0 + 0.01
    if i_d_idx_0 < (N_i_d-1):
        i_d_idx_1 = i_d_idx_0 + 1
        x_1 = I_d[i_d_idx_1] 
        
    I_q_min = I_q[0]
    I_q_max = I_q[-1]
    N_i_q = len(I_q)
    i_q_idx_0 = np.int32(np.floor((i_q - I_q_min)/(I_q_max - I_q_min)*(N_i_q-1)))
    i_q_idx_1 = i_q_idx_0
    y_0 = I_q[i_q_idx_0]
    y_1 = y_0 + 0.01
    if i_q_idx_0 < (N_i_q-1):
        i_q_idx_1 = i_q_idx_0 + 1
        y_1 = I_q[i_q_idx_1] 

    Theta_min = Theta[0]
    Theta_max = Theta[-1]
    N_th = len(Theta)
    #if theta==Theta_max: theta=0.99*theta # avoids array overflow
    theta_idx_0 = np.int32(np.floor(((theta) - Theta_min)/(Theta_max - Theta_min)*(N_th-1)))
    theta_idx_1 = theta_idx_0 + 1
    if theta_idx_1 >= N_th:
        theta_idx_1 = 0       
    z_0 = Theta[theta_idx_0]
    z_1 = Theta[theta_idx_1]
 
    x = i_d
    y = i_q
    z = theta

    # phi_d
    C_000 = Phi_d[i_d_idx_0,i_q_idx_0,theta_idx_0]
    C_100 = Phi_d[i_d_idx_1,i_q_idx_0,theta_idx_0]   
    C_010 = Phi_d[i_d_idx_0,i_q_idx_1,theta_idx_0]  
    C_110 = Phi_d[i_d_idx_1,i_q_idx_1,theta_idx_0]  
    C_001 = Phi_d[i_d_idx_0,i_q_idx_0,theta_idx_1]
    C_101 = Phi_d[i_d_idx_1,i_q_idx_0,theta_idx_1]   
    C_011 = Phi_d[i_d_idx_0,i_q_idx_1,theta_idx_1]  
    C_111 = Phi_d[i_d_idx_1,i_q_idx_1,theta_idx_1]   
       
    phi_d = trilinear(x,y,z,x_0,y_0,z_0,x_1,y_1,z_1,C_000,C_100,C_010,C_110,C_001,C_101,C_011,C_111)

    # phi_q
    C_000 = Phi_q[i_d_idx_0,i_q_idx_0,theta_idx_0]
    C_100 = Phi_q[i_d_idx_1,i_q_idx_0,theta_idx_0]   
    C_010 = Phi_q[i_d_idx_0,i_q_idx_1,theta_idx_0]  
    C_110 = Phi_q[i_d_idx_1,i_q_idx_1,theta_idx_0]  
    C_001 = Phi_q[i_d_idx_0,i_q_idx_0,theta_idx_1]
    C_101 = Phi_q[i_d_idx_1,i_q_idx_0,theta_idx_1]   
    C_011 = Phi_q[i_d_idx_0,i_q_idx_1,theta_idx_1]  
    C_111 = Phi_q[i_d_idx_1,i_q_idx_1,theta_idx_1]   
       
    phi_q = trilinear(x,y,z,x_0,y_0,z_0,x_1,y_1,z_1,C_000,C_100,C_010,C_110,C_001,C_101,C_011,C_111)
    
    # torque
    C_000 = Torque[i_d_idx_0,i_q_idx_0,theta_idx_0]
    C_100 = Torque[i_d_idx_1,i_q_idx_0,theta_idx_0]   
    C_010 = Torque[i_d_idx_0,i_q_idx_1,theta_idx_0]  
    C_110 = Torque[i_d_idx_1,i_q_idx_1,theta_idx_0]  
    C_001 = Torque[i_d_idx_0,i_q_idx_0,theta_idx_1]
    C_101 = Torque[i_d_idx_1,i_q_idx_0,theta_idx_1]   
    C_011 = Torque[i_d_idx_0,i_q_idx_1,theta_idx_1]  
    C_111 = Torque[i_d_idx_1,i_q_idx_1,theta_idx_1]   
       
    torque = trilinear(x,y,z,x_0,y_0,z_0,x_1,y_1,z_1,C_000,C_100,C_010,C_110,C_001,C_101,C_011,C_111)

    
    if i_d_idx_0<(N_i_d-1) and i_q_idx_0<(N_i_q-1):
        L_dd = (Phi_d[i_d_idx_0+1,   i_q_idx_0,0] - Phi_d[i_d_idx_0-1, i_q_idx_0,0])/(I_d[i_d_idx_0+1]-I_d[i_d_idx_0-1])
        L_qq = (Phi_q[  i_d_idx_0, i_q_idx_0+1,0] - Phi_q[i_d_idx_0, i_q_idx_0-1,0])/(I_q[i_q_idx_0+1]-I_q[i_q_idx_0-1])
        L_dq = (Phi_d[  i_d_idx_0, i_q_idx_0+1,0] - Phi_d[i_d_idx_0, i_q_idx_0-1,0])/(I_q[i_q_idx_0+1]-I_q[i_q_idx_0-1])
        L_qd = (Phi_q[i_d_idx_0+1, i_q_idx_0,0] - Phi_q[i_d_idx_0-1, i_q_idx_0,0])/(I_d[i_d_idx_0+1]-I_d[i_d_idx_0-1])
    else:
        L_dd = 1.0
        L_qq = 1.0
        L_dq = 1.0
        L_qd = 1.0      

    L_s = np.array([[L_dd,L_dq],
                    [L_qd,L_qq]])
    
    return phi_d,phi_q,torque,L_s


@numba.njit()
def park_1(abc,theta):
    
    dqz = np.zeros(3,dtype=np.float64)
    
    dqz[0] = 2.0/3.0*(abc[0]*np.cos(theta) + abc[1]*np.cos(theta-2/3*np.pi) + abc[2]*np.cos(theta-4/3*np.pi));
    dqz[1] =-2.0/3.0*(abc[0]*np.sin(theta) + abc[1]*np.sin(theta-2/3*np.pi) + abc[2]*np.sin(theta-4/3*np.pi));
    dqz[2] = 1.0/3.0*np.sum(abc);
    
    return dqz


@numba.njit()
def ipark_1(dqz,theta):
    abc = np.zeros(3,dtype=np.float64)
    
    abc[0] = (dqz[0]*np.cos(theta)               - dqz[1]*np.sin(theta));
    abc[1] = (dqz[0]*np.cos(theta-2.0/3.0*np.pi) - dqz[1]*np.sin(theta-2.0/3.0*np.pi));
    abc[2] = (dqz[0]*np.cos(theta-4.0/3.0*np.pi) - dqz[1]*np.sin(theta-4.0/3.0*np.pi));
    
    return abc
    
@numba.njit()
def park_2(abc,theta):
    
    dqz = np.zeros(3,dtype=np.float64)
    
    dqz[0] = 2.0/3.0*(abc[0]*np.cos(theta) + abc[1]*np.cos(theta-2/3*np.pi) + abc[2]*np.cos(theta-4/3*np.pi));
    dqz[1] = 2.0/3.0*(abc[0]*np.sin(theta) + abc[1]*np.sin(theta-2/3*np.pi) + abc[2]*np.sin(theta-4/3*np.pi));
    dqz[2] = 1.0/3.0*np.sum(abc);
    
    return dqz


@numba.njit()
def ipark_2(dqz,theta):
    abc = np.zeros(3,dtype=np.float64)
    
    abc[0] = (dqz[0]*np.cos(theta)               + dqz[1]*np.sin(theta)               +  dqz[2]);
    abc[1] = (dqz[0]*np.cos(theta-2.0/3.0*np.pi) + dqz[1]*np.sin(theta-2.0/3.0*np.pi) +  dqz[2]);
    abc[2] = (dqz[0]*np.cos(theta-4.0/3.0*np.pi) + dqz[1]*np.sin(theta-4.0/3.0*np.pi) +  dqz[2]);
    
    return abc
