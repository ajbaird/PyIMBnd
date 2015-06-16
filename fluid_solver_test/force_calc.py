###############################################################################################
# Computes the force generated in the fluid when perscribed an 
# initial velocity field. This is just used to check the validity 
# of the fluid solver in the immersed boundary method solver 
# 
#
# Input: f1, f2 (x,y force arrays), mu (dynamic viscosity), dx (mesh width)
# length (domain size), p (fluid density), time (simulation time), M, N (number of grid points)
#
# Output: f1,f2 arrays
#
# Written by: Austin Baird to check Immersed boundary method N-S solver
################################################################################################

def force_calc(f1,f2,mu,dx,length,p,time,M,N, theta_1, theta_2): 
    
    import numpy as np 
    import fluid_constants as f 
    import fluid_const_tester as ft 

    f1 = (p*ft.a1[0]*ft.w_1 + ft.b_1*ft.k1[0])*np.cos(ft.theta_1) + f.mu*(ft.k1[0]*ft.k1[0]+ft.k1[1]*ft.k1[1])*ft.a1[0]*np.sin(ft.theta_1) +(p*ft.a2[0]*ft.w_2 + ft.b_2*ft.k2[0])*np.cos(ft.theta_2) + mu*(ft.k2[0]*ft.k2[0]+ft.k2[1]*ft.k2[1])*ft.a2[0]*np.sin(ft.theta_2) + p*np.sin(ft.theta_1)*np.cos(ft.theta_2)*(ft.a1[0]*ft.k2[0]+ft.a1[1]*ft.k2[1])*ft.a2[0] + p*np.sin(ft.theta_2)*np.cos(ft.theta_1)*(ft.a2[0]*ft.k1[0]+ft.a2[1]*ft.k1[1])*ft.a1[0]

    f2 = (p*ft.a1[1]*ft.w_1 + ft.b_1*ft.k1[1])*np.cos(ft.theta_1) + mu*(ft.k1[0]*ft.k1[0]+ft.k1[1]*ft.k1[1])*ft.a1[1]*np.sin(ft.theta_1) +(p*ft.a2[1]*ft.w_2 + ft.b_2*ft.k2[1])*np.cos(ft.theta_2) + mu*(ft.k2[0]*ft.k2[0]+ft.k2[1]*ft.k2[1])*ft.a2[1]*np.sin(ft.theta_2) + p*np.sin(ft.theta_1)*np.cos(ft.theta_2)*(ft.a1[0]*ft.k2[0]+ft.a1[1]*ft.k2[1])*ft.a2[1] + p*np.sin(ft.theta_2)*np.cos(ft.theta_1)*(ft.a2[0]*ft.k1[0]+ft.a2[1]*ft.k1[1])*ft.a1[1]
#    for i in range(M):
#        for j in range(N):
#            f1[i,j] = (p*ft.a1[0]*ft.w_1 + ft.b_1*ft.k1[0])*np.cos(ft.theta_1[i,j]) + f.mu*(ft.k1[0]*ft.k1[0]+ft.k1[1]*ft.k1[1])*ft.a1[0]*np.sin(ft.theta_1[i,j]) +(p*ft.a2[0]*ft.w_2 + ft.b_2*ft.k2[0])*np.cos(ft.theta_2[i,j]) + mu*(ft.k2[0]*ft.k2[0]+ft.k2[1]*ft.k2[1])*ft.a2[0]*np.sin(ft.theta_2[i,j]) + p*np.sin(ft.theta_1[i,j])*np.cos(ft.theta_2[i,j])*(ft.a1[0]*ft.k2[0]+ft.a1[1]*ft.k2[1])*ft.a2[0] + p*np.sin(ft.theta_2[i,j])*np.cos(ft.theta_1[i,j])*(ft.a2[0]*ft.k1[0]+ft.a2[1]*ft.k1[1])*ft.a1[0]

 #           f2[i,j] = (p*ft.a1[1]*ft.w_1 + ft.b_1*ft.k1[1])*np.cos(ft.theta_1[i,j]) + mu*(ft.k1[0]*ft.k1[0]+ft.k1[1]*ft.k1[1])*ft.a1[1]*np.sin(ft.theta_1[i,j]) +(p*ft.a2[1]*ft.w_2 + ft.b_2*ft.k2[1])*np.cos(ft.theta_2[i,j]) + mu*(ft.k2[0]*ft.k2[0]+ft.k2[1]*ft.k2[1])*ft.a2[1]*np.sin(ft.theta_2[i,j]) + p*np.sin(ft.theta_1[i,j])*np.cos(ft.theta_2[i,j])*(ft.a1[0]*ft.k2[0]+ft.a1[1]*ft.k2[1])*ft.a2[1] + p*np.sin(ft.theta_2[i,j])*np.cos(ft.theta_1[i,j])*(ft.a2[0]*ft.k1[0]+ft.a2[1]*ft.k1[1])*ft.a1[1]



    return f1, f2
