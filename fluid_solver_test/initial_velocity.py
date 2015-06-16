##########################################################################
# this function generates a velocity field with an exact solution to the 
# N-S equations. Input into the fluid solver to determine accuracy of 
# solver. Written by Austin Baird, UNC: Chapel Hill
#########################################################################

def initial_velocity(u,v,M,N,dx,length,time): 

    import numpy as np
    import fluid_constants as f
    import fluid_const_tester as ft 

    for i in range(f.M):
        for j in range(f.N):
            ft.theta_1[i,j] = ft.k1[0]*(f.dx*(i)) + ft.k1[1]*(f.dx*(j)) + ft.w_1*time
            ft.theta_2[i,j] = ft.k2[0]*(f.dx*(i)) + ft.k2[1]*(f.dx*(j)) + ft.w_2*time

    u = ft.a1[0]*np.sin(ft.theta_1) + ft.a2[0]*np.sin(ft.theta_2)    
    v = ft.a1[1]*np.sin(ft.theta_1) + ft.a2[1]*np.sin(ft.theta_2)
    return u, v, ft.theta_1, ft.theta_2 
# for i in range(f.M):
#     for j in range(f.N):
#         u[i,j] = ft.a1[0]*np.sin(ft.theta_1[i,j]) + ft.a2[0]*np.sin(ft.theta_2[i,j])
#         v[i,j] = ft.a1[1]*np.sin(ft.theta_1[i,j]) + ft.a2[1]*np.sin(ft.theta_2[i,j])
  





