################################################################
# This makes the arrays w1, w2 whcih then are transformed to solve
# for fluid velocity
#################################################################

def makew(w1,w2,f1,f2,u,v):
    import numpy as np
    import fluid_constants as f 

    #calculate constants 
    d_dx = 1.0/f.dx
    dt_p = f.dt/f.p
    #print u, v
# Comput v values from the peskin paper

    for i in range(f.M-1):
        for j in range(f.N-1):    #may want to take away the -1 
            # for periodic boundary we need to have special values for our indexing. 
            ip1 = i+1
            jp1 = j+1
            i_1 = i-1
            j_1 = j-1

            if i==f.M:
                ip1 = 0

            if j==f.N:
                jp1 = 0

            if i==0:
                i_1 = f.M-1

            if j==0:
                j_1 = f.N-1

        # we are now ready to compute the v's from peskins paper using u,v,f1,f2
        # First thing we need to determine is whether we use upwind or downwind differencing 

            if u[i,j] < 0:
                du1 = (1/f.dx)*(u[ip1,j] - u[i,j])
                dv1 = (1/f.dx)*(v[ip1,j] - v[i,j])
            else:
                du1 = (1/f.dx)*(u[i,j] - u[i_1,j])
                dv1 = (1/f.dx)*(v[i,j] - v[i_1,j])

            if v[i,j] < 0: 
                du2 = (1/f.dx)*(u[i,jp1] - u[i,j])
                dv2 = (1/f.dx)*(v[i,jp1] - v[i,j])
            else:
                du2 = (1/f.dx)*(u[i,j] - u[i,j_1])
                dv2 = (1/f.dx)*(v[i,j] - v[i,j_1])


# We are now ready to calculate the values of w1, w2 eqn 14.42 

            w1[i,j] = u[i,j] - f.dt*(u[i,j]*du1 + v[i,j]*du2) + dt_p*f1[i,j]
            w2[i,j] = v[i,j] - f.dt*(u[i,j]*dv1 + v[i,j]*dv2) + dt_p*f2[i,j]

    return w1, w2
            
