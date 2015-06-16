########################################################################
# This function moves the boundary at the local fluid velocity 
# b1, b2 are the x,y boundary values, and Q is the number of 
# immersed boundary points. b1, b2 are updated and returned
#######################################################################

def move_boundary(b1, b2, Q, u, v):
    import numpy as np 
    import fluid_constants as f 

    # vectors to hold velocitys of the boundary 

    uarray = np.zeros((Q,1))
    varray = np.zeros((Q,1))


    # determine where you are in cartesian coordinates 
    #print b1, b2
    for i in range(Q):
        x1 = np.ceil(b1[i,0]/f.dx)
        x2 = np.ceil(b2[i,0]/f.dx)

        # loop through points to get avg. velocity, 4 point delta  

        for k in range(4):
            for h in range(4):
                ii = x1-k
                jj = x2-h
               # print ii, jj, x1, k, h, x2
                # check to see if we are within 3 grid points near the boundary

                if ii == -1: 
                    ii = f.M-1 
                if ii == -2: 
                    ii = f.M-2
                if ii == -3:
                    ii = f.M-3 
                if ii == f.M:
                    ii = 0 
                if ii == f.M + 1: 
                    ii = 1
                if  ii == f.M +2: 
                    ii = 2 

                if jj == -1: 
                    jj = f.N-1 
                if jj == -2: 
                    jj = f.N-2
                if jj == -3:
                    jj = f.N-3 
                if jj == f.N:
                    jj = 0 
                if jj == f.N + 1: 
                    jj = 1
                if  jj == f.N +2: 
                    jj = 2 
                #print jj, ii, i
                #print u[ii,jj], b1[i,0]
# calculate the velocity at the boundary by taking a weighted average of the nearby velocities

                uarray[i,0] = uarray[i,0]+u[ii,jj]*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x1-k) - b1[i,0])))*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x2-h)-b2[i,0])))*f.dx**2 
                varray[i,0] = varray[i,0]+v[ii,jj]*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x1-k) - b1[i,0])))*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x2-h)-b2[i,0])))*f.dx**2 
        #print uarray, varray
    
    # now we've computed the local velocity so we can now move the boundary at that speed 

        b1[i,0] = b1[i,0] + uarray[i,0]*f.dt 
        b2[i,0] = b2[i,0] + varray[i,0]*f.dt

    # check to see if we've moved the boundary point out of the domain, if we have move it to the other side

        for i in range(Q):
            if b1[i,0] > f.width: 
                b1[i,0] = b1[i,0] - f.width

            if b2[i,0] > f.width: 
                b2[i,0] = b2[i,0] - f.width 

    return b1, b2
