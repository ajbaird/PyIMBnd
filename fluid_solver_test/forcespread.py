###################################################################
# This code will spread the force generated at the boundary 
# due to the springs, to neighboring fluid gird points 
# array1/2 are arrays of boundary points (x,y points), fb1/2 are the forces 
# and Q are the total points. f1/2 are updated
###################################################################

def forcespread(array1, array2, fb1, fb2, f1, f2, Q):
    
    import numpy as np
    import matplotlib.pylab as plt 


    import fluid_constants as f 

    # first determine which cartesian grid point we are near

    for i in range(Q):
        x1 = np.ceil(array1[i,0]/f.dx)
        x2 = np.ceil(array2[i,0]/f.dx)
       # print 'x1: '+ str(x1) 
       # print 'x2: '+ str(x2)
        # spred forces
        for k in range(4):
            for h in range(4):
                ii = x1+1-k    #may be missing a +1 here .... 
                jj = x2+1-h

                # Check to see if near boundary if it is shift it to other side of the domain (periodic bndry conditions)

                if ii == -1:
                    ii = f.M-1
                if ii == -2:
                    ii = f.M-2 
                if ii == -3:
                    ii = f.M-3
                if ii == f.M:
                    ii = 0 
                if ii == f.M+1:
                    ii = 1
                if ii == f.M+2:
                    ii = 2 

                # now for the jj iterator

                if jj == -1:
                    jj = f.N-1
                if jj == -2:
                    jj = f.N-2
                if jj == -3:
                    jj = f.N-3 
                if jj == f.N:
                    jj = 0 
                if jj == f.N+1:
                    jj = 1
                if jj == f.N+2:
                    jj = 2

                # spread forces with a delta function 
                


                #print 'ii: '+ str(ii) 
                #print 'jj: '+str(jj)
                #taking out a +1 to the f1, f2 indices....
                # there is something here not working! 

                f1[ii,jj] = f1[ii,jj]+fb1[i,0]*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x1-k)-array1[i,0])))*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x2-h)-array2[i,0])))*f.ds  

                f2[ii,jj] = f2[ii,jj]+fb2[i,0]*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x1-k)-array1[i,0])))*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x2-h)-array2[i,0])))*f.ds


    return f1, f2 
