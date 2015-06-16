#########################################################################################################
#
# This is a two dimensional immersed boundary method code. The fluid solver uses an fft, so the domain 
# must be periodic. New fluid solvers will be added as it gets updated. Right now the boundary is moved 
# via target points (springs and beams forthcoming), geometry of target points can be specified. 
# Preformance testing and error testing has been completed and documented (check documentation). 
# this is a first order method near the boundary and a second order method away from the boundary. 
# to run: 
#     Be sure that a .txt file containing your boundary points is present, then simply type: 
#     python immersed_boundary.py  
# 
# movement of the boundary can be manipulated in the move_ib function within the code (see documentation) 
# a rubber band example, and a moving rod are both simple examples included in documentation.            
#
# Written by: Austin Baird, UNC Chapel Hill. Part of a PhD thesis work on immersed boundary and pumping 
# applications. 
##########################################################################################################

import numpy as np
import subprocess
import os 
import sys 
import matplotlib.pylab as plt
import fluid_constants as f 


#########################################################################################################
# Begin with function definitions
#########################################################################################################


##################################################################################################
# function to move the immersed boundary. Right now it is a moving rod. Array 1/2 hold moving points 
# Q is the number of points along the boundary. 
##################################################################################################
def moveIB(array1, array2, Q):
    #import fluid_constants as f

    amp = f.width/8  #from Lauras test problem
    
    for i in range(f.Q):
        array1[i] = f.width/2 + amp*np.sin(2*np.pi*f.time/f.time_final)
        array2[i] = f.width/4 + i*f.ds

    return array1, array2 

############################################################################
#this computes forces on the boundary due to hookes law stiff*displacement
############################################################################

def target_force(b1, b2, bt1, bt2, fb1, fb2, stiff, Q):

    for i in range(Q):
        fb1[i,0] = stiff*(bt1[i,0] - b1[i,0])   #spring constant*displacement between target and boundary points
        fb2[i,0] = stiff*(bt2[i,0] - b2[i,0])

    #print fb1, fb2
    return fb1, fb2 

#######################################################################################################################
# This code will spread the force generated at the boundary  due to the springs,to neighboring fluid gird points 
# array1/2 are arrays of boundary points (x,y points), fb1/2 are the forces  and Q are the total points. f1/2 are updated
#######################################################################################################################

def forcespread(array1, array2, fb1, fb2, f1, f2, Q):
    
    #import numpy as np
    #import matplotlib.pylab as plt 


    #import fluid_constants as f 

    # first determine which cartesian grid point we are near

    for i in range(Q):
        x1 = np.ceil(array1[i,0]/f.dx)
        x2 = np.ceil(array2[i,0]/f.dx)
        # spred forces
        for k in range(4):
            for h in range(4):
                ii = x1+1-k    
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
                

                f1[ii,jj] = f1[ii,jj]+fb1[i,0]*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x1-k)-array1[i,0])))*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x2-h)-array2[i,0])))*f.ds  

                f2[ii,jj] = f2[ii,jj]+fb2[i,0]*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x1-k)-array1[i,0])))*(1/(4*f.dx))*(1+np.cos((np.pi/(2*f.dx))*(f.dx*(x2-h)-array2[i,0])))*f.ds


    return f1, f2 


################################################################
# This makes the arrays w1, w2 whcih then are transformed to solve
# for fluid velocity
#################################################################

def makew(w1,w2,f1,f2,u,v):
    #import numpy as np
    #import fluid_constants as f 

    #calculate constants 
    d_dx = 1.0/f.dx
    dt_p = f.dt/f.p
    #print u, v
# Comput v values from the peskin paper

    for i in range(f.M):
        for j in range(f.N):    #may want to take away the -1, took away -1  
            # for periodic boundary we need to have special values for our indexing. 
            ip1 = i+1
            jp1 = j+1
            i_1 = i-1
            j_1 = j-1

            if i==f.M-1:  #added a -1
                ip1 = 0

            if j==f.N-1:
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
            
####################################################################
# This file takes in the forces from the boundary and updates 
# pressure (press) and fluid velocity (u,v) accordingly
###################################################################

def fluidsolve(a,b,c,d,f1,f2,u,v,uft,vft,w1,w2,pft,time,press):
    #import fluid_constants as f

    # First thing to do is to calculate w1, w2 and take the fft of them 
    w1, w2 = makew(w1,w2,f1,f2,u,v)

    # After calculating w1,w2 we can now begin taking the fft of our arrays
    w1ft = np.fft.fft2(w1)
    w2ft = np.fft.fft2(w2)
    a = a + 0j
    b = b + 0j
    c = c + 0j
    d = d + 0j
    uft = uft + 0j 
    vft = vft + 0j
    
    # calculate uft, vft from w's, from Peskin paper 
    uft = (w1ft - (a*w1ft + b*w2ft))*d
    vft = (w2ft - (b*w1ft + c*w2ft))*d
    # Now calculate the pressure, need it to be able to store complex values
    pft = pft + 0j

    for i in range(f.M):
        for j in range(f.N):
            if i==0 and j==0 or i==0 and j==f.N/2:
                pft[i,j] = 0
            elif i==f.M/2 and j==0 or i ==f.M/2 and j==f.N/2:
                pft[i,j] = 0
            else: 
                pft[i,j] = ((f.p*f.dx/f.dt)*(np.sin(2*np.pi*i/f.M)*np.imag(w1ft[i,j]) + np.sin(2*np.pi*(j)/f.N)*np.imag(w2ft[i,j])))/((np.sin(2*np.pi*i/f.N)*np.sin(2*np.pi*i/f.M)) + np.sin(2*np.pi*j/f.N)*np.sin(2*np.pi*j/f.N))
                pft[i,j] = pft[i,j] - np.sqrt(np.complex(-1))* ((f.p*f.dx/f.dt)*(np.sin(2*np.pi*i/f.M)*np.real(w1ft[i,j]) + np.sin(2*np.pi*(j)/f.N)*np.real(w2ft[i,j])))/((np.sin(2*np.pi*i/f.N)*np.sin(2*np.pi*i/f.M)) + np.sin(2*np.pi*j/f.N)*np.sin(2*np.pi*j/f.N))
    

# We can now take the inverse transform to get u,v, and press 

    u = np.fft.ifft2(uft)
    v = np.fft.ifft2(vft)
    press = np.fft.ifft2(pft)


# we now want to select the real portions of these arrays, ideally the imaginary portion is small 

    u = np.real(u) 
    v = np.real(v) 
    press = np.real(press)

    return u, v, press 

#############################################################################################################
# This function moves the boundary at the local fluid velocity 
# b1, b2 are the x,y boundary values, and Q is the number of 
# immersed boundary points. b1, b2 are updated and returned
############################################################################################################

def move_boundary(b1, b2, Q, u, v):
    #import numpy as np 
    #import fluid_constants as f 

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


#####################################################################################################
# This function will be used to output the data in .vtk format
# INputs are u/v velocities, vorticity, time step counter k, and boundary points. Each of these will be output 
# in a .vtk file for each 'printing' time step. This will include the appropriate header.
####################################################################################################


def record(u,v,vort,press,f1,f2,b1,b2,k): 

    #import numpy as np 
    #import os
    #import fluid_constants as f 

# first we want to name the files appropriately. 

    save_v = str('velv_binary_%03d' % k) + '.vtk'
    save_u = str('velu_binary_%03d' % k) + '.vtk'
    save_vort = str('vort_binary_%03d' % k) + '.vtk'
    save_pressure = str('press_binary_%03d' % k) + '.vtk'
    save_forcesx = str('forcesx_binary_%03d' % k) + '.vtk'
    save_forcesy = str('forcesx_binary_%03d' % k) + '.vtk'

# now open the files to save the data to. 

    output1 = open(os.path.join("./velv",save_v), 'wb')
    output2 = open(os.path.join("./velu",save_u), 'wb')    
    output3 = open(os.path.join("./vorticity",save_vort), 'wb')
    output4 = open(os.path.join("./pressure",save_pressure), 'wb')
    output5 = open(os.path.join("./forcesx",save_forcesx), 'wb')
    output6 = open(os.path.join("./forcesy",save_forcesy), 'wb')

# begin to output the data into their designated files 
    
    mult = f.M*f.N
    multstr = str(mult)
#header files need to be added first: 

    output1.write('# vtk DataFile Version 3.0' + '\n' + 'example test vtk file' + '\n' + 'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'velv' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 

    output2.write('# vtk DataFile Version 3.0' + '\n'+ 'example test vtk file'+ '\n' +'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'velu' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 

    output3.write('# vtk DataFile Version 3.0' + '\n' + 'example test vtk file'+ '\n'+'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'vorticity' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default'+ '\n') 

    output4.write('# vtk DataFile Version 3.0' + '\n'+ 'example test vtk file'+'\n' +'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'pressure' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 

    output5.write('# vtk DataFile Version 3.0' + '\n' + 'example test vtk file' + '\n' +'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'forcesx' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 

    output6.write('# vtk DataFile Version 3.0' + '\n' + 'example test vtk file' + '\n' +'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'forcesy' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 



#print the data points: 

    for i in range(f.N):
        for j in range(f.M): 
            output1.write(str(v[i,j]))
            output1.write(' ')
            output2.write(str(u[i,j]))
            output2.write(' ')
            output3.write(str(vort[i,j])) 
            output3.write(' ')          
            output4.write(str(press[i,j]))
            output4.write(' ')
            output5.write(str(f1[i,j]))
            output5.write(' ')
            output6.write(str(f2[i,j]))
            output6.write(' ')
            if j==f.M-1:
                output1.write('\n')
                output2.write('\n')
                output3.write('\n')
                output4.write('\n')
                output5.write('\n')
                output6.write('\n')

    print 'wrote', save_u

#----------------------------------------------------------------------------------------#
# compute the vorticity to use in vizualization 
#----------------------------------------------------------------------------------------#

def vorticity(u,v,vort):
    #import fluid_constants as f
    #import numpy as np
    # calculate vorticity away from the boundary

    for i in np.arange(1,f.M-1,1):
        for j in np.arange(1,f.N-1,1):

            vort[i,j] = (v[i+1,j]-v[i-1,j])/(2*f.dx) - (u[i,j+1]-u[i,j-1])/(2*f.dx);

    for j in np.arange(1,f.N-1,1):
        vort[0,j]=(v[1,j]-v[f.M-1,j])/(2*f.dx) - (u[i,j+1]-u[i,j-1])/(2*f.dx); 
        vort[f.M-1,j]=(v[1,j]-v[f.M-1,j])/(2*f.dx) - (u[i,j+1]-u[i,j-1])/(2*f.dx); 


    for i in np.arange(1,f.M-2,1):
        vort[i,0]=(v[i+1,j]-v[i-1,j])/(2*f.dx) - (u[i,1]-u[i,f.N-1])/(2*f.dx); 
        vort[i,f.N-1]=(v[i+1,j]-v[i-1,j])/(2*f.dx) - (u[i,1]-u[i,f.N-1])/(2*f.dx); 

    vort[0,0]=(v[1,0]-v[f.M-1,0])/(2*f.dx) - (u[0,1]-u[0,f.N-1])/(2*f.dx); 
    vort[0,f.N-1]=(v[1,f.N-1]-v[f.M-1,f.N-1])/(2*f.dx) - (u[0,0]-u[0,f.N-1])/(2*f.dx); 

    #print vort
    return vort

#------------------------------------------------------------------------------------------------------# 
# Main routine, you shouldn't need to change anything here, only in the fluid constants section and 
# the move boundary section
#-------------------------------------------------------------------------------------------------------# 


if __name__ == "__main__":

#####################################################################
# Initialize parameters needed for computation
# dt, meshsize and other things may be changed here 
#####################################################################
    #import fluid_constants as f
    #ptime = 0
    k = 0  # frame number

    ###############################################################################
    # Array initialization
    ##############################################################################

    press = np.zeros((f.M,f.N))   #pressure size NXM matrix 
    u = np.zeros((f.M,f.N))   #x-directed velocity 
    v = np.zeros((f.M,f.N))   #y_directed velocity
    f1 = np.zeros((f.M,f.N))   #force in x_direction
    f2 = np.zeros((f.M,f.N))
    t = np.zeros(f.T)   #number of time steps
    mag = np.zeros((f.M,f.N)) #magnitude of velocity vector at each grid point

    #-----------------------Precomputed matrices from paper---------------------#

    w1 = np.zeros((f.M,f.N))
    w2 = np.zeros((f.M,f.N))
    a = np.zeros((f.M,f.N))
    b = np.zeros((f.M,f.N))
    c = np.zeros((f.M,f.N))
    d = np.zeros((f.M,f.N))

    #---------------------Matricies to store fourier transforms-----------------#

    uft = np.zeros((f.M,f.N))
    vft = np.zeros((f.M,f.N))
    pft = np.zeros((f.M,f.N))


    #---------------------Cartesian coordinate matrices--------------------------#

    xcoord = np.zeros((f.M,f.N))   #holds the x-value in the (i,j) position
    ycoord = np.zeros((f.M,f.N))   #holds the y-value in the (i,j) position 
    vort = np.zeros((f.M,f.N))
    #--------------------Boundary arrays and tartget arrays---------------------#


    b1 = np.zeros((f.Q,1))
    b2 = np.zeros((f.Q,1))
    b1t = np.zeros((f.Q,1))   #target points for the immersed boundary 
    b2t = np.zeros((f.Q,1))

    fb1 = np.zeros((f.Q,1))   #storing forces for each boundary point
    fb2 = np.zeros((f.Q,1))


    #########################################################################################
    # Creating directories to store the output files
    ###########################################################################################

    os.makedirs("./velu")
    os.makedirs("./velv")
    os.makedirs("./pressure")
    os.makedirs("./vorticity")
    os.makedirs("./forcesx")
    os.makedirs("./forcesy")

    ##########################################################################
    #Store cartesian values into array, for plotting
    #########################################################################

    for i in range(f.M):
        for j in range(f.N):
            xcoord[i,j] = j*f.dx   #storing values of array_position*dx to get x's cartesian values
            ycoord[i,j] = i*f.dx

    #print xcoord

    #################################################################################
    #Precompute matrices used to calculate u_hat, p_hat from Peskin and McQueen
    # 1996 eqns 14.52 and 14.51
    ##################################################################################

    for i in range(f.M):
        for j in range(f.N):
            if i==0 and j==0 or i==0 and j==f.N/2:
                a[i,j]==0
                b[i,j]==0
                c[i,j]==0
            elif i==f.M/2 and j==0 or i == f.M/2 and j==f.N/2:
                a[i,j] = 0
                b[i,j] = 0
                c[i,j] = 0

            else:
                a[i,j] = np.sin(2*np.pi*i/f.M)*np.sin(2*np.pi*i/f.M)/((np.sin(2*np.pi*i/f.M)*np.sin(2*np.pi*i/f.M)+np.sin(2*np.pi*j/f.N)*np.sin(2*np.pi*j/f.N)))
                b[i,j] = np.sin(2*np.pi*i/f.M)*np.sin(2*np.pi*j/f.N)/((np.sin(2*np.pi*i/f.M)*np.sin(2*np.pi*i/f.M)+np.sin(2*np.pi*j/f.N)*np.sin(2*np.pi*j/f.N)))
                c[i,j] = np.sin(2*np.pi*j/f.N)*np.sin(2*np.pi*j/f.N)/((np.sin(2*np.pi*i/f.M)*np.sin(2*np.pi*i/f.M)+np.sin(2*np.pi*j/f.N)*np.sin(2*np.pi*j/f.N)))

    for i in range(f.M):
        for j in range(f.N):
            d[i,j] = 1/(1+(4*f.nu*f.dt/(f.dx*f.dx))*np.sin(np.pi*i/f.M)*np.sin(np.pi*i/f.M)+np.sin(np.pi*j/f.N)*np.sin(np.pi*j/f.N))

    #print d.shape, a.shape, b.shape, c.shape
    ##################################################################################
    #To follow the matlab code I will be setting the inital boundary here, however
    #future code should just load a .vertex file so that exterior programs can be 
    #used to generate the geometry
    ###################################################################################

    #boudary points for a rod

    for i in range(f.Q):
        b1[i,0] = f.width/2
        b2[i,0] = (f.width/4) + i*f.ds

    #Target points are just the boundary points 

    for i in range(f.Q):
        b1t[i,0] = f.width/2
        b2t[i,0] = (f.width/4) + i*f.ds




    ####################################################################################
    # Time stepping routine (this is the main routine for IB method)
    ####################################################################################


    while f.time < f.time_final: 
        #zero forces 
        f1 = np.zeros((f.M,f.N))
        f2 = np.zeros((f.M,f.N))
        fb1 = np.zeros((f.Q,1))
        fb2 = np.zeros((f.Q,1))

        #print fb1  
        #move target points for the plate (for this example the plate oscillates)

        b1t, b2t = moveIB(b1t,b2t,f.Q)

        #calculate forces from the boundary being applied to the fluid: 

        fb1, fb2 = target_force(b1,b2,b1t,b2t,fb1,fb2,f.stiff,f.Q)

        #spread the force to the fluid: 

        f1, f2 = forcespread(b1,b2,fb1,fb2,f1,f2,f.Q)

        # Now solve for the new velocity and pressure: 

        u,v,press = fluidsolve(a,b,c,d,f1,f2,u,v,uft,vft,w1,w2,pft,f.time,press)

        # Move boundary at the local fluid velocity

        b1, b2 = move_boundary(b1,b2,f.Q,u,v)

        # now we need to output the data:
        #print u, v vorticity, pressure, forces:
        if f.time>f.ptime:
            vort = vorticity(u,v,vort)
            record(u,v,vort,press,f1,f2,b1,b2,k)

            # python can make its own movie of the data but I would reccomend using paraview with the outputted data. 
            #       plt.figure()
            #       cmaxx = np.max(vort)
            #       cmxx = np.max(cmaxx)
            #       plt.pcolor(xcoord, ycoord, vort,cmap = cm.hot, vmax=cmxx/16.0, vmin=-cmxx/16.0)
            #       plt.colorbar()
            
            #       plt.plot(b1,b2, 'yo')
            #       plt.plot(b1t, b2t, 'go')
            
            f.ptime = f.graphtime + f.ptime   #this makes 'print time' increase so taht this if statement doesn't run every time step 
            #print('recording output files')
            #       save = str('%03d' % k) + '.png'
            #       plt.savefig(save, dpi = 100, bbox_inches='tight')
            #       print 'wrote file', save
            #       plt.show()
            #       plt.clf()
            k = k+1   #update frame number
            
        f.time = f.time + f.dt    #update time 
        print('time', f.time, 'ptime',  f.ptime, 'graphtime', f.graphtime)
    #command = ('mencoder',
     #          'mf://*.png',
     #          '-mf',
     #          'type=png:w=800:h=600:fps=15',
     #          '-ovc',
     #          'lavc',
     #          '-lavcopts',
     #          'vcodec=mpeg4',
     #          '-oac',
     #          'copy',
     #          '-o',
     #          'output.avi')

    #os.spawnvp(os.P_WAIT, 'mencoder', command)

    #print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
    #subprocess.check_call(command)

    #print "\n\n The movie was written to 'output.avi'"

    #print "\n\n You may want to delete *.png now.\n\n"
