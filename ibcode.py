#####################################################################
# Immersed boundary method python code 
# Written by: Austin Baird, PhD candidate, UNC: Chapel Hill
#####################################################################
 

import numpy as np
import matplotlib.pylab as plt
import move_boundary as mv_b
import moveIB  
import target_force as t_f 
import forcespread as fspread 
import fluidsolve 
import fluid_constants as f 
import force_test as ftest
import record as rec 
import matplotlib.colors as col
import matplotlib.cm as cm 
import vorticity as vrt
import subprocess
import os 
import sys 

#####################################################################
# Initialize parameters needed for computation
# dt, meshsize and other things may be changed here 
#####################################################################

ptime = 0
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

    b1t, b2t = moveIB.moveIB(b1t,b2t,f.Q)

    #calculate forces from the boundary being applied to the fluid 

    fb1, fb2 = t_f.target_force(b1,b2,b1t,b2t,fb1,fb2,f.stiff,f.Q)
    
    #spread the force to the fluid 
    
    f1, f2 = fspread.forcespread(b1,b2,fb1,fb2,f1,f2,f.Q)
    #print fb1
    #print fb2
    # Now solve for the new velocity and pressure 
    
    #ftest.force_test(f1,f2,fb1,fb2)

    u,v,press = fluidsolve.fluidsolve(a,b,c,d,f1,f2,u,v,uft,vft,w1,w2,pft,f.time,press)

    #output the data 



    # Move boundary at the local fluid velocity

    b1, b2 = mv_b.move_boundary(b1,b2,f.Q,u,v)
    #print b1, b2
    # now we need to output the data, right now we are only making a vorticity movie 
# this needs to be changed to output full data of u,v,umag,vort, and forces 
    #print u, v 
    if f.time>ptime:
        rec.record(u,v,vort,press,f1,f2,b1,b2,k)
    #make a vorticity plot
        plt.figure()
        vort = vrt.vorticity(u,v,vort)
        cmaxx = np.max(vort)
        cmxx = np.max(cmaxx)
        plt.pcolor(xcoord, ycoord, vort,cmap = cm.hot, vmax=cmxx/16.0, vmin=-cmxx/16.0)
        plt.colorbar()

        plt.plot(b1,b2, 'yo')
        plt.plot(b1t, b2t, 'go')
    
        ptime = f.graphtime + ptime   #this makes 'print time' increase so taht this if statement doesn't run every time step 
        save = str('%03d' % k) + '.png'
        plt.savefig(save, dpi = 100, bbox_inches='tight')
        print 'wrote file', save
        plt.show()
        plt.clf()
        k = k+1   #update frame number

    f.time = f.time + f.dt    #update time 

command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=15',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           'output.avi')

#os.spawnvp(os.P_WAIT, 'mencoder', command)

print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
subprocess.check_call(command)

print "\n\n The movie was written to 'output.avi'"

print "\n\n You may want to delete *.png now.\n\n"


