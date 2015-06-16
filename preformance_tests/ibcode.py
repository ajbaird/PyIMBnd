#####################################################################
# Immersed boundary method python code 
# Written by: Austin Baird, PhD candidate, UNC: Chapel Hill
#####################################################################
 

import numpy as np
import matplotlib.pylab as plt
#import move_boundary as mv_b
#import moveIB  
#import target_force as t_f 
import initial_velocity as iv
#import forcespread as fspread 
import force_calc as fc
import fluidsolve 
import fluid_constants as f 
#import record as rec 
import matplotlib.colors as col
import matplotlib.cm as cm 
#import vorticity as vrt
import subprocess
import os 
import sys 
import time

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
        d[i,j] = 1/(1+(4*f.nu*f.dt/(f.dx*f.dx))*(np.sin(np.pi*i/f.M)*np.sin(np.pi*i/f.M)+np.sin(np.pi*j/f.N)*np.sin(np.pi*j/f.N)))

# store initial velocity to test later in the code: 

u, v, theta_1, theta_2 = iv.initial_velocity(u,v,f.M,f.N,f.dx,f.length,f.time)

#store a copy to use in an error tester later on in the code: 

#testing to see if the initial velocity is the same as matlab code

#print 'u = '
#print u[1,0]
#print 'v = '
#print v[1,0]

#upon testing initial velocities I've dtermined that it's not the problem, most likely the fluid solver...

####################################################################################
# Time stepping routine (this is the main routine for IB method)
####################################################################################

t0 = time.time()
while f.time < f.time_final: 
    
    f1, f2 = fc.force_calc(f1,f2,f.mu,f.dx,f.length,f.p,f.time,f.M,f.N, theta_1, theta_2)

    u,v,press = fluidsolve.fluidsolve(a,b,c,d,f1,f2,u,v,uft,vft,w1,w2,pft,f.time,press)

    #print 'time step done'
   
  #  if f.time>ptime:
       
  #plt.figure()
  #mag = np.sqrt(u**2 + v**2)  # calculate the magitude of velocity
     #   cmaxx = np.max(mag)
      #  cmxx = np.max(cmaxx)
      #  plt.pcolor(xcoord, ycoord, mag,cmap = cm.hot, vmax=cmxx/16.0, vmin=-cmxx/16.0)
      #  plt.colorbar()

    
       # ptime = f.graphtime + ptime   #this makes 'print time' increase so taht this if statement doesn't run every time step 
       # save = str('%03d' % k) + '.png'
       # plt.savefig(save, dpi = 100, bbox_inches='tight')
       # print 'wrote file', save
        #plt.show()
       # plt.clf()
       # k = k+1   #update frame number

    f.time = f.time + f.dt    #update time 
    #print f.time
    print 'time step done'
#create a movie from the .png images generated above: 
t1 = time.time()
print 'computed time is:',  t1-t0 
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

#print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
#subprocess.check_call(command)

#print "\n\n The movie was written to 'output.avi'"

#print "\n\n You may want to delete *.png now.\n\n"

#compute the test fuctions: 

u_test = np.zeros((f.M,f.N))
v_test = np.zeros((f.M,f.N))

u_test, v_test, theta_1, theta_2 = iv.initial_velocity(u_test,v_test,f.M,f.N,f.dx,f.length,f.time)

#I've determinded that the test functions are right (for now) I think I may have been doing it wrong ealier thought 

#print u[0,0]  #this is totally wrong! After one time step, I think that the initial velocity is wrong....
#print v[0,0]
#for i in range(f.M):
    #for j in range(f.N):
    #    u[i,j] = ft.a1[0]*np.sin(ft.theta_1[i,j]) + ft.a2[0]*np.sin(ft.theta_2[i,j])
# flu_err is going to store the L2 error between the computed and exact solution

flu_err = np.zeros((f.M,f.N))
flu_err1 = 0

for i in range(f.M): 
    for j in range(f.N): 
        flu_err[i,j] =  (f.dx)*(f.dx)*np.sqrt((u[i,j]-u_test[i,j])**2 + (v[i,j] - v_test[i,j])**2)
        #flu_err1 = flu_err1 + (f.dx)*(f.dx)*np.sqrt((u[i,j]-u_test[i,j])**2 + (v[i,j] - v_test[i,j])**2)
    #print u 
    #print u_test
err = np.sum(flu_err)
print err
###############################################################

#print out error:
#print flu_err1 
#print flu_err
#plot analytic solution


plt.figure()
flumag = np.sqrt(u_test**2 + v_test**2)
cmaxx = np.max(flumag)
cmxx = np.max(cmaxx)
plt.pcolor(flumag) 
plt.colorbar()
plt.show()


################################################################
#want to calculate the divergence free error: 
################################################################

diverr=0

for i in range(f.M): 
    for j in range(f.N): 
        #for the edge we compute special values 

        jp1 = j+1
        ip1 = i+1
        i_1 = i-1
        j_1 = j-1
        
        if i==f.M-1: 
            ip1 = 0

        if i==0:
            i_1 = f.M-1

        if j==f.N-1:
            jp1 = 0

        if j==0:
            j_1 = f.N-1

        diverr = (1/(2*f.dx))*(u[ip1,j]-u[i_1,j]+v[i,jp1]-v[i,j-1])

print 'diverr='+str(diverr)
