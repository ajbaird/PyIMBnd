#-------------------------------------publishable graph of pumping data----------------------------#
#
# This python code should produce a publishable .eps graph which will be 
# embedded into the latex code without scaling issues. If things seems to 
# to be scaled strangely use the \showthe\columnwidth to determine the appropriate
# width of the figure. then change 'figsize' accordingly. 
#-------------------------------------------------------------------------------------------------#

from numpy import *
import matplotlib.pylab as plt

##unpack two column comma separated data 

div_err,err, dx, N, dt, time  = loadtxt('error.txt', unpack=True)



#plot data

plt.figure(1)
plt.clf
plt.plot(N,err,'ro-', label='L-1 norm')
plt.xlabel(r"Number of grid points (N) ", fontsize=10)
plt.ylabel(r"$\Sigma_{i=1}^{n}err_i$", fontsize=10)
#plt.subplots_adjust(left=0.16)
#plt.subplots_adjust(top=0.1)
plt.title(r"error as a function of N", fontsize=10)    #This needs to change with every graph
plt.grid(True)
#plt.axis([5,55,25,55])   #Change this to fit the data in the graph
plt.savefig('L1norm.eps', bbox_inches='tight')
#plt.show()

plt.figure(2)
plt.clf
plt.plot(N,time,'ro-', label='short')
plt.xlabel(r"Number of grid points (N)", fontsize=10)
plt.ylabel(r"time to compute 5 iterations ($s$)", fontsize=10)
#plt.subplots_adjust(left=0.16)
#plt.subplots_adjust(top=0.1)
plt.title(r"Computation time as a funciton of N", fontsize=10)    #This needs to change with every graph
plt.grid(True)
#plt.axis([5,55,120,300])   #Change this to fit the data in the graph
plt.savefig('time.eps', bbox_inches='tight')
#plt.show()


plt.figure(3)
plt.clf
plt.plot(N,div_err,'bD-', label='short')
#plt.plot(x1,z1,'bD-', label = 'long')
plt.xlabel(r"Number of grid points (N)", fontsize=10)
plt.ylabel(r"Divergence error", fontsize=10)
#plt.subplots_adjust(left=0.16)
#plt.subplots_adjust(top=0.1)
plt.title(r"Divergence error as a function of grid size (N)", fontsize=10)    #This needs to change with every graph
plt.grid(True)
#plt.axis([5,55,120,300])   #Change this to fit the data in the graph
plt.savefig('div_err.eps', bbox_inches='tight')
#plt.show()

#plot fractional excretion: 


