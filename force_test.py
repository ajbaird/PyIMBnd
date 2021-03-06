# this file will test whether the force is the same before the forcespread, as after. 
import fluid_constants as f
import numpy as np
import forcespread as fspread
import matplotlib.pylab as plt


def force_test(f1,f2,fb1,fb2,Q):
    import fluid_constants as f 
    import numpy as np
    #sum the forces on the boundary: 
    
    #intialize containers for total forces 

    boundary_fx = 0 
    boundary_fy = 0 
    spread_fx = 0
    spread_fy = 0


    #boundary forces 
    boundary_fx = fb1.sum()
    boundary_fy = fb2.sum()

    for i in range(f.M):
        for j in range(f.N):
            spread_fx = spread_fx + f1[i,j]
            spread_fy = spread_fy + f2[i,j]
            
    spread_fx = spread_fx*4*f.ds
    spread_fy = spread_fy*4*f.ds
    
    #print spread_fx, boundary_fx
    print('x-differences' +' \n')
    print  boundary_fx-spread_fx

    print('y-differences' +'\n')
    print boundary_fy-spread_fy 

f1 = np.zeros((f.M, f.N))
f2 = np.zeros((f.M, f.N))

# num of boundary points for the test

Q = 2*f.M 

#print Q
# initialize vectors to use in the test 

force1 = np.zeros((Q,1)) 
force2 = np.zeros((Q,1)) 
b1 = np.zeros((Q,1)) 
b2 = np.zeros((Q,1)) 

#print np.size(b1)
#set arbitrary force values for the boundary forces

for i in range(Q):
    force1[i,0] = np.cos(2*np.pi*(i)/(f.M+2))
    force2[i,0] = np.cos(2*np.pi*(i)/(f.N))
    #print('force1:' + str(force1))
    b1[i,0] = (i)*f.ds
    b2[i,0] = (i)*f.ds
#print(force1)
#plt.plot(force1)
#plt.show()
#print(sum(force1))

#testing for negative positive terms amoung other things 
#countneg = 0
#countpos =0
#count = 0
#for i in range(Q): 
#    if force1[i,0] < 0:
#        countneg = countneg +1
#    if force1[i,0] > 0:
#        countpos = countpos +1
#    count = count + force1[i,0]
#    print force1[i,0]
#    print count 

#print countpos
#print countneg
f1,f2 = fspread.forcespread(b1,b2,force1,force2,f1,f2,Q)

force_test(f1,f2,force1,force2,Q)




