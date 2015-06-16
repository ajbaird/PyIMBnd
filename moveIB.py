#################################################################
# This function will update position of the immersed boundary
# need to make this functional for input data. 
################################################################

def moveIB(array1, array2, Q):
    
    import numpy as np 
    import matplotlib.pylab as plt
    import fluid_constants as f

    amp = f.width/8  #from Lauras test problem
    
    for i in range(f.Q):
        array1[i] = f.width/2 + amp*np.sin(2*np.pi*f.time/f.time_final)
        array2[i] = f.width/4 + i*f.ds

    return array1, array2 
