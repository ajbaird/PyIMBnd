# compute the vorticity to use in vizualization 

def vorticity(u,v,vort):
    import fluid_constants as f
    import numpy as np
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
