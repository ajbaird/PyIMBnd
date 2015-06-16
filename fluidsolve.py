####################################################################
# This file takes in the forces from the boundary and updates 
# pressure (press) and fluid velocity (u,v) accordingly
###################################################################

def fluidsolve(a,b,c,d,f1,f2,u,v,uft,vft,w1,w2,pft,time,press):
    import numpy as np 
    import matplotlib.pylab as plt 
    import fluid_constants as f
    import makew 

    # First thing to do is to calculate w1, w2 and take the fft of them 

    w1, w2 = makew.makew(w1,w2,f1,f2,u,v)
    
    #print w1

# After calculating w1,w2 we can now begin taking the fft of our arrays

    w1ft = np.fft.fft2(w1)
    w2ft = np.fft.fft2(w2)
    #print w1ft, w2ft

# We can now solve for uft, vft 

    for i in range(f.M):
        for j in range(f.N):
            uft[i,j] = (w1ft[i,j] - (a[i,j]*w1ft[i,j] + b[i,j]*w2ft[i,j]))*d[i,j]
            vft[i,j] = (w2ft[i,j] - (b[i,j]*w1ft[i,j] + c[i,j]*w2ft[i,j]))*d[i,j]
#    print uft.shape  
# Now calculate the pressure

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
