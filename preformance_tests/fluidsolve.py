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
    
    #for j in range(f.N):
    #   print 'w1 = '
    #   print w1[1,j]
    #   print 'w2 = '
    #   print w2[1,j]
    #fairly sure that the makew fuction is working properly and that w1 and w2 are right
    # After calculating w1,w2 we can now begin taking the fft of our arrays
     
    #print w1
    w1ft = np.fft.fft2(w1)
    w2ft = np.fft.fft2(w2)
    #test1 = np.real(w1ft)
    #test = np.fft.ifft(w1ft,axis=0)
    #print 'after the transform'
    #print w2ft
    #print test

    # The transformed arrays are now the same at matlab! onto the u and v functions and their inverse transform!
    

    #for j in range(f.N):
     #   print 'w1ft = '
      #  print w1ft[1,j]
      #  print 'w2ft = '
      #  print w2ft[1,j]
      #needd to check after an fft is implimented, something is wrong with the fft...
      #print 'wft1 = '
      #print w1ft[0,0]
      #print 'wft2 = '
      #print w2ft[0,0]
      # We can now solve for uft, vft 
      #need to test every component in this computation: a check, b check,c check, d check, w1ft check, w2ft check
      
      #Missing complex information here, need to convert everything to a complex array...
    a = a + 0j
    b = b + 0j
    c = c + 0j
    d = d + 0j
    uft = uft + 0j 
    vft = vft + 0j
    #print a
    
    #for i in range(f.M):
        #for j in range(f.N):
           # uft[i,j] = (w1ft[i,j] - (a[i,j]*w1ft[i,j] + b[i,j]*w2ft[i,j]))*d[i,j]
           # vft[i,j] = (w2ft[i,j] - (b[i,j]*w1ft[i,j] + c[i,j]*w2ft[i,j]))*d[i,j]

    #testing preformance of numpy arrays

    uft = (w1ft - (a*w1ft + b*w2ft))*d
    vft = (w2ft - (b*w1ft + c*w2ft))*d
    # for now I'm concluding that uft and  vft are correct, if everything still doesnt work in the end, come back to this  
    # there is def a problem with the pressure terms. No doens't seem so....
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
