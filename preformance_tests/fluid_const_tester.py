import numpy as np
import fluid_constants as f

a1=np.zeros(2)
a2=np.zeros(2)
bk1=np.zeros(2)
bk2=np.zeros(2)
k1=np.zeros(2)
k2=np.zeros(2)

a1[0]=-2
a1[1]=1
a2[0]=2
a2[1]=2

bk1[0]=1
bk1[1]=2
bk2[0]=-1
bk2[1]=1

theta_1 = np.zeros((f.M, f.N))
theta_2 = np.zeros((f.M, f.N))

k1[0]=(2*np.pi*bk1[0]/f.length)
k1[1]=(2*np.pi*bk1[1]/f.length)
k2[0]=(2*np.pi*bk2[0]/f.length)
k2[1]=(2*np.pi*bk2[1]/f.length)

b_1=1
b_2=1

w_1=0.1
w_2=0.1

print a1[1]
