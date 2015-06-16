###############################################################################################
#Store fluid variables used in all the definitions Change these variables to control fluid dynamics
############################################################################################

def fluid_constants():
#spatial grid size
#For this version of the code, M should equal N, and length should equal
#width
    N = 32
    M = 32

    Q = M

    #Mechanical variables
    mu = 3.0  #dynamic viscosity
    p= 2.0    #density
    nu = mu/p #kinematic viscosity
    stiff = 3.0   #stiffness of springs between boundary points
    
    #Spatial variables
    length = 0.4  #length of domain
    width = 0.4   #width of domain
    dx = length/N        #absolute length of grid cell  
    ds = length/(2*N)   #boundary step  
    dt = dx*dx/nu      #time step  

    #new timing parameter values
    time_final = 10*dt 	#total time
    L = 5 #number of times to graph velocity and pressure
    graphtime = time_final/L   #graphs every dptime  
    vel_time = time_final/(100*20)   #output forces? 
    time = 0
    T = np.ceil(time_final/dt)
    ptime = 0   #printing time

