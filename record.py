#####################################################################################################
# This function will be used to output the data in .vtf format
# INputs are u/v velocities, vorticity, time step counter k, and boundary points. Each of these will be output 
# in a .vtk file for each 'printing' time step. This will include the appropriate header.
####################################################################################################


def record(u,v,vort,press,f1,f2,b1,b2,k): 

    import numpy as np 
    import os
    import fluid_constants as f 

# first we want to name the files appropriately. 

    save_v = str('velv_binary_%03d' % k) + '.vtk'
    save_u = str('velu_binary_%03d' % k) + '.vtk'
    save_vort = str('vort_binary_%03d' % k) + '.vtk'
    save_pressure = str('press_binary_%03d' % k) + '.vtk'
    save_forcesx = str('forcesx_binary_%03d' % k) + '.vtk'
    save_forcesy = str('forcesx_binary_%03d' % k) + '.vtk'

# now open the files to save the data to. 

    output1 = open(os.path.join("./velv",save_v), 'wb')
    output2 = open(os.path.join("./velu",save_u), 'wb')    
    output3 = open(os.path.join("./vorticity",save_vort), 'wb')
    output4 = open(os.path.join("./pressure",save_pressure), 'wb')
    output5 = open(os.path.join("./forcesx",save_forcesx), 'wb')
    output6 = open(os.path.join("./forcesy",save_forcesy), 'wb')

# begin to output the data into their designated files 
    
    mult = f.M*f.N
    multstr = str(mult)
#header files need to be added first: 

    output1.write('# vtk DataFile Version 3.0' + '\n' + 'example test vtk file' + '\n' + 'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'velv' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 

    output2.write('# vtk DataFile Version 3.0' + '\n'+ 'example test vtk file'+ '\n' +'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'velu' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 

    output3.write('# vtk DataFile Version 3.0' + '\n' + 'example test vtk file'+ '\n'+'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'vorticity' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default'+ '\n') 

    output4.write('# vtk DataFile Version 3.0' + '\n'+ 'example test vtk file'+'\n' +'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'pressure' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 

    output5.write('# vtk DataFile Version 3.0' + '\n' + 'example test vtk file' + '\n' +'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'forcesx' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 

    output6.write('# vtk DataFile Version 3.0' + '\n' + 'example test vtk file' + '\n' +'ASCII' + '\n' + 'DATASET STRUCTURED_POINTS' + '\n' + 'DIMENSIONS' + ' ' + str(f.M) + ' ' + str(f.N) + ' ' + str(1) + '\n' + 'SPACING' + ' ' + str(f.dx) + ' ' + str(f.dx) + ' ' + str(1) + '\n' + 'ORIGIN' + ' ' +  str(0) + ' ' + str(0) + ' ' + str(0) + '\n' + 'POINT_DATA' + ' ' + multstr + '\n' + 'SCALARS' + ' ' + 'forcesy' + ' ' + 'float' + str(1) + '\n' + 'LOOKUP_TABLE' + ' ' + 'default' + '\n') 



#print the data points: 

    for i in range(f.N):
        for j in range(f.M): 
            output1.write(str(v[i,j]))
            output1.write(' ')
            output2.write(str(u[i,j]))
            output2.write(' ')
            output3.write(str(vort[i,j])) 
            output3.write(' ')          
            output4.write(str(press[i,j]))
            output4.write(' ')
            output5.write(str(f1[i,j]))
            output5.write(' ')
            output6.write(str(f2[i,j]))
            output6.write(' ')
            if j==f.M-1:
                output1.write('\n')
                output2.write('\n')
                output3.write('\n')
                output4.write('\n')
                output5.write('\n')
                output6.write('\n')

    print 'wrote', save_u
