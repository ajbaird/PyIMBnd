# this file will plot the computed error in the immersed boundary method: 

import numpy as np 



#load the data 

err = np.loadtxt('err')


##graph parameters 

fig_width_pt = 360.0      #this is the width of the figure in latex, make sure its up to date to avoid scaling issues. 
inches_in_pt = 1.0/72.27      #This computes inches in a single point 
golden_mean = (sqrt(5)-1.0)/2.0    #This makes the height to width ratio pretty
fig_width = fig_width_pt*inches_in_pt    #width of figure
fig_height = golden_mean*fig_width    # height of figure
fig_size = [fig_width,fig_height]    #putting these in an array 

#change parameters

plt.rc('font', family='serif')
params = {'backend': 'ps',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': True,
          'figure.figsize': fig_size}

plt.rcParams.update(params)

#plot data

plt.figure(1)
plt.clf
