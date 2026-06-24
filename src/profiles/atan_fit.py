#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import sys

def func(yperp, rate, centre, w):
    return  (rate/np.pi)*np.atan2( (yperp - centre),  w)

# a priori values
prior = np.zeros((3))
prior[0]=3.5        # strike-slip rate
prior[1]=0.       # centre of the atan2 profile
prior[2]=30         # characteristic width

# parameters profile
xmin=-50; xmax=50;
ymin=-6; ymax=6;

# load file
files=['P1.txt','P2.txt','P3.txt','P4.txt']

i = 1
for file in files:
  print('Load file: ', file)
  dist, v, std = np.loadtxt(file, unpack=True, comments='#')
  dist = dist*1e-3
  
  # plot data
  fig = plt.figure(2,figsize = (12,9))
  ax = fig.add_subplot(4,1,int(i))
  ax.plot(dist,v,color='dodgerblue',lw=2.)
  ax.plot(dist,v-std,color='dodgerblue',lw=.5)
  ax.plot(dist,v+std,color='dodgerblue',lw=.5)
  #ax.plot(dist, func(dist, prior[0], prior[1], prior[2]), '-r')
  #plt.show()
  #sys.exit()
  
  ###Optimisation
  try:
    pars, cova = opt.curve_fit(func, dist, v, sigma = std, p0=[prior[0],prior[1],prior[2]],  ftol=1e-1, bounds=bounds)   
  except:
    print('No solution found for the estimation')
    sys.exit()

  # print and plot model
  print('Strike-slip rate: {0:6.2f} ± {1:6.2f} mm/yr'.format(pars[0],cova[0,0]))
  print('Centre: {0:6.2f} ± {1:6.2f} km'.format(pars[1],cova[1,1]))
  print('Locking depth: {0:6.2f} ± {1:6.2f} km'.format(pars[2],cova[2,2]))
  ax.plot(dist, func(dist, pars[0], pars[1], pars[2]), '-r', label='rate={0:6.2f} mm/yr, width={1:6.2f} km'.format(pars[0],pars[2]))
  ax.set_ylim([ymin,ymax])
  ax.set_xlim([xmin,xmax])
  ax.grid(axis='y', color='0.95')
  ax.legend(loc='best')
  i += 1
  print()

plt.xlabel('Distance (km)')
plt.ylabel('Strike-slip velocity (mm/yr)')
fig.suptitle("Atan profile modelling")
fig.savefig('fit_profiles.pdf', format='PDF',dpi=150)
plt.show()

