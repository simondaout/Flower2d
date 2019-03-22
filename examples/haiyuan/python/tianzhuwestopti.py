#!/usr/bin/env python2.7
from modelopti import *
from flatten import *
from networkopti import *
from readgmt import *
import os

# path to gps insar gmt ... directories
maindir='../'
# path to output files
outdir=maindir+'output/'

# Minimisation Parameters
niter=10000 # number of iterations
nburn=5000  # numbert of burned iterations to not be influenced by itinial values

# Model Parameters (defined in modelopti.py) 
inv=inversion(name='Haiyuan Fault System',
          structures=[
              #######################
              # Option1: maindecol + ramp : simple ramp-decollement system
              #######################
              maindecol(
              # define one main fault and various optional secondary faults
              # ss: strike-slip, short: shortening, w: depth
              name='JQH',x=-36.,y=10.2,ss=10.,sigmass=10.,short=-5.,
              sigmashort=5.,w=15.,sigmaw=15,dip=180.,
              # optional: define prior distribution for each parameters (default: 'Unif')
              distss='Unif',distshort='Unif',distH='Unif'
              ),
              # D: horizontal distance to the main fault, H: vertical distance to the main fault
              ramp(
                name='QT',ss=0.,sigmass=0.,D=30,sigmaD=25.,H=15.,sigmaH=15
                ),

              #######################
              # Option 2: flower structure
              ####################### 
              # mainflower(
              #   name='JQH',x=-36.,y=10.2,sstot=10.,sigmasstot=10.,short=-5.,sigmashort=5.,w=15.,sigmaw=15.,dip=180.,
              #   name2='QT',ss2=0.,sigmass2=0.,D2=30,sigmaD2=25.,H2=15.,sigmaH2=15,
              #   name3='Creeping',ss3=10.,sigmass3=10.,D3=0,sigmaD3=0.,H3=15.,sigmaH3=15,
              #   distss='Unif',distshort='Unif',distH='Unif'
              # ),
                ],
          # azimuth of the 2d model: azimuth of the main fault
          strike=-68,
          # l: lenght, w: width, proj: Envisat proj los to east,north,z coordinates (for a mean incidence of 20.5)
          profiles=profile(name='West',x=-38,y=10,l=300,w=50,proj=[0.318, -0.0827, 0.9396]), 
          )

# GPS data (defined in networkopti.py)              
gpsdata=[
        # network: east, north gps station, reduction: dir to gps vectors, wdir: path to reduction, dim: dimention gps vectors (2/3)
        # weight: weight to gps vectors, scale: convert unity to mm if necessary
        network(network='stations_liang_km.dat',reduction='sblock',wdir=maindir+'gps/',dim=2,weight=1.,scale=1),
        ]

# InSAR data (defined in networkopti.py) 
insardata=[
        network(network='T104_spectrum-square_cleanRMS_noflata_LOSVelocity_mmyr_s200_km.xy-los',reduction='T104',wdir=maindir+'insar/',dim=1,weight=1.),
        # residual of a first inversion to compute covariance
        #network(network='T104residus.xylos',reduction='T104',wdir=maindir+'output/insar/',dim=1,weight=1.),
        ]

# optional volumic deformation (defined in modelopti.py)
volumic = [] 

# optionnal: add data for plots in gmt format (>)
gmtfiles=[
        gmt(name='Fault traces',wdir=maindir+'gmt/',filename='Failles_m_km.xy',color='grey',width=2.),
        gmt(name='1920 Rupture',wdir=maindir+'gmt/',filename='Rupture_km.xy',color='blue',width=4.),
        gmt(name='Seismic Gap',wdir=maindir+'gmt/',filename='Gap_km.xy',width=4.,color='red'),
        ]

#optionnal: add data for plots
plotdata=[
        seismi(name='seismicity',wdir=maindir+'seismicity/',filename='usgs_2000-2014_101-105_34-39.xydm',color='orange',width=0.5),
        ]

#topographic plot
topodata = [
  topo(name='SRTM3',wdir=maindir+'gmt/',filename='topo_km.xy-z',color='black',width=0.1)
]

# Parameters covariance  (defined in networkopti.py) 
fullcov='no' # if no: diagonal only --> try this first
             # if yes: create full cov from [sudhaus&jonhson09] 
# optional mask 
# maskcov = [-10,-10,20,20]
# param=[None,None,None] # if None: compute [sill, sigm, lamb], if not: go to inversion
# rampcov='lin' # lin/quad/cub optional ramp

#los coordinates for figure (deprecated?)
x1,y1=50,-50

