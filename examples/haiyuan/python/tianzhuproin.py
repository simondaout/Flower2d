# Input Parameters
maindir='../'
outdir=maindir+'output/'

# boundaries map
xmin,xmax=-200,200
ymin,ymax=-200,200

#Network class: Load InSAR or GPS data
#Parameters:
#network: name input text file
#reduction: reduction name for plot
#wdir: relative path input file
#dim: 1=InSAR, 2,3=GPS
#color: plot option, default: 'black'
#scale: scale option, default: 1
#theta: load insicence angle in 4th column and project los to average incidence angle
#assuming horizontal displacements, default: False
#samp: subsample option, default:1
#perc: cleaning outliers option within bins profile, default: percentile=95
#lmin,lmax: min max options for plots
insardata=[
        network(network='T104_spectrum-square_clean_noflata_LOSVelocity_mmyr_s8_km.xy-los',reduction='T104',wdir=maindir+'insar/',dim=1,color='blue',lmin=-3, lmax=3),
        network(network='T333_LOSVelocity_mmyr_s500_km.xy-los',reduction='T333',wdir=maindir+'insar/',dim=1,color='red',lmin=-3, lmax=3),
        ]
gpsdata=[
        network(network='stations_liang_km.dat',reduction='sblock',wdir=maindir+'gps/',dim=2,scale=1),
        ]

# prof class: Load profiles
# Parameters:
# name: name profile
# x,y: reference point
# l,w: length, width progile
# strike: strike fault
# type:  std - plot mean and standard deviation InSAR;
# distscale - scatter plot with color scale function of the profile-parallel distance;
# stdscat - plot scatter + standar deviation.
profiles=[
         prof(name='West',x=0,y=0,l=300,w=200,strike=-68,type='distscale'),
        ]

#Class gmt: read GMT file: text file with object separated with >
#filename: name input file
#wdir: path input file
#name: name for plot
#color, width options for plot
gmtfiles=[
        gmt(name='Fault traces',wdir=maindir+'gmt/',filename='Failles_m_km.xy',color='grey',width=2.),
        gmt(name='1920 Rupture',wdir=maindir+'gmt/',filename='Rupture_km.xy',color='blue',width=4.),
        gmt(name='Seismic Gap',wdir=maindir+'gmt/',filename='Gap_km.xy',width=4.,color='red'),
        ]

# topo class: Load topographic file
# Parameters:
# filename: name input file
# name: name file for plot
# wdir: path input file
# scale: scale values
# color
# topomin,topomax
# plotminmax: option to also plot min max topo within bins
topodata=[
        topo(name='SRTM3',wdir=maindir+'gmt/',filename='topo_km.xy-z',color='black'),
        ]

# fault2d class: Load 2D fault for plot only
# help to position fault for futher modeling
# Parameters:
# name: name fault
# x,y:  position
fmodel=[
        fault2d(name='Haiyuan',x=-36.,y=10.2),
        ] 


