#!/opt/local/bin/python2.7

from network2d import *
from model2d import *
from readgmt import *

# Input Parameters

maindir='../'
outdir=maindir+'output/'


losmin,losmax = -3., 3.
xmin,xmax=-200,200
ymin,ymax=-200,200

fmodel=[
        fault2d(name='bid',x=-36.,y=10.2,strike=-68),
        ] 

insardata=[
        network(network='T104_spectrum-square_clean_noflata_LOSVelocity_mmyr_s8_km.xy-los',reduction='T104',wdir=maindir+'insar/',dim=1,weight=1.,color='blue'),
        network(network='T333_LOSVelocity_mmyr_s500_km.xy-los',reduction='T333',wdir=maindir+'insar/',dim=1,weight=0.3,color='red'),
        ]

gpsdata=[
        network(network='stations_liang_km.dat',reduction='sblock',wdir=maindir+'gps/',dim=2,weight=1.,scale=1),
        ]


gmtfiles=[
        gmt(name='Fault traces',wdir=maindir+'gmt/',filename='Failles_m_km.xy',color='grey',width=2.),
        gmt(name='1920 Rupture',wdir=maindir+'gmt/',filename='Rupture_km.xy',color='blue',width=4.),
        gmt(name='Seismic Gap',wdir=maindir+'gmt/',filename='Gap_km.xy',width=4.,color='red'),
        ]

profiles=[
         prof(name='West',x=0,y=0,l=300,w=200,strike=-68,type='distscale'),
        ]

topodata=[
        topo(name='SRTM3',wdir=maindir+'gmt/',filename='topo_km.xy-z',color='black',width=0.1),
        ]

