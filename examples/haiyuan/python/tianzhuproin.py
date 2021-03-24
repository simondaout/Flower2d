# Input Parameters
maindir='../'
outdir=maindir+'output/'

# boundaries map
xmin,xmax=-200,200
ymin,ymax=-200,200

insardata=[
        network(network='T104_spectrum-square_clean_noflata_LOSVelocity_mmyr_s8_km.xy-los',reduction='T104',wdir=maindir+'insar/',dim=1,color='blue',lmin=-3, lmax=3),
        network(network='T333_LOSVelocity_mmyr_s500_km.xy-los',reduction='T333',wdir=maindir+'insar/',dim=1,color='red',lmin=-3, lmax=3),
        ]

gpsdata=[
        network(network='stations_liang_km.dat',reduction='sblock',wdir=maindir+'gps/',dim=2,scale=1,lmin=-6,lmax=6),
        ]

profiles=[
         profile(name='West',x=0,y=0,l=300,w=200,strike=-68,type='distscale'),
        ]

gmtfiles=[
        gmt(name='Fault traces',wdir=maindir+'gmt/',filename='Failles_m_km.xy',color='grey',width=2.),
        gmt(name='1920 Rupture',wdir=maindir+'gmt/',filename='Rupture_km.xy',color='blue',width=4.),
        gmt(name='Seismic Gap',wdir=maindir+'gmt/',filename='Gap_km.xy',width=4.,color='red'),
        ]

topodata=[
        topo(name='SRTM3',wdir=maindir+'gmt/',filename='topo_km.xy-z',color='black',scale=-1),
        ]

fmodel=[
        fault2d(name='Haiyuan',x=-36.,y=10.2),
        ] 


