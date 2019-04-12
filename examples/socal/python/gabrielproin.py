maindir='../'
outdir=maindir+'output/profiles/' 
# boundaries map
xmin,xmax=-300,-100
ymin,ymax=0,250

# type profile:
# default: scatter plot
# std: plot mean and standard deviation insar
# distscale: scatter plot with color scale fct of the profile-parallel distance 
profiles=[
        # prof(name='Big Bend',x=-200,y=137,l=200,w=100,strike=-62.5,type='distscale'),
        # prof(name='Big Bend',x=-200,y=137,l=200,w=100,strike=-62.5,type='std'),
        prof(name='Big Bend',x=-230,y=160,l=200,w=20,strike=-62.5),
        prof(name='Big Bend',x=-196,y=144,l=200,w=20,strike=-62.5),
        prof(name='Big Bend',x=-160,y=125,l=200,w=20,strike=-62.5),
        ]

# InSAR data (defined in networkopti.py)
insardata=[
        network(network='t170_cmyr_km.xylos',reduction='T170',wdir=maindir+'insar/',dim=1,scale=-10.,color='red',samp=2)
        ]

# Plot GPS data not implemented yet: not co
gpsdata=[
#         network(network='cgps_stations_km.dat',reduction='PBO',wdir=maindir+'gps/',dim=3,scale=1.),
#         #network(network='cgps_stations_km.dat',reduction='PBO',wdir=maindir+'gps/',dim=3,weight=1.,scale=1.,\
#         # errorfile='../../output/gps/PBOresidus_0.psvelo'),
         ]

# topographic file
topodata=[
        topo(name='SRTM3',wdir=maindir+'gmt/',filename='srtm_ell_nan_km.xyz',color='black'),
     ]

# optionnal: add data for plots in gmt format (>)
gmtfiles=[
        gmt(name='fault traces',wdir=maindir+'gmt/',filename='ca_faults_north_km.xyz',color='grey',width=0.5),
        gmt(name='coast lines',wdir=maindir+'gmt/',filename='ca_coasts_north_km.xyz',color='grey',width=1.),
        ]

fmodel=[
        fault2d(name='SAF',x=-177.9,y=134.9),
        ] 

