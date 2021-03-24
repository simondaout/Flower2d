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
        profile(name='Big Bend',x=-200,y=137,l=200,w=100,strike=-62.5,type='std',lbins=10.),
        ]

# InSAR data (defined in networkopti.py)
insardata=[
        network(network='t170_cmyr_km.xylos',reduction='T170',wdir=maindir+'insar/',dim=1,scale=-10.,color='blue',samp=10)
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
        fault2d(name='SAF',x=-177.9,y=134.9,strike=-62.5),
        ] 

