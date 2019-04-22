# path to gps insar gmt ... directories
maindir='../'
# path to output files
outdir=maindir+'output/'

# Minimisation Parameters
niter=10000# number of iterations
nburn=5000# numbert of burned iterations to not be influenced by itinial values

# Model Parameters (defined in modelopti.py) 
inv=inversion(name='SAF',
        structures=[
            # define one main fault and various optional secondary faults
            # ss: strike-slip, short: shortening, w: depth
            # D: horizontal distance to the main fault, 
            # H: vertical distance to the main fault

            # Simple ramp-decollement structure
            # mainfault(
            #       name='SAF', ss=-26.,sigmass=10,
            #       short=7.,sigmashort=3.,w=20.,sigmaw=20.,dip=0.,
            #       distss='Unif',distds='Unif',distH='Unif',
            #       ),
            # ramp(
            #       name='Ramp',ss=0.,sigmass=0.,D=-50,sigmaD=40.,H=20.,sigmaH=20
            #     )

            # Main flower structure
            mainflower(
                name='SAF',ss=-27.,sigmass=10,short=7.,sigmashort=3.,w=15.,sigmaw=10.,dip=0.,
                name3='Kink',ss3=0.,sigmass3=0.,D3=0.,sigmaD3=20.,H3=10.,sigmaH3=10., 
                name2='SG',ss2=0.,sigmass2=0.,D2=-20.,sigmaD2=20.,H2=10.,sigmaH2=10.,
                # optional: define prior distribution for each parameters (default: 'Unif')
                # distss='Gaus',distds='Gaus',distH='Unif',distD='Unif',
                ),
            # We add another flower struct. to the main one
            flower(
                name2='SM',ss2=0,sigmass2=0,D2=-5.,sigmaD2=5.,H2=5.,sigmaH2=5.,
                name1='Decol',ss1=0.,sigmass1=0.,D1=-25.,sigmaD1=20.,H1=0.,sigmaH1=0.,
                ),
            # We add a frontal ramp 
            #ramp(name='PHT',ss=-6.,sigmass=5.,D=-15.,sigmaD=10.,H=5.,sigmaH=5.) 
            ],
            # azimuth of the 2d model: azimuth of the main fault
            strike=-62.5,
            # l: lenght, w: width, proj: Envisat proj los to east,north,z coordinates (for a mean incidence angle)
            profile=profile(name='San Gabriel',x=-177.9,y=134.9,l=200,w=60,proj=[0.382, -0.0811, 0.92]),
	    fullcov=False,
            )

# GPS data (defined in networkopti.py)
gpsdata=[
        network(network='cgps_stations_km.dat',reduction='PBO',wdir=maindir+'gps/',dim=3,weight=.02,scale=1.,plotName=False),
        ]

# InSAR data (defined in networkopti.py)
insardata=[
        network(network='t170_cmyr_s100_km.xy-los',reduction='T170',wdir=maindir+'insar/',dim=1,weight=0.4,scale=10.,color='blue',cov=[1.3702141838,1.40570896897,4.29813868204])
        ]

# optionnal: add data for plots in gmt format (>)
gmtfiles=[
        #gmt(name='fault traces',wdir=maindir+'gmt/',filename='ca_faults_north_km.xyz',color='grey',width=0.5),
        gmt(name='coast lines',wdir=maindir+'gmt/',filename='ca_coasts_north_km.xyz',color='grey',width=1.),
        ]

# optionnal: add data for plots in gmt format (>)
plotdata=[
        #seismi(name='2002-2011 seismicity (Hauksson, 2012)',wdir=maindir+'seismicity/',filename='sc_02_11_km.xydm',color='orange',width=.5),
        #moho(name='Moho (Tape 2012)',wdir=maindir+'gmt/',filename='cal_moho01_q02_q08_ir02_id01_km.xyz',color='red',width=2.),
     ]

#topographic plot
topodata = [
        topo(name='SRTM3',wdir=maindir+'gmt/',filename='srtm_ell_nan_km.xyz',color='black',width=0.001),
        ]
