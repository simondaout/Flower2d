import numpy as np
import math

class network:
    """ 
    Network class: Load InSAR or GPS data 
    Parameters: 
    network: name input text file
    reduction: reduction name for plot
    wdir: relative path input file
    dim: 1=InSAR, 2,3=GPS
    color: plot option, default: 'black' 
    scale: scale option, default: 1
    theta: load insicence angle in 4th column and project los to average incidence angle
    assuming horizontal displacements, default: False
    samp: subsample option, default:1 
    perc: cleaning outliers option within bins profile, default: percentile=95
    lmin,lmax: min max options for plots
    utm_proj: EPSG UTM projection. If not None, project data from WGS84 to EPSG.
    ref: [lon, lat] reference point (default: None). 
    prof=[east, north, up] optional projection into average LOS vector
    """

    def __init__(self,network,reduction,wdir,dim,color='black',scale=1.,theta=False,\
        samp=1,perc=95,lmin=None,lmax=None,plotName=None, utm_proj=None, ref=None, cst=0, proj=None):

        self.network=network
        self.reduction=reduction
        self.wdir=wdir
        self.dim=dim
        self.color=color
        self.scale=scale
       
        self.Npoint=0.
        self.x,self.y=[],[]
        # Data
        self.ux,self.uy=[],[]
        self.sigmax,self.sigmay=[],[]
        self.ulos=[]
        self.upar,self.uperp=[],[]
        #self.d=[]
        self.theta=theta
        self.samp=samp
        self.perc=perc

        self.lmin = lmin
        self.lmax = lmax
        self.plotName = plotName

        self.cst = cst

        # projection
        self.utm_proj=utm_proj
        self.ref=ref
        if self.utm_proj is not None:
            import pyproj
            self.UTM = pyproj.Proj("EPSG:{}".format(self.utm_proj))
            if self.ref is not None:
                self.ref_x,self.ref_y =  self.UTM(self.ref[0],self.ref[1])
            else:
                self.ref_x,self.ref_y = 0,0

        # projection to LOS
        self.proj = proj

    def loadgps(self):
        gpsf=self.wdir+self.network
        if self.utm_proj is None:
            self.name,self.x,self.y=np.loadtxt(gpsf,comments='#',unpack=True,dtype='S4,f,f')
            # convert to meters
            self.x,self.y = self.x*1e3,self.y*1e3 
        else:
            self.name,self.lon,self.lat=np.loadtxt(gpsf,comments='#',unpack=True,dtype='S4,f,f')
            self.x, self.y = self.UTM(self.lon, self.lat) 
            self.x, self.y = (self.x - self.ref_x), (self.y - self.ref_y)

        self.Npoint=len(self.name)
        
        if self.dim == 2:
            self.ux,self.uy=np.zeros(self.Npoint),np.zeros(self.Npoint)
            self.sigmax,self.sigmay=np.zeros(self.Npoint),np.zeros(self.Npoint)
            #self.d=np.zeros(self.Npoint*self.dim)
            for j in range(self.Npoint):
                station=self.wdir+self.reduction+'/'+self.name[j].decode('utf-8')
                dated,east,north,esigma,nsigma=np.loadtxt(station,comments='#',usecols=(0,1,2,3,4),unpack=True,dtype='f,f,f,f,f')
                self.ux[j],self.uy[j]=east*self.scale,north*self.scale
                self.sigmax[j],self.sigmay[j]=esigma*self.scale,nsigma*self.scale
        elif self.dim == 3:
            self.ux,self.uy, self.uv=np.zeros(self.Npoint),np.zeros(self.Npoint),np.zeros(self.Npoint)
            self.sigmax,self.sigmay,self.sigmav=np.zeros(self.Npoint),np.zeros(self.Npoint),np.zeros(self.Npoint)
            self.ulos,self.sigmalos = np.zeros(self.Npoint),np.zeros(self.Npoint)
            for j in range(self.Npoint):
                station=self.wdir+self.reduction+'/'+self.name[j].decode('utf-8')
                dated,east,north,up,esigma,nsigma,upsigma=np.loadtxt(station,comments='#',usecols=(0,1,2,3,4,5,6),unpack=True,dtype='f,f,f,f,f,f,f')
                self.ux[j],self.uy[j],self.uv[j]=east*self.scale,north*self.scale,up*self.scale
                self.sigmax[j],self.sigmay[j],self.sigmav[j]=esigma*self.scale,nsigma*self.scale,upsigma*self.scale
                if self.proj is not None:
                    self.ulos[j] = self.ux[j]*self.proj[0]+self.uy[j]*self.proj[1]+self.uv[j]*self.proj[2]
                    self.sigmalos[j] = self.sigmax[j]*self.proj[0]+self.sigmay[j]*self.proj[1]+self.sigmav[j]*self.proj[2]

        else:
            print('Error: GNSS dimension is not 2 or 3! Exit')
            sys.exit()

        if ((self.lmin or self.lmax) is None):
             self.lmin = np.min(np.array([self.ux,self.uy])) - 1
             self.lmax = np.max(np.array([self.ux,self.uy])) + 1

    def loadinsar(self):
        insarf=self.wdir+self.network
        if self.utm_proj is None:
            if self.theta is False:
                self.x,self.y,ulos=np.loadtxt(insarf,comments='#',unpack=True,usecols=(0,1,2),dtype='f,f,f')
                # convert to meters
                self.x,self.y,ulos=self.x[::self.samp]*1e3,self.y[::self.samp]*1e3,ulos[::self.samp] 
            else:
                self.x,self.y,ulos,self.los=np.loadtxt(insarf,comments='#',usecols=(0,1,2,3),unpack=True,dtype='f,f,f,f')
                self.x,self.y,ulos,self.los=self.x[::self.samp],self.y[::self.samp],ulos[::self.samp],self.los[::self.samp]
        else:
            if self.theta is False:
                self.lon,self.lat,ulos=np.loadtxt(insarf,comments='#',unpack=True,usecols=(0,1,2),dtype='f,f,f')
                self.lon,self.lat,ulos=self.lon[::self.samp],self.lat[::self.samp],ulos[::self.samp] 
            else:
                self.lon,self.lat,ulos,self.los=np.loadtxt(insarf,comments='#',usecols=(0,1,2,3),unpack=True,dtype='f,f,f,f')
                self.lon,self.lat,ulos,self.los=self.lon[::self.samp],self.lat[::self.samp],ulos[::self.samp],self.los[::self.samp]
            
            self.x, self.y = self.UTM(self.lon, self.lat)
            self.x, self.y = (self.x - self.ref_x), (self.y - self.ref_y)

        # ulos[np.logical_or(ulos==0.0,ulos>9990.)] = np.float('NaN')
        self.ulos=ulos*self.scale + self.cst
        self.Npoint=len(self.ulos)   
    


