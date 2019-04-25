#!/opt/local/bin/python2.7

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
    """

    def __init__(self,network,reduction,wdir,dim,color='black',scale=1.,theta=False,samp=1,perc=95,lmin=None,lmax=None,
        weight=None):
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

    def loadgps(self):
        gpsf=file(self.wdir+self.network)
        self.name,self.x,self.y=np.loadtxt(gpsf,comments='#',unpack=True,dtype='S4,f,f')
        self.Npoint=len(self.name)
        self.ux,self.uy=np.zeros(self.Npoint),np.zeros(self.Npoint)
        self.sigmax,self.sigmay=np.zeros(self.Npoint),np.zeros(self.Npoint)
        #self.d=np.zeros(self.Npoint*self.dim)
        for j in xrange(self.Npoint):
            station=self.wdir+self.reduction+'/'+self.name[j]
            dated,east,north,esigma,nsigma=np.loadtxt(station,comments='#',usecols=(0,1,2,3,4),unpack=True,dtype='f,f,f,f,f')
            self.ux[j],self.uy[j]=east*self.scale,north*self.scale
            self.sigmax[j],self.sigmay[j]=esigma*self.scale,nsigma*self.scale
        
        if (self.lmin or self.lmax) is None:
           self.lmin = np.nanpercentile(np.array([self.ux,self.uy]), 2)
           self.lmin = np.nanpercentile(np.array([self.ux,self.uy]), 98)

    def loadinsar(self):
        insarf=file(self.wdir+self.network)
        if self.theta is False:
            self.x,self.y,ulos=np.loadtxt(insarf,comments='#',unpack=True,usecols=(0,1,2),dtype='f,f,f')
            self.x,self.y,ulos=self.x[::self.samp],self.y[::self.samp],ulos[::self.samp] 
        else:
            self.x,self.y,ulos,self.los=np.loadtxt(insarf,comments='#',usecols=(0,1,2,3),unpack=True,dtype='f,f,f,f')
            self.x,self.y,ulos,self.los=self.x[::self.samp],self.y[::self.samp],ulos[::self.samp],self.los[::self.samp]
        
        ulos[np.logical_or(ulos==0.0,ulos>9990.)] = np.float('NaN')
        self.ulos=ulos*self.scale
        self.Npoint=len(self.ulos)
        if (self.lmin or self.lmax) is None:
            self.lmin = np.nanpercentile(self.ulos, 2)    
            self.lmax = np.nanpercentile(self.ulos, 98)    
    


