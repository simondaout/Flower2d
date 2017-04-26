#!/opt/local/bin/python2.7

import numpy as np
import math

class network:
    def __init__(self,network,reduction,wdir,dim,color='black',weight=1.,scale=1.,theta=False,samp=1):
        self.network=network
        self.reduction=reduction
        self.wdir=wdir
        self.dim=dim
        self.color=color
        self.sigmad=1./weight
        self.scale=scale
       
        self.Npoint=0.
        self.x,self.y=[],[]
        # Data
        self.ux,self.uy=[],[]
        self.sigmax,self.sigmay=[],[]
        self.ulos=[]
        self.upar,self.uperp=[],[]
        #self.d=[]
        #self.sigmad=[]
        self.theta=theta
        self.samp=samp

        # Model
        self.mx,self.my=[],[]
        self.mlos=[]
        self.mpar,self.mperp=[],[]

    def loadgps(self):
        gpsf=file(self.wdir+self.network)
        name,x,y=np.loadtxt(gpsf,comments='#',unpack=True,dtype='S4,f,f')
        index=np.nonzero((x<xlim[0])|(x>xlim[1])|(y<ylim[0])|(y>ylim[1]))
        self.name,self.x,self.y=np.delete(name,index),np.delete(x,index),np.delete(y,index)
        self.Npoint=len(self.name)
        self.ux,self.uy=np.zeros(self.Npoint),np.zeros(self.Npoint)
        self.sigmax,self.sigmay=np.zeros(self.Npoint),np.zeros(self.Npoint)
        #self.d=np.zeros(self.Npoint*self.dim)
        for j in xrange(self.Npoint):
            station=self.wdir+self.reduction+'/'+self.name[j]
            dated,east,north,esigma,nsigma=np.loadtxt(station,comments='#',usecols=(0,1,2,3,4),unpack=True,dtype='f,f,f,f,f')
            self.ux[j],self.uy[j]=east*self.scale,north*self.scale
            self.sigmax[j],self.sigmay[j]=esigma*self.scale,nsigma*self.scale
            #self.d[self.dim*j]=east
            #self.d[self.dim*j+1]=north

    def loadinsar(self):
        insarf=file(self.wdir+self.network)
        if self.theta is False:
            self.x,self.y,self.ulos=np.loadtxt(insarf,comments='#',unpack=True,dtype='f,f,f')
            self.x,self.y,self.ulos=self.x[::self.samp],self.y[::self.samp],self.ulos[::self.samp] 
        else:
            self.x,self.y,self.ulos,self.los=np.loadtxt(insarf,comments='#',unpack=True,dtype='f,f,f,f')
            self.x,self.y,self.ulos,self.los=self.x[::self.samp],self.y[::self.samp],self.ulos[::self.samp],self.los[::self.samp]
        self.ulos=self.ulos*self.scale
        self.Npoint=len(self.ulos)
    


