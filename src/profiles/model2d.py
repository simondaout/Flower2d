#!/opt/local/bin/python2.7

import numpy as np
import math

class fault2d:
    def __init__(self,name,x,y,strike):
        self.name=name
        self.x=x
        self.y=y
        self.strike=strike

class prof:
    def __init__(self,name,x,y,l,w,strike,losmin,losmax,type=None):
        self.name=name
        self.x=x
        self.y=y
        self.l=l
        self.w=w
        self.strike=strike
        self.typ=type
	self.losmin=losmin
	self.losmax=losmax

class topo:
    def __init__(self,name,wdir,filename,color,width,scale=1,topomin=None,topomax=None):
        self.name=name
        self.wdir=wdir
        self.filename=filename
        self.color=color
        self.width=width
        self.scale=scale
        self.topomin=topomin
        self.topomax=topomax
        
        self.yp=[]
        self.xp=[]

    def load(self,xlim=[-1000,1000],ylim=[-1000,1000]):
        fname=file(self.wdir+self.filename)
        x,y,z=np.loadtxt(fname,comments='#',unpack=True,dtype='f,f,f')
        index=np.nonzero((x<xlim[0])|(x>xlim[1])|(y<ylim[0])|(y>ylim[1]))
        self.x,self.y,self.z=np.delete(x,index),np.delete(y,index),np.delete(z,index)*self.scale

class seismi:
    def __init__(self,name,wdir,filename,color,width,scale=1):
        self.name=name
        self.wdir=wdir
        self.filename=filename
        self.color=color
        self.width=width
        self.scale=scale
        
        self.yp=[]
        self.xp=[]

    def load(self,xlim=[-1000,1000],ylim=[-1000,1000]):
        fname=file(self.wdir+self.filename)
        x,y,z,mw=np.loadtxt(fname,comments='#',unpack=True,dtype='f,f,f,f')
        index=np.nonzero((x<xlim[0])|(x>xlim[1])|(y<ylim[0])|(y>ylim[1]))
        self.x,self.y,self.z,self.mw=np.delete(x,index),np.delete(y,index),np.delete(z,index)*self.scale,np.delete(mw,index)

class moho:
    def __init__(self,name,wdir,filename,color,width,scale=1):
        self.name=name
        self.wdir=wdir
        self.filename=filename
        self.color=color
        self.width=width
        self.scale=scale
        
        self.yp=[]
        self.xp=[]

    def load(self,xlim=[-1000,1000],ylim=[-1000,1000]):
        fname=file(self.wdir+self.filename)
        x,y,z=np.loadtxt(fname,comments='#',unpack=True,dtype='f,f,f')
        index=np.nonzero((x<xlim[0])|(x>xlim[1])|(y<ylim[0])|(y>ylim[1]))
        self.x,self.y,self.z=np.delete(x,index),np.delete(y,index),np.delete(z,index)*self.scale


