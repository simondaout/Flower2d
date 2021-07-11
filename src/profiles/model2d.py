import numpy as np
import math
import pyproj
import sys
import pandas

class fault2d:
    """ 
    fault2d class: Load 2D fault for plot only
    help to position fault for futher modeling
    Parameters: 
    name: name fault
    x,y: position east, north
    """

    def __init__(self,name,x,y,strike=None,utm_proj=None, ref=None):
        self.name=name
        lon,lat = x,y
        
        # projection
        self.utm_proj=utm_proj
        self.ref=ref
        if self.utm_proj is not None:
            self.UTM = pyproj.Proj("EPSG:{}".format(self.utm_proj))
            if self.ref is not None:
                self.ref_x,self.ref_y =  self.UTM(self.ref[0],self.ref[1])
            else:
                self.ref_x,self.ref_y = 0,0
        
        if self.utm_proj is None:
            self.x,self.y = lon*1e3,lat*1e3 
            if (lon is None) or (lat is None):
                print('utm_proj is not defined, you must defined ref points (x,y) in UTM. Exit!')
                sys.exit()
        else:
            print('Read reference point profile in lat/lon')
            x, y = self.UTM(lon, lat) 
            self.x,self.y=(x-self.ref_x),(y-self.ref_y)
 
        if strike is not None:
          if strike > 0 :
            self.strike=strike-180
          else:
            self.strike=strike

class profile:
    """ 
    profile class: Load profiles 
    Parameters: 
    name: name profile
    x,y: reference point 
    l,w: length, width progile
    strike: strike profile
    type:  std - plot mean and standard deviation InSAR;
    distscale - scatter plot with color scale function of the profile-parallel distance;
    stdscat - plot scatter + standar deviation. 
    flat: if not None estimate a ramp along profile. lin: linear ramp, quad: quadratic, cub: cubic.
    If number InSAR network is 2 then estimate ramp within the overlaping area (Default: None)
    lbins: larger bins for profile
    loc_ramp: location ramp estimation. Can be positive (for postive distances along profile) or negative. Default: None
    """

    def __init__(self,name,l,w,strike,type=None,
        flat=None,lbins=None,loc_ramp=None,x=None,y=None,lat=None,lon=None,utm_proj=None, ref=None):
        self.name=name
        self.x=x
        self.y=y
        self.l=l*1e3
        self.w=w*1e3
        self.flat=flat
        self.lbins=lbins*1e3
        self.loc_ramp=loc_ramp

        if strike > 0:
            self.strike=strike-180
        else:
            self.strike=strike

        self.typ=type
        # lmin,lmax needs to be an attribute of network because different plots for both

        # projection
        self.utm_proj=utm_proj
        self.ref=ref
        if self.utm_proj is not None:
            self.UTM = pyproj.Proj("EPSG:{}".format(self.utm_proj))
            if self.ref is not None:
                self.ref_x,self.ref_y =  self.UTM(self.ref[0],self.ref[1])
            else:
                self.ref_x,self.ref_y = 0,0
        
        if self.utm_proj is None:
            self.x,self.y = x*1e3,y*1e3 
            if (x is None) or (y is None):
                print('utm_proj is not defined, you must defined ref points (x,y) in UTM. Exit!')
                sys.exit()
        else:
            print('Read reference point profile in lat/lon')
            x, y = self.UTM(lon, lat) 
            self.x,self.y=(x-self.ref_x),(y-self.ref_y)
    
class topo:
    """ 
    topo class: Load topographic file 
    Parameters: 
    filename: name input file
    name: name file for plot
    wdir: path input file
    scale: scale values
    color
    topomin,topomax
    plotminmax: option to also plot min max topo within bins
    """

    def __init__(self,name,filename,wdir,color='black',scale=1,topomin=None,topomax=None,plotminmax=False, 
        width=1.,utm_proj=None, ref=None, axis=None):
        self.name=name
        self.filename=filename
        self.wdir=wdir
        self.color=color
        self.scale=scale
        self.topomin=topomin
        self.topomax=topomax
        self.plotminmax=plotminmax
        self.width = width
        self.axis = axis
        
        self.yp=[]
        self.xp=[]

        # projection
        self.utm_proj=utm_proj
        self.ref=ref
        if self.utm_proj is not None:
            self.UTM = pyproj.Proj("EPSG:{}".format(self.utm_proj))
            if self.ref is not None:
                self.ref_x,self.ref_y =  self.UTM(self.ref[0],self.ref[1])
            else:
                self.ref_x,self.ref_y = 0,0

    def load(self,xlim=None,ylim=None):
        fname=self.wdir+self.filename
        if self.utm_proj is None:
            x,y,z=np.loadtxt(fname,comments='#',unpack=True,dtype='f,f,f')
            # convert to meter
            self.x,self.y,self.z = x*1e3,y*1e3,z*self.scale
        else:
            self.lon,self.lat,z=np.loadtxt(fname,comments='#',unpack=True,dtype='f,f,f')
            x, y = self.UTM(self.lon, self.lat) 
            self.x,self.y,self.z=(x-self.ref_x),(y-self.ref_y),z*self.scale
        
        # remove data outside map
        if (xlim is not None) and (ylim is not None):
          index=np.nonzero((self.x<xlim[0])|(self.x>xlim[1])|(self.y<ylim[0])|(self.y>ylim[1]))
          self.x,self.y,self.z = np.delete(self.x,index),np.delete(self.y,index),np.delete(self.z,index)

class shapefile:
    """
    shapefiel class
    Parameters:
    name,filename: name input file, given name
    wdir: path input file
    edgecolor, color, linewidth: plot option
    utm_proj: EPSG UTM projection. If not None, project data from to EPSG.
    """
    
    def __init__(self,name,wdir,filename,color='black',edgecolor='black',linewidth=2.,utm_proj=None):
        self.name=name
        self.filename=filename
        self.wdir=wdir
        self.color=color
        self.edgecolor=edgecolor
        self.linewidth=linewidth
        self.crs=utm_proj

class seismicity:
    """
    seismicity class: read usgs csv files
    Parameters:
    name, filename : give name, fiel name 
    wdir: path input file
    color, width: plot option
    utm_proj: EPSG UTM projection. If not None, project data from to EPSG.
    if fmt = 'csv':
    Column attributes: time,latitude,longitude,depth,mag,magType,nst,gap,dmin,rms,net,id,updated,place,type,horizontalError,depthError,magError,magNst,status,locationSource,magSource
    if fmt = 'txt':
    Column attributes: latitude,longitude,depth,mag
    """    
    
    def __init__(self,name,wdir,filename,color='black',width=2.,utm_proj=None,ref=None,fmt='csv'):
        self.name=name
        self.filename=filename
        self.wdir=wdir
        self.color=color
        self.width=width
        self.utm_proj=utm_proj
        self.ref=ref
        self.fmt = fmt

        # projection
        if self.utm_proj is not None:
            self.UTM = pyproj.Proj("EPSG:{}".format(self.utm_proj))
            if self.ref is not None:
                self.ref_x,self.ref_y =  self.UTM(self.ref[0],self.ref[1])
            else:
                self.ref_x,self.ref_y = 0,0

    def load(self,xlim=None,ylim=None):
        fname=self.wdir+self.filename
        if self.fmt == 'csv':
          df = pandas.read_csv(fname)
          lat,lon,depth,mag=df['latitude'][:].to_numpy(),df['longitude'][:].to_numpy(),df['depth'][:].to_numpy(),df['mag'][:].to_numpy() 
          if self.utm_proj is None:
            self.x,self.y,self.z,self.mag = lon,lat,depth,mag
          else:
            x, y = self.UTM(lon, lat)
            self.x,self.y,self.z,self.mag=(x-self.ref_x),(y-self.ref_y),depth,mag
        elif  self.fmt == 'txt':
          date,self.mag,y,x,depth = np.loadtxt(fname,comments = '#',unpack = True,dtype = 'S,f,f,f,f')
          if self.utm_proj is not None:
             x, y = self.UTM(x, y)
             self.x, self.y = (x - self.ref_x), (y - self.ref_y)
          else:
             self.x, self.y = x*1e3, y*1e3 
        if np.nanmean(abs(depth)) < 100:
            self.depth = depth*1000
           



