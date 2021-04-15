import numpy as np
import math,sys

#GMT files
class gmt:
    def __init__(self,name,wdir,filename,color='black',width=2.,utm_proj=None, ref=None):
        self.name=name
        self.wdir=wdir
        self.filename=filename
        self.color=color
        self.width=width
        
        # projection
        self.utm_proj=utm_proj
        self.ref=ref
        if self.utm_proj is not None:
            import pyproj
            self.UTM = pyproj.Proj("+init=EPSG:{}".format(self.utm_proj))
            if self.ref is not None:
                self.ref_x,self.ref_y =  self.UTM(self.ref[0],self.ref[1])
            else:
                self.ref_x,self.ref_y = 0,0

        self.x=[]
        self.y=[]

        self.xp=[]
        self.yp=[]

    #load gmt segments
    def load(self,delimiter=' ',xlim=[-1000,1000],ylim=[-1000,1000]):
        x=[[]]
        y=[[]]
        i=0
        
        infile = open(self.wdir+self.filename,"r")
        for line in infile:
            if '>' in line:
                i=i+1
                x.append([])
                y.append([])
            else:
                temp = list(map(float, line.split()))
                xt, yt = temp[0],temp[1]
                if self.utm_proj is not None:
                    lon,lat = np.float(xt),np.float(yt)
                    xt, yt = self.UTM(lon, lat)
                    xt, yt = xt-self.ref_x,yt-self.ref_y
                else:
                    # convert km to m
                     xt, yt = xt*1e3, yt*1e3
                if xlim is not None:
                  if (xt>xlim[0]) and (xt<xlim[1]) and (yt>ylim[0]) and (yt<ylim[1]):
                    x[i].append(xt)
                    y[i].append(yt)
                else:
                    x[i].append(xt)
                    y[i].append(yt)

        return x,y
