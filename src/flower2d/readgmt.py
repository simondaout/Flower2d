import numpy as np
import math,sys
import csv


#GMT files
class gmt:
    """
    Class gmt: read GMT file: text file with object separated with >
    filename: name input file
    wdir: path input file
    name: name for plot
    color, width options for plot
    """
    def __init__(self,name,wdir,filename,color='black',width=1.):
        self.name=name
        self.wdir=wdir
        self.filename=filename
        self.color=color
        self.width=width

        self.x=[]
        self.y=[]

        self.xp=[]
        self.yp=[]

    #load gmt segments
    def load(self,delimiter=' ',ybounds=None,xbounds=None):
        n=0
        x=[[]]
        y=[[]]
        i=0

        infile=csv.reader(open(self.wdir+self.filename,"rb"),delimiter=delimiter)
        for line in infile:
            if line[0] is '>':
                i=i+1
                x.append([])
                y.append([])
            else:
                x[i].append(float(line[0]))
                y[i].append(float(line[1]))

        return x,y
