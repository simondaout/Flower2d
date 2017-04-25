import numpy as np
import math,sys
import csv


#GMT files
class gmt:
    def __init__(self,name,wdir,filename,color='black',width=2.):
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
    def load(self,delimiter=' '):
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
                x[i].append(line[0])
                y[i].append(line[1])

        return x,y
