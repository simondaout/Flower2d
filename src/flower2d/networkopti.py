import numpy as np
import math
import copy
import scipy.optimize as sp
import matplotlib.pyplot as plt

# from flatten import *
from modelopti import *
import time

def gaussian(x, mu, sig):
    return 1/(np.sqrt(2*np.pi)*sig)*np.exp(-np.power((x-mu), 2.) / (2 * np.power(sig, 2.)))

def expCov(t, sil, lam, sig):
    return sil - (sig**2)*np.exp(-t/lam)

class network(object):
    def __init__(self,network,reduction,wdir,dim,weight=1.,scale=1.,errorfile=None,\
        los=None,heading=None,color='black',mask=None,perc=100,cov=None):
        
        # network name
        self.network = network
        self.reduction = reduction
        self.wdir = wdir
        # data weighting
        self.wd = 1./weight 
        #self.sigmad = 1./weight
        # scaling between network
        self.scale = scale
        # number of stations
        self.Npoint = 0 
        # number of displacement components
        self.dim = dim
        self.d = []
        self.N = 0
        self.x,self.y,self.z = [],[],[]
        self.upar,self.uperp = [],[]
        # model
        self.mx,self.my = [],[]
        self.mlos = []
        self.mpar,self.mperp = [],[]
        # base
        #self.a,self.b = base[0],base[1]
        self.a,self.b,self.c = 0,0,0
        # optional file with errorbar
        self.errorfile=errorfile
        # optional los file
        self.los=los
        self.heading=heading
        self.color=color
        self.mask=mask
        self.perc=perc
        self.cov=cov

    def load(self,flt):
        f = file(self.wdir+self.network,'r')

        self.fmodel = flt.fmodel
        self.struc = flt.structures 
        self.volum = flt.volum
        self.Mseg = flt.Mseg 
        self.Mstruc = flt.Mstruc
        self.Mvol = flt.Mvol
        self.Mdis = flt.Mdis
        self.str = flt.str
        self.profile = flt.profiles
        #self.proj =  flt.profiles.proj

        if 1 == self.dim:     # InSAR case
          if self.los is not None:
            x,y,los,theta=np.loadtxt(f,comments='#',unpack=True,dtype='f,f,f,f')
            xp=(x-self.fmodel[0].x)*self.profile.s[0]+(y-self.fmodel[0].y)*self.profile.s[1]
            yp=(x-self.fmodel[0].x)*self.profile.n[0]+(y-self.fmodel[0].y)*self.profile.n[1]
            index=np.nonzero((xp>self.profile.xpmax)|(xp<self.profile.xpmin)|(yp>self.profile.ypmax)|(yp<self.profile.ypmin))
            self.ulos,self.x,self.y,self.xp,self.yp,self.theta=np.delete(los,index),np.delete(x,index),np.delete(y,index),np.delete(xp,index),np.delete(yp,index),np.delete(theta,index)
            self.phi = np.deg2rad(-90-self.heading)
            self.theta = np.deg2rad(90.-self.theta)
            phim,thetam=np.mean(self.phi),np.mean(self.theta)
            print np.rad2deg(phim),np.rad2deg(thetam)
            self.proj=[np.cos(self.theta)*np.cos(self.phi),
            np.cos(self.theta)*np.sin(self.phi),
            np.sin(self.theta)
            ]
            self.projm=[np.cos(thetam)*np.cos(phim),
            np.cos(thetam)*np.sin(phim),
            np.sin(thetam)]
            print self.projm
            flt.profiles.proj = self.projm
            # sys.exit()
          else:
            x,y,los=np.loadtxt(f,comments='#',unpack=True,dtype='f,f,f')
            xp=(x-self.fmodel[0].x)*self.profile.s[0]+(y-self.fmodel[0].y)*self.profile.s[1]
            yp=(x-self.fmodel[0].x)*self.profile.n[0]+(y-self.fmodel[0].y)*self.profile.n[1]
            index=np.nonzero((xp>self.profile.xpmax)|(xp<self.profile.xpmin)|(yp>self.profile.ypmax)|(yp<self.profile.ypmin))
            self.ulos,self.x,self.y,self.xp,self.yp=np.delete(los,index),np.delete(x,index),np.delete(y,index),np.delete(xp,index),np.delete(yp,index)
            if self.profile.proj is not None:
                self.proj= self.profile.proj
            else:
                print 'No average projection look angle set in the profile class'
                sys.exit(2)

          if self.mask is not None:
            uu = np.flatnonzero(np.logical_or(self.yp<=self.mask[0], self.yp>=self.mask[1]))
            self.yp = self.yp[uu]
            self.xp = self.xp[uu]
            self.x = self.x[uu]
            self.y = self.y[uu]
            self.ulos = self.ulos[uu]

          # optional cleaning
          bins = np.arange(self.profile.ypmin,self.profile.ypmax,2)
          inds = np.digitize(self.yp,bins)
          # print inds
          distance = []
          ypt = []
          xpt = []
          xt = []
          yt = []
          ulost = []

          # print len(self.ulos)

          for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            if len(uu)>0:
                distance.append(bins[j] + (bins[j+1] - bins[j])/2.)
                indice = np.flatnonzero(np.logical_and(self.ulos[uu]>=np.percentile(\
                  self.ulos[uu],100-self.perc),self.ulos[uu]<=np.percentile(self.ulos[uu],self.perc)))
                ypt.append(self.yp[uu][indice])
                xpt.append(self.xp[uu][indice])
                yt.append(self.y[uu][indice])
                xt.append(self.x[uu][indice])
                ulost.append(self.ulos[uu][indice])

          self.ulos = np.concatenate(np.array(ulost))
          # print len(self.ulos)
          # sys.exit()
          self.yp = np.concatenate(np.array(ypt))
          self.xp = np.concatenate(np.array(xpt))
          self.x = np.concatenate(np.array(xt))
          self.y = np.concatenate(np.array(yt))

          # optional scaling 
          self.ulos = self.ulos*self.scale
          # Remove points that are at the same location (to avoid the covariance matrix being non definite positive)
          self.yp, uu = np.unique(self.yp, return_index = True)
          self.xp = self.xp[uu]
          self.x = self.x[uu]
          self.y = self.y[uu]
          self.ulos = self.ulos[uu]


          if self.los is not None:
            self.proj[0],self.proj[1],self.proj[2] = \
            self.proj[0][uu],self.proj[1][uu],self.proj[2][uu]

          self.Npoint,self.N = len(self.ulos),len(self.ulos)
	
          # number of data
          self.N = self.Npoint*self.dim
          # data vector
          self.d = np.atleast_1d(self.ulos)
          self.sigmad = np.ones((self.Npoint))*self.wd

        elif 2 == self.dim: # GPS case
            name,x,y = np.loadtxt(f,comments = '#',unpack = True,dtype = 'S4,f,f')
            # perp and par composante ref to the first fault 
            yp = (x-self.fmodel[0].x)*self.profile.n[0]+(y-self.fmodel[0].y)*self.profile.n[1]
            xp = (x-self.fmodel[0].x)*self.profile.s[0]+(y-self.fmodel[0].y)*self.profile.s[1]
            # data selection
            index = np.nonzero((xp>self.profile.xpmax)|(xp<self.profile.xpmin)|(yp>self.profile.ypmax)|(yp<self.profile.ypmin))

            self.name,self.x,self.y,self.xp,self.yp = np.delete(name,index),np.delete(x,index),np.delete(y,index),np.delete(xp,index),np.delete(yp,index)
            self.Npoint = len(self.name)
            
            # proj vor gps vectors
            self.proj= self.profile.proj

            self.ux,self.uy = np.zeros(self.Npoint),np.zeros(self.Npoint)
            self.upar,self.uperp = np.zeros(self.Npoint),np.zeros(self.Npoint)
            self.sigmax,self.sigmay = np.ones(self.Npoint),np.ones(self.Npoint)
            self.sigmapar,self.sigmaperp = np.zeros(self.Npoint),np.zeros(self.Npoint)
            self.ulos,self.sigmalos = np.zeros(self.Npoint),np.zeros(self.Npoint)
            
            for i in xrange(self.Npoint):
                station = self.wdir+self.reduction+'/'+self.name[i]
                dated,east,north,esigma,nsigma = np.loadtxt(station,comments = '#',usecols = (0,1,2,3,4),unpack = True,dtype = 'f,f,f,f,f')
                self.ux[i],self.uy[i] = east*self.scale,north*self.scale
                self.sigmax[i],self.sigmay[i] = esigma*self.scale*self.wd,nsigma*self.scale*self.wd
                
                # Extract perpendicular en parallele components for GPS data
                self.upar[i] = self.ux[i]*self.profile.s[0]+self.uy[i]*self.profile.s[1]
                self.uperp[i] = self.ux[i]*self.profile.n[0]+self.uy[i]*self.profile.n[1]
                # attention faux si erreur ellipsoidale
                #self.sigmapar[i]=self.sigmax[i]
                #self.sigmaperp[i]=self.sigmay[i]
                self.sigmaperp[i]=((self.sigmax[i]*np.cos(self.str))**2 + (self.sigmay[i]*np.sin(self.str))**2)**0.5
                self.sigmapar[i]=((self.sigmax[i]*np.sin(self.str))**2 + (self.sigmay[i]*np.cos(self.str))**2)**0.5
                
                self.ulos[i] = self.ux[i]*self.proj[0]+self.uy[i]*self.proj[1]
                self.sigmalos[i] = self.sigmax[i]*self.proj[0]+self.sigmay[i]*self.proj[1]

            #self.sig = np.column_stack([self.sigmapar,self.sigmaperp])

            if self.errorfile is not None:
                errorfile=self.errorfile
                east,north,self.sigmapar,self.sigmaperp=np.loadtxt(errorfile,comments = '#',usecols = (0,1,2,3),unpack = True,dtype = 'f,f,f,f')
                self.sigmapar,self.sigmaperp=2*abs(self.sigmapar)*self.scale*self.wd,2*abs(self.sigmaperp)*self.wd*self.scale

            self.u = np.column_stack([self.upar,self.uperp])
            uncertainties = np.column_stack([self.sigmax,self.sigmay,self.sigmapar,self.sigmaperp])
            self.sig = np.column_stack([self.sigmapar,self.sigmaperp])
            print 'Uncertainties:'
            print 'East_Dev_(mm/yr)    North_Dev_(mm/yr) Par_Dev(mm/yr) Perp_Dev(mm/yr)'
            print uncertainties
            
            # number of data
            self.N = self.Npoint*self.dim
            # data vector
            self.d = np.zeros((self.N))
            # uncertainty vector
            self.sigmad = np.zeros((self.N))
            for i in xrange(self.Npoint):
                for k in xrange(self.dim):
                    self.d[self.dim*i+k] = self.u[i,k]
                    self.sigmad[self.dim*i+k] = self.sig[i,k]

        else :
            name,x,y = np.loadtxt(f,comments = '#',unpack = True,dtype = 'S4,f,f')
            # perp and par composante ref to the first fault 
            yp = (x-self.fmodel[0].x)*self.profile.n[0]+(y-self.fmodel[0].y)*self.profile.n[1]
            xp = (x-self.fmodel[0].x)*self.profile.s[0]+(y-self.fmodel[0].y)*self.profile.s[1]
            # data selection
            index = np.nonzero((xp>self.profile.xpmax)|(xp<self.profile.xpmin)|(yp>self.profile.ypmax)|(yp<self.profile.ypmin))

            self.name,self.x,self.y,self.xp,self.yp = np.delete(name,index),np.delete(x,index),np.delete(y,index),np.delete(xp,index),np.delete(yp,index)
            self.Npoint = len(self.name)
            
            # proj vor gps vectors
            self.proj= self.profile.proj
            
            self.ux,self.uy,self.uv = np.zeros(self.Npoint),np.zeros(self.Npoint),np.zeros(self.Npoint)
            self.upar,self.uperp = np.zeros(self.Npoint),np.zeros(self.Npoint)
            self.sigmax,self.sigmay,self.sigmav = np.ones(self.Npoint),np.ones(self.Npoint),np.zeros(self.Npoint)
            self.sigmapar,self.sigmaperp = np.zeros(self.Npoint),np.zeros(self.Npoint)
            self.ulos,self.sigmalos = np.zeros(self.Npoint),np.zeros(self.Npoint)
            for i in xrange(self.Npoint):
                station = self.wdir+self.reduction+'/'+self.name[i]
                dated,east,north,up,esigma,nsigma,upsigma = np.loadtxt(station,comments = '#',usecols = (0,1,2,3,4,5,6),unpack = True,dtype = 'f,f,f,f,f,f,f')
                self.ux[i],self.uy[i],self.uv[i] = east*self.scale,north*self.scale,up*self.scale
                self.sigmax[i],self.sigmay[i],self.sigmav[i] = esigma*self.scale*self.wd,nsigma*self.scale*self.wd,upsigma*self.scale*self.wd
                
                # Extract perpendicular en parallele components for GPS data
                self.upar[i] = self.ux[i]*self.profile.s[0]+self.uy[i]*self.profile.s[1]
                self.uperp[i] = self.ux[i]*self.profile.n[0]+self.uy[i]*self.profile.n[1]
                # attention faux si erreur ellipsoidale
                #self.sigmapar[i]=self.sigmax[i]
                #self.sigmaperp[i]=self.sigmay[i]
                self.sigmaperp[i]=((self.sigmax[i]*np.cos(self.str))**2 + (self.sigmay[i]*np.sin(self.str))**2)**0.5
                self.sigmapar[i]=((self.sigmax[i]*np.sin(self.str))**2 + (self.sigmay[i]*np.cos(self.str))**2)**0.5

                self.ulos[i] = self.ux[i]*self.proj[0]+self.uy[i]*self.proj[1]+self.uv[i]*self.proj[2]
                self.sigmalos[i] = self.sigmax[i]*self.proj[0]+self.sigmay[i]*self.proj[1]+self.sigmav[i]*self.proj[2]
           
            print 'Load error file for gps..'
            print
            if self.errorfile is not None:
                errorfile=self.errorfile
                east,north,self.sigmapar,self.sigmaperp,self.sigmav=np.loadtxt(errorfile,comments = '#',usecols = (0,1,2,3,4),unpack = True,dtype = 'f,f,f,f,f')
                self.sigmapar,self.sigmaperp,self.sigmav=2*abs(self.sigmapar)*self.wd*self.scale,2*abs(self.sigmaperp)*self.wd*self.scale,2*abs(self.sigmav)*self.wd*self.scale

            self.u = np.column_stack([self.upar,self.uperp,self.uv])
            uncertainties = np.column_stack([self.sigmax,self.sigmay,self.sigmapar,self.sigmaperp,self.sigmav])
            self.sig = np.column_stack([self.sigmapar,self.sigmaperp,self.sigmav])
            print 'Uncertainties:'
            print 'East_Dev_(mm/yr)    North_Dev_(mm/yr) Par_Dev(mm/yr) Perp_Dev(mm/yr) Up_Dev_(mm/yr)'
            print uncertainties

            # number of data
            self.N = self.Npoint*self.dim
            # data vector
            self.d = np.zeros((self.N))
            # uncertainty vector
            self.sigmad = np.zeros((self.N))
            for i in xrange(self.Npoint):
                for k in xrange(self.dim):
                    self.d[self.dim*i+k] = self.u[i,k]
                    self.sigmad[self.dim*i+k] = self.sig[i,k]
        # All done
        return
    
    
    def computeField(self):
        # compute velocity vield
        print
        if self.dim is 3:
            disp = np.column_stack([self.yp,self.uperp,self.upar,self.uv])
            sortd = disp[disp[:,0].argsort()]
            name = self.name[disp[:,0].argsort()] 
            np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            conv = np.column_stack([sortd,sortd[:,1]-sortd[-1,1],sortd[:,2]-sortd[-1,2],sortd[:,3]])
            print ' #Name   #Distance   #PerpendicularV  #ParallelV   #Shortening   #Shearing   '     
            for i in xrange(self.Npoint):
                #print self.name[i], conv[i,:]
                print name[i], conv[i,:]
            print
        else:
            disp = np.column_stack([self.yp,self.uperp,self.upar])
            sortd = disp[disp[:,0].argsort()]
            name = self.name[disp[:,0].argsort()] 
            np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            conv = np.column_stack([sortd,sortd[:,1]-sortd[-1,1],sortd[:,2]-sortd[-1,2]])
            print ' #Name   #Distance   #PerpendicularV  #ParallelV   #Shortening   #Shearing   '     
            for i in xrange(self.Npoint):
                #print self.name[i], conv[i,:]
                print name[i], conv[i,:]
            print

        
    def computeSVD(self):
        deplac = np.zeros((self.dim,self.Npoint))
        for i in xrange(self.dim):
            deplac[i,:] = self.u[:,i]
        print

        print 'Compute SVD of the Velocity field:'
        m = np.mean(deplac,axis=1)
        #print m
        for i in xrange(self.dim):
            deplac[i,:] = deplac[i,:] - m[i] 
        U,eignv,Vt = np.linalg.svd(deplac,full_matrices=False)
        print 'Amplitude: ',eignv
        #print 'Base:', U
        
        print
        print 'Compute the velocity field in the fault base:'
        F = np.vstack((self.profile.s[:2],self.profile.n[:2]))
        proj = np.dot(np.dot(np.linalg.inv(F),U[:2,:2]),eignv[:2])
        print 'amplitude of displacement in the fault base: ', proj
    

    def computeFullCovariance(self, distmax=70., every=1., ramp='lin',maskcov=None, xbounds=[None, None], plot=False, outdir=None):
        '''
        Computes the covariance matrix for the dataset:
        Args:
            * full      : if True, estimates a full covariance matrix by computing the empirical covariance of the set (only makes sense for the InSAR).
                          if False, returns a diagonal matrix with sigmad on the diagonal
        '''

        if self.cov is None:
            sill,sigm,lamb = None, None, None
        else:
            sill,sigm,lamb = self.cov[0], self.cov[1], self.cov[2]
        if sigm is None or lamb is None:
            # get data
            if xbounds[0] is None:
                xmin = self.yp.min()
                xmax = self.yp.max()
            else:
                xmin = xbounds[0]
                xmax = xbounds[1]
        
            if maskcov is not None:
                kk = np.flatnonzero(np.logical_or(np.logical_or(np.logical_and(self.yp>=xmin, self.yp<=maskcov[0]),np.logical_and(self.yp>=maskcov[1], self.yp<=maskcov[2])),np.logical_and(self.yp>=maskcov[3], self.yp<=xmax)))
            else:
                kk = np.arange(len(self.xp))

            x = self.xp[kk]
            y = self.yp[kk]
            z = copy.deepcopy(self.d[kk])
            Nsamp = y.shape[0]
            # Remove a ramp
            kk = np.nonzero((y<1000))
            x_temp = x[kk]
            y_temp = y[kk]
            z_temp = z[kk]
        
            if ramp is 'cub':
                ## ramp ay**3 + b*y**2 + cy + dx**2 + ex + f 
                G = np.zeros((len(y_temp),6))
                G[:,0] = y_temp**3
                G[:,1] = y_temp**2
                G[:,2] = y_temp
                G[:,3] = x_temp**2
                G[:,4] = x_temp
                G[:,5] = 1.
                pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),z_temp)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]; f = pars[5] 
                z = z - (a*y**3 + b*y**2 + c*y + d*x**2 + e*x + f )
                print 'Remove ramp: {}*y**3 + {}*y**2 + {}*y + {}*x**2 + {}*x + {})'.format(a,b,c,d,e,f)
            if ramp is 'quad':
                ## ramp ay**2 + by + cx**2 + dx + e 
                G = np.zeros((len(y_temp),5))
                G[:,0] = y_temp**2
                G[:,1] = y_temp
                G[:,2] = x_temp**2
                G[:,3] = x_temp
                G[:,4] = 1.
                pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),z_temp)
                a = pars[0]; b = pars[1]; c = pars[2]; d = pars[3]; e = pars[4]  
                z = z - (a*y**2 + b*y + c*x**2 + d*x + e )
                print 'Remove ramp: {}*y**2 + {}*y + {}*x**2 + {}*x + {})'.format(a,b,c,d,e)
       
            else:
                # ramp ay + cx + d 
                G = np.zeros((len(y_temp),3))
                G[:,0] = y_temp
                G[:,1] = x_temp
                G[:,2] = 1.
                pars = np.dot(np.dot(np.linalg.inv(np.dot(G.T,G)),G.T),z_temp)
                a = pars[0]; b = pars[1]; c = pars[2];   
                z = z - (a*y + b*x + c)
                print 'Remove ramp: {}*y + {}*x + {})'.format(a,b,c)
       
            # plot profile
            fig = plt.figure(0)
            ax = fig.add_subplot(111)
            ax.scatter(y,z,s=0.1,marker='o',color='blue')
            
            # plot histo
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            histo = ax.hist(z,alpha=0.5,bins=50, normed=True, histtype='stepfilled',range=(-3,3))
            distribution = histo[0]
            bins = histo[1]
            dist = []
            for i in range(len(bins)-1):
                dist.append(bins[i] + (bins[i+1] - bins[i])/2.)
            try:
                pars, cova = sp.curve_fit(gaussian, dist, distribution)
                mu = pars[0]
                sig = pars[1]
                #a = pars[1]
                #print('amplitude: {}'.format(a))
                #print('mu: {}'.format(mu))
                #print('sig: {}'.format(sig))
                plt.axvline(x=2*sig, c="red", linestyle='--')
                plt.axvline(x=-2*sig, c="red", linestyle='--')
                exmin,exmax = ax.get_xlim()
                ex = np.arange(exmin,exmax,(exmax-exmin)/100)
                plt.plot(ex, gaussian(ex,mu,sig), '-r')
                plt.xlabel("mu: {:0.3f} \n2-sigma: {:0.3f} (mm/yr)".format(mu,sig))
                fig.savefig(outdir+'/insar/'+self.reduction+'_dist.eps', format='EPS')
            except:
                print('No solution found for gaussian estimation')

            #sys.exit()
            plt.show()
            
            # compute variance
            var = np.std(z)**2
            print
            print 'variance: {} mm/yr**2'.format(var)
            print
            
            # Build the permutations
            ii, jj = np.meshgrid(range(Nsamp), range(Nsamp))
            ii = ii.flatten()
            jj = jj.flatten()
            uu = np.flatnonzero(ii>jj)
            ii = ii[uu]
            jj = jj[uu]
            # Compute the distances
            #dis = np.abs(y[ii] - y[jj])
            dis = np.sqrt(np.power((y[ii] - y[jj]),2) + np.power((x[ii] - x[jj]),2))
            # Compute the semivariogram
            dr = z[ii] - z[jj]
            dv = (dr)**2
            # Digitize
            bins = np.arange(0., distmax, every)
            inds = np.digitize(dis, bins)
            # Average
            distance = []
            semivariogram = []
            sigma = []
            for i in range(len(bins)-1):
                uu = np.flatnonzero(inds==i)
                if len(uu)>0:
                    distance.append(bins[i] + (bins[i+1] - bins[i])/2.)
                    semivariogram.append(np.mean(dv[uu])/2)
                    #sigma.append(np.std(dv[uu]))
            distance = np.array(distance)
            #sigma = np.array(sigma)
            semivariogram = np.array(semivariogram)
            # Fit an exponential
            try:
                pars, cova = sp.curve_fit(expCov, distance, semivariogram)                  
            except:
                try:                                                                  
                    pars, cova = sp.curve_fit(expCov, distance, semivariogram, ftol=1e-5)   
                except:
                    self.semivariogram = {'distance':distance,
                                      'semivariogram':semivariogram}
                    print('No solution found for covariance estimation')
                    sys.exit(1)
            sill = pars[0]          
            lamb = pars[1]
            sigm = pars[2]
            # sigma**2=sill
            #sigm=np.sqrt(sill) 
            # Show me
            if plot:
                print('Sill: {}'.format(sill))
                print('Sigma: {} mm/yr'.format(sigm))
                print('Lambda: {} km'.format(lamb))
                
                fig = plt.figure(2,figsize = (10,6))
                ax = fig.add_subplot(122)
                ## Covario
                ax.plot(distance, sill - semivariogram, '.k')
                ax.plot(distance, sill - expCov(distance, sill, lamb,sigm), '-r')
                plt.xlabel('Distance between observations (km)')
                plt.ylabel('Covariance (mm/yr)**2')
               
                ax = fig.add_subplot(121)
                ## Semivario
                ax.plot(distance, semivariogram, '.k')
                ax.plot(distance, expCov(distance, sill, lamb,sigm), '-r')
                #ax.errorbar(distance, sill - expCov(distance, sill, lamb),yerr=sigma,fmt=None)
                ymin,ymax=plt.ylim()
                plt.ylim([0,ymax])
                plt.xlabel('Distance between observations (km)')
                plt.ylabel('Semivariogram (mm/yr)**2')
                fig.savefig(outdir+'/insar/'+self.reduction+'_semivar.eps', format='EPS')
                #time.sleep(10)
                plt.show()
            
        # Build the covariance matrix
        Cd = np.zeros((self.xp.shape[0], self.xp.shape[0]))
        Xh, Xv = np.meshgrid(self.xp, self.xp.T)
        Yh, Yv = np.meshgrid(self.yp, self.yp.T)
        Cd = sigm**2 * np.exp(-np.sqrt((Xh-Xv)**2+(Yh-Yv)**2)/lamb)
        if sill > sigm**2:
            # fill the diagonal with the variance value
            np.fill_diagonal(Cd,sill)
        return Cd
        #import pymc
        #Cd = sigm**2 * np.exp(-np.abs(Xh-Xv)/lamb)
        #Cd = pymc.gp.Covariance(eval_fun=pymc.gp.cov_funs.exponential.euclidean, amp=sigm, scale=lamb)
        # All done
        #return Cd(self.yp, self.yp)
   
    def residual(self,m,det,invCd):
        g=np.asarray(self.g(m))
        return np.sqrt(1/(det*2*np.pi))*np.exp(-0.5*invCd*(self.d-g)**2)
    
    def g(self,m):
        # initialize g matrix
        g = np.zeros((self.N))

        self.fmodel[0].ss = m[0]
        # Compute SS on the Shear zone 
        tot_ss = 0
        for l in range(1,self.Mseg):
                tot_ss += m[l*3]
        self.fmodel[0].sst = m[0] + tot_ss

        # Get model parameters for the main fault
        self.fmodel[0].ds, self.fmodel[0].H = m[1:3]
       
        if self.struc[0].Mseg == 1:
            self.struc[0].conservation()
            u = self.fmodel[0].displacement(self.yp)
        else:
            # Get model parameters for the back-thrust
            self.fmodel[1].ss, self.fmodel[1].D, self.fmodel[1].H = m[3], m[4], m[5]
            # Get model parameters for the ramp 
            self.fmodel[2].ss,self.fmodel[2].D, self.fmodel[2].H = m[6], m[7], m[8]
            # conservation of motion
            self.struc[0].conservation()
            #self.fmodel[1].info()
            #self.fmodel[2].info()
            u = self.fmodel[0].displacement(self.yp)

            # control on locking depth: depth cannot be negatif
            if (self.fmodel[1].w < 0) or (self.fmodel[2].w < 0) or (abs(self.fmodel[2].vh) > abs(m[1])):
                return np.ones((self.N,))*1e14
        
            # control on lenght
            #if (self.fmodel[1].L < 5) or (self.fmodel[2].L < 5):
            #    return np.ones((self.N,))*1e14

            #print m[1], self.fmodel[2].ds, self.fmodel[2].ds*math.cos(self.fmodel[2].dipr)
            u = u + self.fmodel[1].displacement(self.yp)
            u = u + self.fmodel[2].displacement(self.yp)
        
        # Index
        start = self.struc[0].Mker
        Mtemp = self.struc[0].Mseg 
        # Iterate over the segments
        for j in xrange(1,self.Mstruc):
            # loop on all seg of the struc
            for k in xrange(self.struc[j].Mseg):
                ramp = self.struc[j].segments[k]
                Mker = ramp.Mker 
                ramp.ss,ramp.D,ramp.H = m[start:start+Mker]
                start += Mker

            # conservation of motion
            self.struc[j].conservation(self.fmodel[Mtemp-1],self.fmodel[0])
           
            # compute dispalcements
            for k in xrange(self.struc[j].Mseg):
                #print
                #print ramp.info()
                ramp = self.struc[j].segments[k]
                
                # control on locking depths: depth cannot be negatif
                if (ramp.w < 0) : # !!!! put after conservation
	            return np.ones((self.N,))*1e14
            
                #control on lenght
                #if (ramp.L < 5):
                #    return np.ones((self.N,))*1e14

                u = u + self.fmodel[Mtemp+k].displacement(self.yp)

            # control on conservation of motion
            #if abs(self.struc[j].segments[0].vh) > abs(m[1]) or (np.sign(self.struc[j].segments[0].ds) != np.sign(m[1])):
            #    return np.ones((self.N,))*1e14

            Mtemp += self.struc[j].Mseg
            #start = start + self.struc[j].Mker 
            
        # loop on volumic deformation structures
        for j in xrange(self.Mvol):
            self.volum[j].ds,self.volum[j].D = m[self.Mdis+j],m[self.Mdis+j+1]
            #print m[self.Mdis+j] 
            u = u + self.volum[j].displacement(self.yp)

        if 1 == self.dim:
            self.mlos = u[:,0]*self.profile.s[0]*self.proj[0]+\
              u[:,0]*self.profile.s[1]*self.proj[1]+\
              u[:,1]*self.profile.n[0]*self.proj[0]+\
              u[:,1]*self.profile.n[1]*self.proj[1]+\
              u[:,2]*self.proj[2]
            self.dp = self.mlos + self.a+self.b*self.yp # baseline term
        
        elif 2 == self.dim:
            base = np.ones((u.shape[0], 2))
            base[:,0] *=  self.a
            base[:,1] *=  self.b
            self.dp = (u[:,:2] + base).flatten()

            # mpar/mperp without baseline term
            self.mpar = u[:,0]
            self.mperp = u[:,1]
            self.mx=u[:,0]*self.profile.s[0]+u[:,1]*self.profile.n[0]
            self.my=u[:,0]*self.profile.s[1]+u[:,1]*self.profile.n[1]
            self.mz=u[:,2]
            self.mlos=self.mx*self.proj[0]+self.my*self.proj[1]+self.mz*self.proj[2]
        
        else:
            base = np.ones((u.shape[0], 3))
            base[:,0] *=  self.a
            base[:,1] *=  self.b
            base[:,2] *=  self.c
            self.dp = (u[:,:3] + base).flatten()

            self.mpar = u[:,0]
            self.mperp = u[:,1]
            self.mx=u[:,0]*self.profile.s[0]+u[:,1]*self.profile.n[0]
            self.my=u[:,0]*self.profile.s[1]+u[:,1]*self.profile.n[1]
            self.mz=u[:,2]
            self.mlos=self.mx*self.proj[0]+self.my*self.proj[1]+self.mz*self.proj[2]

        return self.dp

    def info(self):
        print "%s Data: " %self.reduction, self.d

