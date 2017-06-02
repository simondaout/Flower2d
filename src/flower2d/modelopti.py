import numpy as np
import math
import sys
# from angles import normalize 

class profile:
    def __init__(self,name,x,y,l,w,proj=None):
        # profile parameters
        self.name = name
        self.x = x
        self.y = y
        self.l = l
        self.w = w
        self.proj = proj

class inversion:
    def __init__(self,name,structures,strike,profiles):
        self.name = name
        self.fmodel = []
        self.structures = structures
        self.strike = strike
        self.profiles = profiles
        
        d2r =  math.pi/180
        self.str = strike*d2r
        self.s = [math.sin(self.str),math.cos(self.str),0]
        self.n = [math.cos(self.str),-math.sin(self.str),0]
        # strike of the fault is the strike of the profile
        self.profiles.s = self.s
        self.profiles.n = self.n

class segment:
    def __init__(self,name,ss,sigmass,H,sigmaH,distss,distH):
        self.name = name
        # ss: strike-slip
        self.ss = ss
        self.sigmass = sigmass
        # H: vertical distance
        self.H = H
        self.sigmaH = sigmaH
        # w: locking depth
        self.w, self.sigmaw = 0,0
        # horizontal distance to the main fault
        self.fperp = 0
        self.distss=distss
        self.distshort=distss
        self.distH=distH

    def displacement(self,yp):
        
        #######################################################
        # Surface Displacements Due to Edge Dislocation
        # Paul Segall 1996
        # Input:
        #       yp = coordinates of the observation points
        #       w1 = depth of updip end;
        #       L = length downdip;
        #######################################################

        L = self.L
        dipr = self.dipr
        # D: distante to the main fault
        shift = -self.fperp*np.ones(yp.shape)
        z = 0
        w1 = self.w
        w2 = self.w+L*math.sin(dipr)
        zeta1 = (yp+shift)/w1
        zeta2 = (yp+shift-L*math.cos(dipr))/w2
        denom1 = 1+zeta1**2
        denom2 = 1+zeta2**2
        u = np.zeros((yp.size,3))

        # vertical displacement
        uv =  (1./math.pi) * (math.sin(dipr)*np.arctan(zeta1) + (math.cos(dipr)+math.sin(dipr)*zeta1)/denom1 - math.sin(dipr)*np.arctan(zeta2) - (math.cos(dipr)+math.sin(dipr)*zeta2)/denom2 )
    
        # horizontal displacement
        uperp =  -(1./math.pi) * (math.cos(dipr)*np.arctan(zeta1) + (math.sin(dipr) -math.cos(dipr)*zeta1)/denom1 - math.cos(dipr)*np.arctan(zeta2) - (math.sin(dipr) - math.cos(dipr)*zeta2)/denom2)

        #######################################################
        # Surface Displacements Due to an Half-infinite dipping dislocation
        # Paul Segall 1996
        #######################################################

        # half-infinite dislocation
        upar =  (1./(2*math.pi)) * (np.arctan2(math.sin(dipr)*(yp+shift)+math.cos(dipr)*(self.w+z),math.sin(dipr)*(self.w+z)-math.cos(dipr)*(yp+shift)) + \
            np.arctan2(math.sin(dipr)*(yp+shift)-math.cos(dipr)*(z-self.w),-math.cos(dipr)*(yp+shift)-math.sin(dipr)*(z-self.w))) 
        
        # if not creeping segment, half-infinite dislocation for strike-slip
        #if abs(self.D) > 10:
            #L = 660
            #print self.name,self.L,self.fperp

        # Linf = 660
        # upar = upar -  (1./(2*math.pi)) * (np.arctan2(math.sin(dipr)*(yp+shift)+math.cos(dipr)*(self.w+z),math.sin(dipr)*(self.w+z)-math.cos(dipr)*(yp+shift)) + \
        #     np.arctan2(math.sin(dipr)*(yp+shift)-math.cos(dipr)*(z-self.w),-math.cos(dipr)*(yp+shift)-math.sin(dipr)*(z-self.w))) - \
        #     ( (1./(2*math.pi)) * (np.arctan2(math.sin(dipr)*(yp+shift+Linf*math.cos(dipr))+math.cos(dipr)*(self.w+z+Linf*math.sin(dipr)),math.sin(dipr)*(self.w+z+Linf*math.sin(dipr))-math.cos(dipr)*(yp+shift+Linf*math.cos(dipr))) +
        #     np.arctan2(math.sin(dipr)*(yp+shift+Linf*math.cos(dipr))-math.cos(dipr)*(z-(self.w+Linf*math.sin(dipr))),-math.cos(dipr)*(yp+shift+Linf*math.cos(dipr))-math.sin(dipr)*(z-(self.w+Linf*math.sin(dipr))))) )

        u[:,0],u[:,1],u[:,2] = self.ss*upar,self.ds*uperp,self.ds*uv

        return u

# one decollement as first structure 
class maindecol:
    def __init__(self,name,x,y,ss,sigmass,short,sigmashort,w,sigmaw,dip,\
        distss='Unif',distshort='Unif',distH='Unif'):
       
        self.segments = [
                main(name = name,x = x,y = y,ss = ss,sigmass = sigmass,ds = short,sigmads = sigmashort, H = w,sigmaH = sigmaw,dip = dip,distss=distss,distshort=distshort,distH=distH),
                ]
        self.Mseg = 1
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))
        self.segments[0].w = w
        self.segments[0].sigmaw = sigmaw
        
    def conservation(self):
        self.segments[0].w = self.segments[0].H
        # print  self.segments[0].H, self.segments[0].w
        self.segments[0].vh = self.segments[0].ds
        self.segments[0].gamma = 0

# flower structure as first structure
class mainflower:
    def __init__(self,name,x,y,ss,sigmass,short,sigmashort,w,sigmaw,dip,\
        name2,H2,sigmaH2,ss2,sigmass2,D2,sigmaD2,\
        name3,H3,sigmaH3,ss3,sigmass3,D3,sigmaD3,\
        distss='Unif',distshort='Unif',distH='Unif',distD='Unif'):
       
        self.segments = [
                main(name = name,x = x,y = y,ss = ss,sigmass = sigmass,ds = short,sigmads = sigmashort,H = w,sigmaH = sigmaw,dip = dip, distss=distss, distshort=distshort, distH=distH),
                second(name = name3,ss = ss3,sigmass = sigmass3,D = D3,sigmaD = sigmaD3,H = H3,sigmaH = sigmaH3, distss=distss, distH=distH, distD=distD), # kink
                second(name = name2,ss = ss2,sigmass = sigmass2,D = D2,sigmaD = sigmaD2,H = H2,sigmaH = sigmaH2, distss=distss, distH=distH, distD=distD), # ramp
                ]
        self.Mseg = len(self.segments)
        self.segments[0].w,self.segments[1].w,self.segments[2].w = self.segments[0].H, self.segments[0].H - self.segments[1].w, self.segments[0].H - self.segments[2].w
        self.segments[0].sigmaw = sigmaw
        self.segments[1].fperp,self.segments[2].fperp =  self.segments[1].D, self.segments[2].D  
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))

    def conservation(self):
        # update LD !!
        self.segments[0].w = self.segments[0].H
        self.quadrant = self.segments[0].quadrant
        # print 
        # print  self.quadrant
        # Get model parameters for the main fault
        v1, H1 = self.segments[0].ds, self.segments[0].H
        # print
        # print 'v1:', v1
        # print 'dip decol:', self.segments[0].dip

        self.segments[1].fperp =  self.segments[1].D 
        self.segments[1].w = self.segments[0].H - self.segments[1].H
        self.segments[1].L = math.sqrt(self.segments[1].D**2+(self.segments[1].H)**2)
        # True Dip
        self.segments[1].dipr = math.atan2((self.segments[1].H),-self.segments[1].D)
        self.segments[1].dip = self.segments[1].dipr*180/math.pi
        # si dip0 = 0, D>0, H>0 --> -D pour que beta entre 90 et 180
        # si dip0 = 180, D<0, H>0 --> D pour que beta entre 90 et 180
        self.beta = math.atan2((self.segments[1].H),-self.segments[1].D*self.quadrant)
        # print 'D,H , beta: ',  self.segments[1].D, self.segments[1].H, np.rad2deg(self.beta)

        self.segments[2].fperp = self.segments[2].D
        self.segments[2].w = self.segments[0].H - self.segments[2].H
        self.segments[2].L = math.sqrt(self.segments[2].D**2+(self.segments[2].H)**2)
        # si dip0 = 0, D<0, H>0 --> -D pour que gamma entre 0 et 90
        # si dip0 = 180, D>0, H>0 --> D pour que beta entre 0 et 90
        self.segments[2].dipr = math.atan2((self.segments[2].H),-self.segments[2].D)
        self.segments[2].dip = self.segments[2].dipr*180/math.pi
        # gamma betzeen 0 and pi/2
        # save gamma in segment attribute for next structure
        self.segments[2].gamma = math.atan2((self.segments[2].H),-self.segments[2].D*self.quadrant)
        self.gamma = self.segments[2].gamma
        # print 'D, H, gamma: ',  self.segments[2].D, self.segments[2].H, np.rad2deg(self.gamma)

        # dip0=0 --> alpha = gamma
        # kink
        self.segments[1].ds = -v1*math.sin(self.gamma)/math.sin(self.beta-self.gamma)
        # ramp
        self.segments[2].ds = v1*math.sin(self.beta)/math.sin(self.beta-self.gamma)
        # print 'v3:', self.segments[1].ds 
        # print 'v2:', self.segments[2].ds
      
        # compute the horizontal component of the dip-slip
        self.segments[0].vh = self.segments[0].ds
        self.segments[1].vh = -self.segments[1].ds*math.cos(self.beta)
        self.segments[2].vh = self.segments[2].ds*math.cos(self.gamma)
        # print 'v1h:', self.segments[0].vh 
        # print 'v3h:', self.segments[1].vh 
        # print 'v2h:', self.segments[2].vh
        # sys.exit()

# flower structure
class flower:
    def __init__(self,name1,ss1,sigmass1,H1,sigmaH1,D1,sigmaD1,\
        name2,ss2,sigmass2,H2,sigmaH2,D2,sigmaD2,\
        distss='Unif',distH='Unif',distD='Unif'):

        self.segments = [
                second(name = name2,ss = ss2,sigmass = sigmass2,D = D2, sigmaD = sigmaD2,H = H2,sigmaH = sigmaH2, distss=distss, distH=distH, distD=distD), # kink
                second(name = name1,ss = ss1,sigmass = sigmass1,D = D1, sigmaD = sigmaD1,H = H1,sigmaH = sigmaH1, distss=distss, distH=distH, distD=distD), # ramp
                ]
        self.Mseg = len(self.segments)
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))

    def conservation(self,seg,decol):
       
        #print
        #print 'ramp: ', self.segments[0].name
        v1 = seg.ds 
        gamma0 = seg.gamma
        self.quadrant = decol.quadrant
        
        # geoemtry of the ramp
        self.segments[1].fperp = self.segments[1].D + seg.fperp
        # Give a lenght to the segment
        self.segments[1].w = seg.w - self.segments[1].H 
        self.segments[1].L = math.sqrt(self.segments[1].D**2+(self.segments[1].H)**2)
        self.segments[1].dipr = math.atan2((self.segments[1].H),-self.segments[1].D)
        self.segments[1].dip = self.segments[1].dipr*180/math.pi
        self.gamma = math.atan2((self.segments[1].H),-self.segments[1].D*self.quadrant)
        self.segments[1].gamma = self.gamma
        # compute change of dip in between the two segments
        self.alpha = self.gamma - gamma0
        dip1,delta_dip = gamma0*180/math.pi, self.alpha*180/math.pi
        #print
        #print 'previous dip', dip1
        #print 'gamma, alpha', self.gamma*180/math.pi, delta_dip
        
        # geometry of the kink
        self.segments[0].fperp = self.segments[0].D + seg.fperp
        self.segments[0].w = seg.w - self.segments[0].H 
        self.segments[0].L = math.sqrt(self.segments[0].D**2+(self.segments[0].H)**2)
        self.segments[0].dipr = math.atan2((self.segments[0].H),-self.segments[0].D)
        self.segments[0].dip = self.segments[0].dipr*180/math.pi
        self.beta = math.atan2((self.segments[0].H),-self.segments[0].D*self.quadrant)
        #print 'D,H , beta: ',  self.segments[1].D, self.segments[1].H,  np.rad2deg(self.beta)
       
        # kink
        self.segments[0].ds = -v1*math.sin(self.alpha)/math.sin(self.beta-self.gamma)
        # ramp
        self.segments[1].ds = v1*math.sin(self.beta-self.gamma+self.alpha)/math.sin(self.beta-self.gamma)
        #print 'v2:', self.segments[1].ds 
        #print 'v3:', self.segments[0].ds
        
        # compute the horizontal component of the dip-slip
        self.segments[0].vh = -self.segments[0].ds*math.cos(self.beta)
        self.segments[1].vh = self.segments[1].ds*math.cos(self.gamma)
        #print 'v3h:', self.segments[0].vh 
        #print 'v2h:', self.segments[1].vh

# ramp segment
class ramp:
    def __init__(self,name,ss,sigmass,D,sigmaD,H, sigmaH,\
        distss='Unif',distH='Unif',distD='Unif'):

        self.segments = [
                second(name = name,ss = ss,sigmass = sigmass,D = D,sigmaD = sigmaD,H = H,sigmaH = sigmaH, distss=distss, distH=distH, distD=distD),    
                ]

        self.Mseg = len(self.segments)
        self.segments[0].fperp = self.segments[0].D 
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))
    
    def conservation(self,seg,decol):
        #print
        #print 'ramp: ', self.segments[0].name
        # get parameters previous segment
        v1 = seg.ds 
        gamma0 = seg.gamma
        self.quadrant = decol.quadrant
        # print self.quadrant
        
        # geometry of the ramp
        self.segments[0].fperp = self.segments[0].D + seg.fperp
        # Give a lenght to the segment
        self.segments[0].w = seg.w - self.segments[0].H 
        self.segments[0].L = math.sqrt(self.segments[0].D**2+(self.segments[0].H)**2)
        # Compute True Dip
        self.segments[0].dipr = math.atan2((self.segments[0].H),-self.segments[0].D)
        self.segments[0].dip = self.segments[0].dipr*180/math.pi
        # Compute dip in conservation of motion convention
        self.gamma = math.atan2((self.segments[0].H),-self.segments[0].D*self.quadrant)
        # print
        # print 'dip ramp:', self.segments[0].dip
        # print 'gamma:',self.gamma*180/math.pi
        # compute change of dip in between the two segments
        self.alpha = self.gamma - gamma0
        dip1,delta_dip = gamma0*180/math.pi, self.alpha*180/math.pi
        # print 'dip previous ramp', dip1
        # print 'alpha', delta_dip
        self.segments[0].ds = v1*math.cos(self.gamma-self.alpha)/math.cos(self.gamma)
        self.segments[0].vh = self.segments[0].ds*math.cos(self.gamma)
        # print 'ramp:', self.segments[0].ds, self.segments[0].vh
        # sys.exit()

# popup structure
class popup:
    def __init__(self,name4,H4,sigmaH4,D4,sigmaD4,\
        name5,H5,sigmaH5,D5,sigmaD5,\
        name6,ss6,sigmass6,D6,sigmaD6,\
        distss='Unif',distH='Unif',distD='Unif'):

        self.segments = [
                # ramp poppup
                second(name = name4,ss = 0,sigmass = 0,D = D4, sigmaD = sigmaD4,H = H4,sigmaH = sigmaH4), 
                # kink popup
                second(name = name5,ss = 0,sigmass = 0,D = D5, sigmaD = sigmaD5,H = H5,sigmaH = sigmaH5), 
                # decol
                second(name = name6,ss = ss6,sigmass = sigmass6,D = D6, sigmaD = sigmaD6,H = 0,sigmaH = 0)  
                ]
        self.Mseg = len(self.segments)
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))
    
    def conservation(self,seg,decol):
       
        #print
        #print 'ramp: ', self.segments[0].name
        v1 = seg.ds 
        gamma0 = seg.gamma
        self.quadrant = decol.quadrant
        
        # geoemtry of the ramp
        self.segments[0].fperp = self.segments[0].D + seg.fperp
        # Give a lenght to the segment
        self.segments[0].w = seg.w - self.segments[0].H 
        self.segments[0].L = math.sqrt(self.segments[0].D**2+(self.segments[0].H)**2)
        self.segments[0].dipr = math.atan2((self.segments[0].H),-self.segments[0].D)
        self.segments[0].dip = self.segments[0].dipr*180/math.pi
        self.gamma = math.atan2((self.segments[0].H),-self.segments[0].D*self.quadrant)
        self.segments[0].gamma = self.gamma 
        # compute change of dip in between the two segments
        self.alpha = self.gamma - gamma0
        dip1,delta_dip = gamma0*180/math.pi, self.alpha*180/math.pi
        #print 'dip1', dip1
        #print 'gamma, alpha', np.rad2deg(self.gamma), delta_dip
        
        # geometry of the kink
        self.segments[1].fperp = self.segments[1].D + seg.fperp
        self.segments[1].w = seg.w - self.segments[1].H 
        self.segments[1].L = math.sqrt(self.segments[1].D**2+(self.segments[1].H)**2)
        self.segments[1].dipr = math.atan2((self.segments[1].H),-self.segments[1].D)
        self.segments[1].dip = self.segments[1].dipr*180/math.pi
        self.beta = math.atan2((self.segments[1].H),-self.segments[1].D*self.quadrant)
        #print 'D,H , beta: ',  self.segments[1].D, self.segments[1].H,  self.segments[1].dip
        
        # geometry of the ramp
        self.segments[2].fperp = self.segments[2].D + seg.fperp
        self.segments[2].w = seg.w - self.segments[2].H 
        self.segments[2].L = math.sqrt(self.segments[2].D**2+(self.segments[2].H)**2)
        self.segments[2].dipr = math.atan2((self.segments[2].H),-self.segments[2].D)
        self.segments[2].dip = self.segments[2].dipr*180/math.pi
        self.segments[2].ds = 0.
        self.segments[2].vh = 0.
        #print 'D,H , beta: ',  self.segments[1].D, self.segments[1].H,  self.segments[1].dip
       
        # kink
        self.segments[1].ds = -v1*math.sin(self.alpha)/math.sin(self.beta-self.gamma)
        # ramp
        self.segments[0].ds = v1*math.sin(self.beta-self.gamma+self.alpha)/math.sin(self.beta-self.gamma)
        #print 'v3:', self.segments[1].ds 
        #print 'v2:', self.segments[0].ds
        
        # compute the horizontal component of the dip-slip
        self.segments[1].vh = -self.segments[1].ds*math.cos(self.beta)
        self.segments[0].vh = self.segments[0].ds*math.cos(self.gamma)
        #print 'v3h:', self.segments[1].vh 
        #print 'v2h:', self.segments[0].vh
        #sys.exit()

# creeping segment: problem bc half-infinite segment
class creeping:
    def __init__(self,name,ss,sigmass,H,sigmaH,D,sigmaD,\
        distss='Unif',distH='Unif',distD='Unif'):

        self.segments = [
                    second(name = name,ss = ss,sigmass = sigmass,D = D,sigmaD = sigmaD,H = H,sigmaH = sigmaH, distss=distss,distH=distH, distD=distD)
                    ]
        self.Mseg = len(self.segments)
        self.segments[0].fperp = self.segments[0].D 
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))
    
    def conservation(self,seg,decol):
        w0 = decol.w

        self.segments[0].fperp = self.segments[0].D 
        self.segments[0].w = w0 - self.segments[0].H
        self.segments[0].L = self.segments[0].H
        self.segments[0].dipr = math.pi/2
        self.segments[0].dip = 90
        # fix the ds on the creeping seg to be zero  
        self.segments[0].ds = decol.ds/10e14
        self.segments[0].vh = decol.ds/10e14

class bookshelf:
    def __init__(self,name,ds,sigmads,D,sigmaD,distshort='Unif',distD='Unif'):
        self.name = name
        self.D = D
        self.sigmaD = sigmaD
        self.ds = ds
        self.sigmads = sigmads
        self.distshort = distshort
        self.distD = distD
        
        # 2 parametres ici
        self.Mker = 2 

    def displacement(self,yp):
        # def heaviside fic
        uperp = (yp-self.D)/self.D
        uperp[yp<self.D], uperp[yp>0] = 0,-1
        u = np.zeros((yp.size,3))
        u[:,0],u[:,1],u[:,2] = np.zeros(len(yp)),self.ds*uperp,np.zeros(len(yp))
        return u
    
    def info(self):
        print self.name
        print "#      ds       D        "
        print ('%4.1f mm/yr     %4i km  '%(self.ds,self.D))
        print "#     sigmads      "
        print ('   %4.3f mm/yr %4i km' %(self.sigmads, self.sigmaD))
        print

#class main structure
class main(segment):
    def __init__(self,name,ss,sigmass,H,sigmaH,ds,sigmads,x,y,dip,\
        distss,distH,distshort):
        segment.__init__(self,name,ss,sigmass,H,sigmaH,distss,distH)
        # ds: shortening 
        self.ds = ds
        self.sigmads = sigmads
        self.distshort = distshort
        # E,W position in km
        self.x = x
        self.y = y
        # D: distance to the main fault
        self.D = 0
        # L fixe and large for main structure: half-infinite dislocation
        self.L = 660
        # diping fix horizontal to model the shortening
        self.dip = dip
        self.dipr = self.dip*math.pi/180
        # 3 unknow parameter for the main segment
        self.Mker = 3
        if math.cos(self.dipr) >= 0:
            self.quadrant = 1
        else:
            self.quadrant = -1
        # initializt sst
        self.sst=0
    
    def info(self):
        print self.name
        print "#    ss          ds              w    "
        print ('%4.1f mm/yr         %4.1f mm/yr         %4i km'%(self.sst,self.ds,self.w))
        print "#    sigmass     sigmads    sigmaw    "
        print ('%4.3f mm/yr         %4.3f mm/yr         %4.3f km'%(self.sigmass,self.sigmads,self.sigmaw))
        print

    def write(self,fid):
        fid.write('{}\n'.format(self.name))
        fid.write('# ss ds w\n')
        np.savetxt(fid, np.vstack([self.ss,self.ds,self.w]).T, fmt = '%.6f' ,delimiter = '\t', newline = '\n' )
        np.savetxt(fid,np.vstack([self.sigmass,self.sigmads,self.sigmaH]).T,fmt = '%.6f',delimiter = '\t', newline = '\n')

# class secondary structure
class second(segment):
    def __init__(self,name,ss,sigmass,H,sigmaH,D,sigmaD,distss,distH,distD):
        segment.__init__(self,name,ss,sigmass,H,sigmaH,distss,distH)
        # D: distance to the main fault
        self.D = D
        self.sigmaD = sigmaD
        # prior distribution for D
        self.distD=distD
        # initialise at 0, ll be then define by the conservation of motion
        self.ds, self.dip, self.L  = 0, 0, 0 
        self.sigmads = 0
        # 3 free parameters for a secondary segment
        self.Mker = 3
        self.sigmadip=0
        self.sigmaL=0

    def info(self):
        print self.name
        print "#   ss  ds  Width  dip L  w"
        print ('%4.1f mm/yr     %4.1f mm/yr     %4i km      %4.1f   %4.1f km    %4i km' %(self.ss,self.ds,self.fperp,self.dip,self.L,self.w))
        print "#   sigmass  sigmads  sigmaD   sigmaw"
        print ('%4.3f mm/yr   %4.3f mm/yr    %4.3f km       %4.3f km' %(self.sigmass,self.sigmads,self.sigmaD,self.sigmaH))
        print

    def write(self,fid):
        # save inversion results
        fid.write('{}\n'.format(self.name))
        fid.write('#  ss ds D fperp dip L H w\n')
        np.savetxt(fid, np.vstack([self.ss,self.ds,self.D,self.fperp,self.dip,self.L,self.H,self.w]).T, fmt = '%.6f' ,delimiter = '\t', newline = '\n' )
        np.savetxt(fid,np.vstack([self.sigmass,self.sigmads, self.sigmaD, self.sigmaD,self.sigmadip,self.sigmaL,self.sigmaH,self.sigmaH]).T,fmt = '%.6f',delimiter = '\t', newline = '\n')

class topo:
    def __init__(self,name,wdir,filename,color,width,scale=1):
        self.name = name
        self.wdir = wdir
        self.filename = filename
        self.color = color
        self.width = width
        self.scale=scale

    def load(self,flt):
        fmodel = flt.fmodel
        profile = flt.profiles
        fname = file(self.wdir+self.filename)
        x,y,z = np.loadtxt(fname,comments = '#',unpack = True,dtype = 'f,f,f')
        xp = (x-fmodel[0].x)*profile.s[0]+(y-fmodel[0].y)*profile.s[1]
        yp = (x-fmodel[0].x)*profile.n[0]+(y-fmodel[0].y)*profile.n[1]
        index = np.nonzero((xp>profile.xpmax)|(xp<profile.xpmin)|(yp>profile.ypmax)|(yp<profile.ypmin))
        
        self.x,self.y,self.z = np.delete(x,index),np.delete(y,index),np.delete(z,index)*self.scale

class seismi:
    def __init__(self,name,wdir,filename,color,width,scale=1):
        self.name = name
        self.wdir = wdir
        self.filename = filename
        self.color = color
        self.width = width
        self.scale=scale

    def load(self,flt):
        fmodel = flt.fmodel
        profile = flt.profiles
        fname = file(self.wdir+self.filename)
        x,y,z,mw = np.loadtxt(fname,comments = '#',unpack = True,dtype = 'f,f,f,f')
        xp = (x-fmodel[0].x)*profile.s[0]+(y-fmodel[0].y)*profile.s[1]
        yp = (x-fmodel[0].x)*profile.n[0]+(y-fmodel[0].y)*profile.n[1]
        index = np.nonzero((xp>profile.xpmax)|(xp<profile.xpmin)|(yp>profile.ypmax)|(yp<profile.ypmin))
        self.x,self.y,self.z,self.mw = np.delete(x,index),np.delete(y,index),np.delete(z,index)*self.scale,np.delete(mw,index)

class moho:
    def __init__(self,name,wdir,filename,color,width,scale=1):
        self.name = name
        self.wdir = wdir
        self.filename = filename
        self.color = color
        self.width = width
        self.scale=scale
        
    def load(self,flt):
        fmodel = flt.fmodel
        profile = flt.profiles
        fname = file(self.wdir+self.filename)
        x,y,z = np.loadtxt(fname,comments = '#',unpack = True,dtype = 'f,f,f')
        xp = (x-fmodel[0].x)*profile.s[0]+(y-fmodel[0].y)*profile.s[1]
        yp = (x-fmodel[0].x)*profile.n[0]+(y-fmodel[0].y)*profile.n[1]
        index = np.nonzero((xp>profile.xpmax)|(xp<profile.xpmin)|(yp>profile.ypmax)|(yp<profile.ypmin))
        self.x,self.y,self.z = np.delete(x,index),np.delete(y,index),np.delete(z,index)*self.scale
