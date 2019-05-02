import numpy as np
import math
import sys
import time, logging

class profile:
    """ 
    profile class: Load profile 
    Parameters: 
    name: name profile
    x,y: reference point: everything is refered to this point 
    l,w: length, width profile
    strike: strike profile (default: defined by inversion class)
    proj=[east, notth, up]: average LOS projection into east, north, up used for plots [default: None]
    type:  std - plot mean and standard deviation InSAR;
    distscale - scatter plot with color scale function of the profile-parallel distance;
    stdscat - plot scatter + standar deviation. 
    """

    def __init__(self,x,y,l,w,strike=None,proj=None,type=None,name=''):
        # profile parameters
        self.name = name
        self.x = x
        self.y = y
        self.l = l
        self.w = w
        self.proj = proj
        if strike > 0:
            self.strike=strike-180
        else:
            self.strike=strike
        self.typ=type

class inversion:
    """
    Inversion class
    structure: define all structures
    strike: azimuth inversion profile (ie perpandicular to the main fault)
    profile: profile class 
    fullcov      :   if True, estimates a full covariance matrix for InSAR data
                     if False, returns a diagonal matrix with sigmad on the diagonal
    maskcov: give a mask for covariance estimation on sub-sampled data in perpendicular 
    distance to the center of profile km eg. [0,20,40,50] (default: None)
    rampcov: remove ramp before covariance estimation. eg. lin, quad, cub (default: lin)
    name: give a name (Optional)
    depthmax: maximum depth for plot (Optinala)
    """
    
    def __init__(self,structures,strike,profile, 
        fullcov=False, maskcov=None, rampcov=False, name='', depthmax=None):
        self.name = name
        self.fmodel = []
        self.structures = structures
        self.profile = profile
        # ref point inversion is ref point profile
        self.x, self.y = self.profile.x, self.profile.y

        # define north to the right of the profile 
        if strike > 0:
            self.strike=strike-180
        else:
            self.strike=strike
        self.str = np.deg2rad(self.strike)
        self.s = [math.sin(self.str),math.cos(self.str),0]
        self.n = [math.cos(self.str),-math.sin(self.str),0]
        self.profile.s = self.s
        self.profile.n = self.n
        self.profile.str=self.str
        self.profile.strike=self.strike

        self.fullcov = fullcov
        self.maskcov = maskcov
        self.rampcov = rampcov

        self.depthmax=depthmax

class segment:
    def __init__(self,name,ss,sigmass,H,sigmaH,distss,distH):
        self.name = name
        # ss: strike-slip
        self.ss = ss
        self.sigmass = sigmass
        # H: vertical distance
        self.H = H
        self.sigmaH = sigmaH
        # # w: locking depth
        # self.w, self.sigmaw = 0, 0
        # horizontal distance to the main fault
        self.fperp = 0
        self.distss=distss
        self.distH=distH

    def displacement(self,yp):
        """
        #######################################################
        # Surface Displacements Due to Edge Dislocation
        # Paul Segall 1996
        # Input:
        #       yp = coordinates of the observation points
        #       w1 = depth of updip end;
        #       L = length downdip;
        #######################################################
        """

        L = self.L
        self.dipr = np.deg2rad(self.dip)
        shift = -self.fperp*np.ones(yp.shape)
        z = 0
        w1 = self.w
        w2 = self.w+L*math.sin(self.dipr)
        zeta1 = (yp+shift)/w1
        zeta2 = (yp+shift-L*math.cos(self.dipr))/w2
        denom1 = 1+zeta1**2
        denom2 = 1+zeta2**2
        u = np.zeros((yp.size,3))

        # vertical displacement
        uv =  (1./math.pi) * (math.sin(self.dipr)*np.arctan(zeta1) + (math.cos(self.dipr)+math.sin(self.dipr)*zeta1)/denom1 - math.sin(self.dipr)*np.arctan(zeta2) - (math.cos(self.dipr)+math.sin(self.dipr)*zeta2)/denom2 )
    
        # horizontal displacement
        uperp =  -(1./math.pi) * (math.cos(self.dipr)*np.arctan(zeta1) + (math.sin(self.dipr) -math.cos(self.dipr)*zeta1)/denom1 - math.cos(self.dipr)*np.arctan(zeta2) - (math.sin(self.dipr) - math.cos(self.dipr)*zeta2)/denom2)

        #######################################################
        # Surface Displacements Due to an Half-infinite dipping dislocation
        # Paul Segall 1996
        #######################################################

        # half-infinite dislocation
        upar =  (1./(2*math.pi)) * (np.arctan2(math.sin(self.dipr)*(yp+shift)+ \
            math.cos(self.dipr)*(self.w), math.sin(self.dipr)*(self.w) - math.cos(self.dipr)*(yp+shift)) + \
            np.arctan2(math.sin(self.dipr)*(yp+shift) + math.cos(self.dipr)*self.w, \
                -math.cos(self.dipr)*(yp+shift) + math.sin(self.dipr)*self.w ))  

        if L < 660:

            shift = shift - L*math.cos(self.dipr)
            w1 = self.w + L*math.sin(self.dipr)

            uparf = (1./(2*math.pi)) * (np.arctan2(math.sin(self.dipr)*(yp+shift)+ \
            math.cos(self.dipr)*(w1), math.sin(self.dipr)*(w1) - math.cos(self.dipr)*(yp+shift)) + \
            np.arctan2(math.sin(self.dipr)*(yp+shift) + math.cos(self.dipr)*w1, \
                -math.cos(self.dipr)*(yp+shift) + math.sin(self.dipr)*w1 )) 

            upar = upar - uparf


        u[:,0],u[:,1],u[:,2] = self.ss*upar,self.ds*uperp,self.ds*uv

        return u

#class main structure
class main(segment):
    def __init__(self,name,ss,sigmass,sigmaH,short,sigmashort,dip,\
        distss,distH,distshort,distD,distL,distdip,
        H=0,D=0,sigmaD=0,L=660.,sigmaL=0.,sigmadip=0):

        # H is initially zero bc main fault 
        segment.__init__(self,name,ss,sigmass,H,sigmaH,distss,distH)

        self.vh = short
        self.sigmavh = sigmashort
        self.distvh = distshort

        # position to the center of profile
        self.D, self.fperp = D, D
        self.sigmaD = sigmaD
        self.distD = distD
        # L fixe and large for main structure: half-infinite dislocation
        self.L = L
        self.sigmaL = sigmaL
        self.distL = distL
        self.dip = dip
        self.sigmadip = sigmadip
        self.distdip = distdip
        self.dipr = np.deg2rad(self.dip)

        self.Mker = 6
        if math.cos(self.dipr) >= 0:
            self.quadrant = 1
            self.gamma = np.deg2rad(self.dip)
        else:
            self.quadrant = -1
            self.gamma = np.deg2rad(180 - self.dip)

        self.ds = self.vh/np.cos(self.gamma)
        self.sst = ss

        if self.ds > 0. and self.L > 50. and ((self.dip == 0) or  (self.dip == 180)):
            print('Warning! Setting a long dipping thrust will create long-wavelength vertical Displacements')
            print('Warning! Set dip to 0 or 180 to model a large scale shortening')
    
    def info(self):
        print(self.name)
        line=['ss', 'ss_total',  'short',  'Depth', 'D' , 'Length', 'dip']
        print('{:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*line))
        line=[self.ss, self.sst, self.vh, self.w,self.D,self.L,self.dip]
        print('{:8.2f} {:8.2f} {:8.2f} {:8.0f} {:8.0f} {:8.0f} {:8.0f}'.format(*line))
        line=['sigmass', 'sigmashort',  'sigmaDepth',  'sigmaD', 'sigmaL','sigmadip']
        print('{:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*line))
        line=[self.sigmass, self.sigmavh, self.sigmaw,self.sigmaD,self.sigmaL,self.sigmadip]
        print('{:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}'.format(*line))
        # print("#    ss      ss_total         ds              w      D       L    dip")
        # print('%4.1f mm/yr   %4.1f mm/yr      %4.1f mm/yr      %4i km     %4i km      %4i km    %4.1f'%(self.ss, self.sst, self.ds, self.w,self.D,self.L,self.dip))
        # print("#    sigmass     sigmads    sigmaw    sigmaD    sigmaL    sigmadip")
        # print('%4.3f mm/yr         %4.3f mm/yr         %4.3f km         %4.3f km         %4.3f km       %4.3f'%(self.sigmass, self.sigmads, self.sigmaw,self.sigmaD,self.sigmaL,self.sigmadip))

    def write(self,fid):
        fid.write('{}\n'.format(self.name))
        fid.write('# ss sst short depth\n')
        np.savetxt(fid, np.vstack([self.ss,self.sst,self.vh,self.w]).T, fmt = '%.6f' ,delimiter = '\t', newline = '\n' )
        fid.write('# sigma_ss sigma_short sigma_depth\n')
        np.savetxt(fid,np.vstack([self.sigmass,self.sigmavh,self.sigmaH]).T,fmt = '%.6f',delimiter = '\t', newline = '\n')

# class secondary structure
class second(segment):
    def __init__(self,name,ss,sigmass,H,sigmaH,D,sigmaD,distss,distH,distD):
        segment.__init__(self,name,ss,sigmass,H,sigmaH,distss,distH)
        # D: distance to the main fault
        self.D = D
        self.sigmaD = sigmaD
        self.distD=distD
        # initialise at 0, will be then define by the conservation of motion
        self.ds, self.dip, self.L  = 0, 0, 0 
        self.sigmads = 0
        # 3 free parameters for a secondary segment
        self.Mker = 3
        self.sigmadip=0
        self.sigmaL=0

    def info(self):
        print(self.name)
        line=['ss', 'ds', 'D', 'Width', 'H', 'dip', 'Length' , 'Depth']
        print('{:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}'.format(*line))
        line=[self.ss,self.ds,self.D,self.fperp,self.H,self.dip,self.L,self.w]
        print('{:8.2f} {:8.2f} {:8.2f} {:8.0f} {:8.0f} {:8.0f} {:8.0f} {:8.0f}'.format(*line))
        line=['sigmass', 'sigmads',  'sigmaD',  'sigmaH']
        print('{:>8} {:>8} {:>8} {:>8}'.format(*line))
        line=[self.sigmass,self.sigmads,self.sigmaD,self.sigmaH]
        print('{:8.2f} {:8.2f} {:8.2f} {:8.2f}'.format(*line))
        # print('%4.1f mm/yr     %4.1f mm/yr     %4i km      %4.1f   %4.1f km    %4i km' %(self.ss,self.ds,self.fperp,self.dip,self.L,self.H))
        # print("#   sigmass  sigmads  sigmaD   sigmaH")
        # print('%4.3f mm/yr   %4.3f mm/yr    %4.3f km       %4.3f km' %(self.sigmass,self.sigmads,self.sigmaD,self.sigmaH))

    def write(self,fid):
        # save inversion results
        fid.write('{}\n'.format(self.name))
        fid.write('#  ss ds D fperp dip L H w\n')
        np.savetxt(fid, np.vstack([self.ss,self.ds,self.D,self.fperp,self.dip,self.L,self.H,self.w]).T, fmt = '%.6f' ,delimiter = '\t', newline = '\n' )
        np.savetxt(fid,np.vstack([self.sigmass,self.sigmads, self.sigmaD, self.sigmaD,self.sigmadip,self.sigmaL,self.sigmaH,self.sigmaH]).T,fmt = '%.6f',delimiter = '\t', newline = '\n')


# one decollement as first structure 
class mainfault:
    def __init__(self,name,ss,sigmass,short,sigmashort,w,sigmaw,dip,\
        distss='Unif',distshort='Unif',distH='Unif',distL='Unif',distD='Unif',distdip='Unif',
        D=0.,sigmaD=0.,L=660.,sigmaL=0.,sigmadip=0.):
       
        # init H to zero because main fault
        self.segments = [
                main(name = name, D = D, sigmaD = sigmaD, ss = ss, sigmass = sigmass, short = short, 
                    sigmashort = sigmashort, H = 0, sigmaH = sigmaw, dip = dip, sigmadip=sigmadip, L = L, sigmaL = sigmaL, 
                    distss=distss,distshort=distshort,distH=distH,distL=distL,distD=distD,distdip=distdip),
                ]

        self.Mseg = 1
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))
        self.segments[0].w = w
        self.segments[0].sigmaw = sigmaw
        self.winit,self.segments[0].winit = w, w
        self.segments[0].w = w
        
    def conservation(self):
        self.segments[0].w = self.winit - self.segments[0].H
        
        # self.segments[0].gamma = self.segments[0].dipr
        if math.cos(self.segments[0].dipr) >= 0:
            self.segments[0].gamma = np.deg2rad(self.segments[0].dip)
        else:
            self.segments[0].gamma = np.deg2rad(180 - self.segments[0].dip)
        self.segments[0].ds = self.segments[0].vh/np.cos(self.segments[0].gamma)

        self.segments[0].fperp = self.segments[0].D
        # print('main F:{}, dip:{}, ds: {}, vh:{}'.format(self.segments[0].name,self.segments[0].dip,self.segments[0].ds,self.segments[0].vh))

# ramp segment
class ramp:
    def __init__(self,name,ss,sigmass,D,sigmaD,H, sigmaH,\
        distss='Unif',distH='Unif',distD='Unif'):

        self.segments = [
                second(name = name,ss = ss,sigmass = sigmass,D = D,
                    sigmaD = sigmaD,H = H,sigmaH = sigmaH, 
                    distss=distss, distH=distH, distD=distD),    
                ]

        self.Mseg = len(self.segments)
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))
    
    def conservation(self,seg,decol):
        v1 = seg.ds 
        self.quadrant = decol.quadrant
        self.segments[0].fperp = self.segments[0].D + seg.fperp
        
        self.segments[0].w = seg.w - self.segments[0].H 
        self.segments[0].L = math.sqrt(self.segments[0].D**2+(self.segments[0].H)**2)
        self.segments[0].dipr = math.atan2((self.segments[0].H),-self.segments[0].D)
        self.segments[0].dip = np.rad2deg(self.segments[0].dipr)
        # self.gamma = math.atan2((self.segments[0].H),-self.segments[0].D*self.quadrant)
        self.gamma = math.atan2((self.segments[0].H),-self.segments[0].D*self.quadrant)
        self.segments[0].gamma = self.gamma
        # compute change of dip in between the two segments
        # # dip0=0 --> alpha = gamma
        self.alpha = self.gamma - seg.gamma
        self.segments[0].ds = v1*math.cos(self.gamma-self.alpha)/math.cos(self.gamma)
        self.segments[0].vh = self.segments[0].ds*math.cos(self.gamma)
        # print('MF:{}, dip0:{}, gamma0:{}, ds0:{}, vh0:{}'.format(seg.name, seg.dip, np.rad2deg(seg.gamma), v1, seg.vh))
        # print('ramp:{}, gamma:{}, alpha:{}, ds:{}, vh:{}:'.format(self.segments[0].name, np.rad2deg(self.gamma), np.rad2deg(self.alpha), self.segments[0].ds, self.segments[0].vh))
        # sys.exit()

# flower structure as first structure
class mainflower:
    def __init__(self,name,ss,sigmass,short,sigmashort,w,sigmaw,dip,\
        name2,H2,sigmaH2,ss2,sigmass2,D2,sigmaD2,\
        name3,H3,sigmaH3,ss3,sigmass3,D3,sigmaD3,\
        distss='Unif',distshort='Unif',distH='Unif',distD='Unif',distL='Unif',distdip='Unif',
        D=0.,sigmaD=0.,L=660.,sigmaL=0.,sigmadip=0.):

        self.segments = [
                main(name = name, D = D, sigmaD = sigmaD, ss = ss, sigmass = sigmass, short = short,
                   sigmashort = sigmashort, H = 0, sigmaH = sigmaw, dip = dip, sigmadip= sigmadip, L = L, sigmaL = sigmaL,
                   distss=distss, distshort=distshort, distH=distH, distL=distL,distD=distD, distdip=distdip),
                second(name = name3, ss = ss3, sigmass = sigmass3, D = D3, 
                   sigmaD = sigmaD3, H = H3, sigmaH = sigmaH3, 
                   distss= distss, distH=distH, distD=distD), 
                second(name = name2,ss = ss2, sigmass = sigmass2, D = D2,
                   sigmaD = sigmaD2,H = H2,sigmaH = sigmaH2, 
                   distss=distss, distH=distH, distD=distD), 
                ]

        self.winit,self.segments[0].winit = w, w
        self.segments[0].ds = self.segments[0].vh/np.cos(self.segments[0].dipr)
        self.Mseg = len(self.segments)
        self.segments[0].w,self.segments[1].w,self.segments[2].w = w, w - self.segments[1].H, w - self.segments[2].H
        self.segments[0].sigmaw = sigmaw
        self.segments[0].fperp, self.segments[1].fperp,self.segments[2].fperp = self.segments[0].D, self.segments[1].D, self.segments[2].D  
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))
        
    def conservation(self):

        # update LD
        self.segments[0].w = self.winit - self.segments[0].H
        self.quadrant = self.segments[0].quadrant
        self.segments[0].fperp =  self.segments[0].D 

        # Get model parameters for the main fault
        self.segments[0].ds = self.segments[0].vh/np.cos(self.segments[0].gamma)
        # v1 is ds on the main fault. shortening if flat seg only
        v1 = self.segments[0].ds

        self.segments[1].fperp =  self.segments[1].D + self.segments[0].fperp
        self.segments[1].w = self.segments[0].w - self.segments[1].H
        self.segments[1].L = math.sqrt(self.segments[1].D**2+(self.segments[1].H)**2)
        self.segments[1].dipr = math.atan2((self.segments[1].H),-self.segments[1].D)
        self.segments[1].dip = np.rad2deg(self.segments[1].dipr)
        # si dip0 = 0, D>0, H>0 --> -D pour que beta entre 90 et 180
        # si dip0 = 180, D<0, H>0 --> D pour que beta entre 90 et 180
        self.beta = math.atan2((self.segments[1].H),-self.segments[1].D*self.quadrant)

        self.segments[2].fperp = self.segments[2].D + self.segments[0].fperp
        self.segments[2].w = self.segments[0].w - self.segments[2].H
        self.segments[2].L = math.sqrt(self.segments[2].D**2+(self.segments[2].H)**2)
        # si dip0 = 0, D<0, H>0 --> -D pour que gamma entre 0 et 90
        # si dip0 = 180, D>0, H>0 --> D pour que beta entre 0 et 90
        self.segments[2].dipr = math.atan2((self.segments[2].H),-self.segments[2].D)
        self.segments[2].dip = np.rad2deg(self.segments[2].dipr)
        # save gamma in segment attribute for next structure
        self.segments[2].gamma = math.atan2((self.segments[2].H),-self.segments[2].D*self.quadrant)
        self.gamma = self.segments[2].gamma
        # alpha between 0 and 90
        self.alpha = self.gamma - self.segments[0].gamma
        # kink
        self.segments[1].ds = -self.segments[0].vh*math.sin(self.alpha)/math.sin(self.beta-self.gamma)
        # ramp
        self.segments[2].ds = self.segments[0].vh*math.sin(self.beta-self.gamma+self.alpha)/math.sin(self.beta-self.gamma)
        # # dip0=0 --> alpha = gamma
        # self.segments[1].ds = -v1*math.sin(self.gamma)/math.sin(self.beta-self.gamma)
        # self.segments[2].ds = v1*math.sin(self.beta)/math.sin(self.beta-self.gamma)      
        
        # compute the horizontal component of the dip-slip
        self.segments[1].vh = -self.segments[1].ds*math.cos(self.beta)
        self.segments[2].vh = self.segments[2].ds*math.cos(self.gamma)
        # print('MP:{}, dip: {}, v1: {}, V1h:{}:'.format(self.segments[0].name,self.segments[0].dip,v1,self.segments[0].vh ))
        # print('Kink:{}, D3:{}, H3:{}, beta:{}, V3:{}, Vh3:{}'.format(self.segments[1].name,self.segments[1].D, self.segments[1].H, np.rad2deg(self.beta),self.segments[1].ds,self.segments[1].vh))
        # print('Ramp:{}, D2:{}, H2:{}, alpha:{}, gamma:{}, V2:{}, Vh2:{}'.format(self.segments[2].name,self.segments[2].D, self.segments[2].H, np.rad2deg(self.alpha), np.rad2deg(self.gamma),self.segments[2].ds,self.segments[2].vh))

# flower structure
class flower:
    def __init__(self,name1,ss1,sigmass1,H1,sigmaH1,D1,sigmaD1,\
        name2,ss2,sigmass2,H2,sigmaH2,D2,sigmaD2,\
        distss='Unif',distH='Unif',distD='Unif'):

        self.segments = [
                second(name = name2,ss = ss2,sigmass = sigmass2,D = D2, 
                sigmaD = sigmaD2,H = H2,sigmaH = sigmaH2, 
                distss=distss, distH=distH, distD=distD), # kink
                second(name = name1,ss = ss1,sigmass = sigmass1,D = D1, 
                sigmaD = sigmaD1,H = H1,sigmaH = sigmaH1, 
                distss=distss, distH=distH, distD=distD), # ramp
                ]
        self.Mseg = len(self.segments)
        self.Mker = sum(map((lambda x: getattr(x,'Mker')),self.segments))

    def conservation(self,seg,decol):
       
        #print
        #print 'ramp: ', self.segments[0].name
        v1 = seg.ds 
        self.quadrant = decol.quadrant
        self.segments[1].fperp = self.segments[1].D + seg.fperp
        self.segments[1].w = seg.w - self.segments[1].H 
        self.segments[1].L = math.sqrt(self.segments[1].D**2+(self.segments[1].H)**2)
        self.segments[1].dipr = math.atan2((self.segments[1].H),-self.segments[1].D)
        self.segments[1].dip = np.rad2deg(self.segments[1].dipr)
        self.gamma = math.atan2((self.segments[1].H),-self.segments[1].D*self.quadrant)
        self.segments[1].gamma = self.gamma
        # compute change of dip in between the two segments
        self.alpha = self.gamma - seg.gamma
        # print(np.rad2deg(self.gamma), np.rad2deg(seg.gamma))
        # geometry of the kink
        self.segments[0].fperp = self.segments[0].D + seg.fperp
        self.segments[0].w = seg.w - self.segments[0].H 
        self.segments[0].L = math.sqrt(self.segments[0].D**2+(self.segments[0].H)**2)
        self.segments[0].dipr = math.atan2((self.segments[0].H),-self.segments[0].D)
        self.segments[0].dip = np.rad2deg(self.segments[0].dipr)
        self.beta = math.atan2((self.segments[0].H),-self.segments[0].D*self.quadrant)
        # kink
        self.segments[0].ds = -v1*math.sin(self.alpha)/math.sin(self.beta-self.gamma)
        # ramp
        self.segments[1].ds = v1*math.sin(self.beta-self.gamma+self.alpha)/math.sin(self.beta-self.gamma)
        # compute the horizontal component of the dip-slip
        self.segments[0].vh = -self.segments[0].ds*math.cos(self.beta)
        self.segments[1].vh = self.segments[1].ds*math.cos(self.gamma)
        # print('Previous F:{}, dip:{}, ds:{}'.format( seg.name, np.rad2deg(seg.gamma), seg.ds))
        # print('Kink:{}, D:{}, H:{}, beta:{}, ds:{}, vh:{}'.format(self.segments[0].name,self.segments[0].D, self.segments[0].H, np.rad2deg(self.beta), self.segments[0].ds,self.segments[0].vh ))
        # print('Ramp:{}, D:{}, H:{}, gamma:{}, alpha:{}, ds:{}, vh:{}'.format(self.segments[1].name,self.segments[1].D, self.segments[1].H, np.rad2deg(self.gamma),np.rad2deg(self.alpha), self.segments[1].ds,self.segments[1].vh ))
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
        self.segments[0].dip = np.rad2deg(self.segments[0].dipr)
        self.gamma = math.atan2((self.segments[0].H),-self.segments[0].D*self.quadrant)
        self.segments[0].gamma = self.gamma 
        # compute change of dip in between the two segments
        self.alpha = self.gamma - gamma0
        dip1,delta_dip = np.rad2deg(gamma0), np.rad2deg(self.alpha)
        #print 'dip1', dip1
        #print 'gamma, alpha', np.rad2deg(self.gamma), delta_dip
        
        # geometry of the kink
        self.segments[1].fperp = self.segments[1].D + seg.fperp
        self.segments[1].w = seg.w - self.segments[1].H 
        self.segments[1].L = math.sqrt(self.segments[1].D**2+(self.segments[1].H)**2)
        self.segments[1].dipr = math.atan2((self.segments[1].H),-self.segments[1].D)
        self.segments[1].dip = np.rad2deg(self.segments[1].dipr)
        self.beta = math.atan2((self.segments[1].H),-self.segments[1].D*self.quadrant)
        #print 'D,H , beta: ',  self.segments[1].D, self.segments[1].H,  self.segments[1].dip
        
        # geometry of the ramp
        self.segments[2].fperp = self.segments[2].D + seg.fperp
        self.segments[2].w = seg.w - self.segments[2].H 
        self.segments[2].L = math.sqrt(self.segments[2].D**2+(self.segments[2].H)**2)
        self.segments[2].dipr = math.atan2((self.segments[2].H),-self.segments[2].D)
        self.segments[2].dip = np.rad2deg(self.segments[2].dipr)
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
        print(self.name)
        print("#      ds       D        ")
        print('%4.1f mm/yr     %4i km  '%(self.ds,self.D))
        print("#     sigmads      ")
        print ('   %4.3f mm/yr %4i km' %(self.sigmads, self.sigmaD))

class topo:
    """ 
    topo class: Load topographic file for plot
    Parameters: 
    filename: name input file
    name: name file for plot
    wdir: path input file
    scale: scale values
    color
    topomin,topomax
    plotminmax: option to also plot min max topo within bins
    """

    def __init__(self,name,wdir,filename,color='black',width=1.,
        scale=1,topomin=None,topomax=None,plotminmax=False):

        self.name = name
        self.wdir = wdir
        self.filename = filename
        self.color = color
        self.scale=scale
        self.plotminmax = plotminmax
        self.topomin = topomin
        self.topomax = topomax
        self.width = width


    def load(self,flt):
        logger = flt.logger
        try:
            fmodel = flt.fmodel
            profile = flt.profile
            fname = file(self.wdir+self.filename)
            x,y,z = np.loadtxt(fname,comments = '#',unpack = True,dtype = 'f,f,f')
            xp = (x-profile.x)*profile.s[0]+(y-profile.y)*profile.s[1]
            yp = (x-profile.x)*profile.n[0]+(y-profile.y)*profile.n[1]
            index = np.nonzero((xp>profile.xpmax)|(xp<profile.xpmin)|(yp>profile.ypmax)|(yp<profile.ypmin))
            self.x,self.y,self.z = np.delete(x,index),np.delete(y,index),np.delete(z,index)*self.scale
            # print(np.nanmin(self.z),np.nanmax(self.z))
            # sys.exit()
        except Exception as e: 
            logger.critical(e)
            print(topo.__doc__)
            sys.exit()

class seismi:
    """ 
    seismi class: Load seimicity file for plot 
    Parameters: 
    filename: name input file
    name: name file for plot
    wdir: path input file
    scale: scale values
    color, width: plot options
    """

    def __init__(self,name,wdir,filename,width, color='orange',scale=1):
        self.name = name
        self.wdir = wdir
        self.filename = filename
        self.color = color
        self.width = width
        self.scale=scale

    def load(self,flt):
        logger = flt.logger
        try:
            fmodel = flt.fmodel
            profile = flt.profile
            fname = file(self.wdir+self.filename)
            x,y,z,mw = np.loadtxt(fname,comments = '#',unpack = True,dtype = 'f,f,f,f')
            xp = (x-profile.x)*profile.s[0]+(y-profile.y)*profile.s[1]
            yp = (x-profile.x)*profile.n[0]+(y-profile.y)*profile.n[1]
            index = np.nonzero((xp>profile.xpmax)|(xp<profile.xpmin)|(yp>profile.ypmax)|(yp<profile.ypmin))
            self.x,self.y,self.z,self.mw = np.delete(x,index),np.delete(y,index),np.delete(z,index)*self.scale,np.delete(mw,index)
        except Exception as e: 
            logger.critical(e)
            print(seismi.__doc__)
            sys.exit()

class moho:
    """ 
    moho class: Load moho file for plot 
    Parameters: 
    filename: name input file
    name: name file for plot
    wdir: path input file
    scale: scale values
    color, width: plot options
    """

    def __init__(self,name,wdir,filename,width,color='red',scale=1):
        self.name = name
        self.wdir = wdir
        self.filename = filename
        self.color = color
        self.width = width
        self.scale=scale
        
    def load(self,flt):
        logger = flt.logger
        try:
            fmodel = flt.fmodel
            profile = flt.profile
            fname = file(self.wdir+self.filename)
            x,y,z = np.loadtxt(fname,comments = '#',unpack = True,dtype = 'f,f,f')
            xp = (x-profile.x)*profile.s[0]+(y-profile.y)*profile.s[1]
            yp = (x-profile.x)*profile.n[0]+(y-profile.y)*profile.n[1]
            index = np.nonzero((xp>profile.xpmax)|(xp<profile.xpmin)|(yp>profile.ypmax)|(yp<profile.ypmin))
            self.x,self.y,self.z = np.delete(x,index),np.delete(y,index),np.delete(z,index)*self.scale
        except Exception as e: 
            logger.critical(e)
            print(moho.__doc__)
            sys.exit()

class fault2d:
    """ 
    fault2d class: Load 2D fault for plot only
    help to position fault for futher modeling
    Parameters: 
    name: name fault
    x,y: position east, north
    """

    def __init__(self,name,x,y,strike=None):
        self.name=name
        self.x=x
        self.y=y
        self.strike=strike
