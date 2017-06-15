#!/usr/bin/env python2.7

print 
print '# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #'
print '#                                                                   #'
print '#         Bayesian inversion for fault geometry and slip           #'
print '#                                                                   #'
print '# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #'

from flatten import *
from plot2d import *

import numpy as np
#from numpy.lib.stride_tricks import as_strided
import math

from os import path, environ
from sys import argv,exit,stdin,stdout
# from scipy import optimize
from pylab import *
import getopt
#import multiprocessing
import pymc

def usage():
    print ''
    print 'optimize.py : Nonlinear 2 dimentional inversion of the fault slip and geometry'
    print 'usage: optmize.py input_file'
    print ''

# Read input file
try:
    opts,args = getopt.getopt(sys.argv[1:], "h", ["help"])
except getopt.GetoptError, err:
    print str(err)
    print "for help use --help"
    sys.exit()

for o in opts:
    if o in ("-h","--help"):
       usage()
       sys.exit()
    else:
       assert False, "unhandled option"
       sys.exit()

if 1 == len(argv):
  usage()
  assert False, "no input file"

elif 2 == len(argv):
  fname = sys.argv[1]
  print
  print 'input file :', fname
  sys.path.append(path.dirname(path.abspath(fname)))
  exec ("from " + path.basename(fname)+ " import *")
else:
  assert False, "too many arguments"

# Create directories for output files
inv.outdir = outdir
if not os.path.exists(inv.outdir):
        os.makedirs(inv.outdir)
outgps = inv.outdir+'gps/'
if not os.path.exists(outgps):
        os.makedirs(outgps)
outinsar = inv.outdir+'insar/'
if not os.path.exists(outinsar):
        os.makedirs(outinsar)
outmap = inv.outdir+'map/'
if not os.path.exists(outmap):
        os.makedirs(outmap)
outpro = inv.outdir+'profiles/'
if not os.path.exists(outpro):
        os.makedirs(outpro)
outstat = inv.outdir+'stat/'
if not os.path.exists(outstat):
        os.makedirs(outstat)

print
print "---------------------------------------------------------------------------"
print '   PRIOR MODEL  '
print "---------------------------------------------------------------------------"
print

# data
inv.gpsdata, inv.insardata = gpsdata, insardata
manifolds = inv.gpsdata+inv.insardata
manifolds = flatten(manifolds)

# fault segment model
fmodel = []
fmodel.append(map((lambda x: getattr(x,'segments')),inv.structures))
inv.fmodel = flatten(fmodel)

# profile parameters
profile = inv.profiles
profile.xp0 = (profile.x-inv.fmodel[0].x)*inv.s[0]+(profile.y-inv.fmodel[0].y)*inv.s[1]
profile.yp0 = (profile.x-inv.fmodel[0].x)*inv.n[0]+(profile.y-inv.fmodel[0].y)*inv.n[1]
profile.ypmax,profile.ypmin = profile.yp0+profile.l/2,profile.yp0-profile.l/2
profile.xpmax,profile.xpmin = profile.xp0+profile.w/2,profile.xp0-profile.w/2
    
# Load data plot
inv.gmtfiles = gmtfiles
inv.plotdata = plotdata
inv.topodata = topodata
for i in xrange(len(inv.plotdata)):
    plot = inv.plotdata[i]
    plot.load(inv)
            
    # perp and par composante ref to the first fault 
    plot.yp = (plot.x-inv.fmodel[0].x)*inv.n[0]+(plot.y-inv.fmodel[0].y)*inv.n[1]
    plot.xp = (plot.x-inv.fmodel[0].x)*inv.s[0]+(plot.y-inv.fmodel[0].y)*inv.s[1]

for i in xrange(len(inv.topodata)):
    plot = inv.topodata[i]
    plot.load(inv)
            
    # perp and par composante ref to the first fault 
    plot.yp = (plot.x-inv.fmodel[0].x)*inv.n[0]+(plot.y-inv.fmodel[0].y)*inv.n[1]
    plot.xp = (plot.x-inv.fmodel[0].x)*inv.s[0]+(plot.y-inv.fmodel[0].y)*inv.s[1]
    # print plot.z
    # sys.exit()


# get size of the green function 
# number of structures
inv.Mstruc = len(inv.structures)
# number of fault segments
inv.Mseg = sum(map((lambda x: getattr(x,'Mseg')),inv.structures))
inv.Mbase = 0
# baseline term for gps data depends on dimension
for j in xrange(len(inv.gpsdata)):
    inv.Mbase = inv.Mbase+inv.gpsdata[j].dim
# baseline term for insar equal 2 (linear ramp)
inv.Mbase = inv.Mbase+len(inv.insardata)*2
# volumic deformation
inv.volum = volumic
inv.Mvol =  len(volumic)

print 'Number of structures: ' , inv.Mstruc
print 'Number of segments: ' , inv.Mseg
print 'Number of data network: ', len(manifolds)

# get number of free parameters
inv.Mdis = sum(map((lambda x: getattr(x,'Mker')),inv.structures))
inv.Minv = inv.Mdis+inv.Mvol*2
inv.M = inv.Minv+inv.Mbase
print 'M: ', inv.M

fidf = open(outstat+'{}.txt'.format(fname), 'w')
fidf.write('Prior Model:\n')
print
print 'Prior model:'
for j in xrange(inv.Mseg):
    inv.fmodel[j].info()
    inv.fmodel[j].write(fidf)

# samplinh
inv.nsample=niter-nburn
# Ini traces
for i in xrange(inv.Mseg):
    inv.fmodel[i].traceD = inv.fmodel[i].D*np.ones((inv.nsample))
    inv.fmodel[i].traceH = inv.fmodel[i].H*np.ones((inv.nsample))
    inv.fmodel[i].tracew = inv.fmodel[i].w*np.ones((inv.nsample))

def buildm():
    m_name,m_init,sigmam,pdist = np.zeros(inv.M).tolist(),np.zeros(inv.M),np.zeros(inv.M),['Unif'] * inv.M
    # initial fault parameter
    start = 3
    # get parameters for secondary segments
    for j in xrange(1,inv.Mseg):
        m_init[start],sigmam[start],pdist[start] = inv.fmodel[j].ss,inv.fmodel[j].sigmass,inv.fmodel[j].distss
        m_name[start] = '{} Strike Slip'.format(inv.fmodel[j].name)
        m_init[start+1],sigmam[start+1],pdist[start+1] = inv.fmodel[j].D,inv.fmodel[j].sigmaD,inv.fmodel[j].distss
        m_name[start+1] = '{} D'.format(inv.fmodel[j].name)
        m_init[start+2],sigmam[start+2],pdist[start+2] = inv.fmodel[j].H,inv.fmodel[j].sigmaH,inv.fmodel[j].distH
        m_name[start+2] = '{} H'.format(inv.fmodel[j].name)
        start+= 3
    
    # get parameters for the main segment
    m_init[0],sigmam[0],pdist[0] = inv.fmodel[0].ss,inv.fmodel[0].sigmass,inv.fmodel[0].distss
    m_name[0] = '{} Strike Slip'.format(inv.fmodel[0].name)
    m_init[1],sigmam[1],pdist[1] = inv.fmodel[0].ds,inv.fmodel[0].sigmads,inv.fmodel[0].distshort
    m_name[1] = '{} Shortening'.format(inv.fmodel[0].name)
    m_init[2],sigmam[2],pdist[2] = inv.fmodel[0].H,inv.fmodel[0].sigmaH,inv.fmodel[0].distH
    m_name[2] = '{} H'.format(inv.fmodel[0].name)

    # get parameter for volumic deformations
    for j in xrange(inv.Mvol):
        print m_init[inv.Mdis+j]
        m_init[inv.Mdis+j], sigmam[inv.Mdis+j], pdist[inv.Mdis+j] =inv.volum[j].ds, inv.volum[j].sigmads, inv.fmodel[j].distshort
        m_name[inv.Mdis+j] = '{} Shortening'.format(inv.volum[j].name)
        m_init[inv.Mdis+j+1], sigmam[inv.Mdis+j+1], pdist[inv.Mdis+j+1] =inv.volum[j].D, inv.volum[j].sigmaD, inv.fmodel[j].distD
        m_name[inv.Mdis+j+1] = '{} D'.format(inv.volum[j].name)
    
    # print m_init
    # sys.exit()

    # !!! hardcoding for baseline uncertainties !!!
    M = 0 
    for i in xrange(len(manifolds)):
        # define baseline network
        if manifolds[i].dim>1:
            m_init[inv.Minv+M:inv.Minv+M+manifolds[i].dim],sigmam[inv.Minv+M:inv.Minv+M+manifolds[i].dim] = 0, 50
            M += manifolds[i].dim
        else:
            m_init[inv.Minv+M],sigmam[inv.Minv+M] = 0., 10.
            m_init[inv.Minv+M+1],sigmam[inv.Minv+M+1] = 0, 0.1
            M += 2

    #print sigmam
    #sys.exit()
    #m_init[inv.Minv:],sigmam[inv.Minv:] = 0,50
    uu = 1
    for i in range(inv.Minv, inv.M):
        m_name[i] = '{} Baseline {}'.format(profile.name,uu)
        uu +=  1
    
    # Define the bounds
    m_max = (m_init+sigmam)
    m_min = (m_init-sigmam)

    # export m_init and sigmam
    for i in xrange(len(manifolds)):
        manifolds[i].sigmam = sigmam
        manifolds[i].m_init = m_init
            
    b = np.column_stack((m_min,m_max))
    return m_name,m_init.tolist(),sigmam.tolist(),m_min.tolist(),m_max.tolist(),pdist
# get data 
def data():
    d = np.zeros((N))
    start = 0
    for i in xrange(len(manifolds)):
        d[start:start+manifolds[i].N] = manifolds[i].d
        start+= manifolds[i].N
    return d

m_name,m_init,sigmam,m_min,m_max,pdist = buildm()

# manifolds: targets 
for i in xrange(len(manifolds)):
    manifolds[i].load(inv)
N = sum(map((lambda x: getattr(x,'N')),manifolds))
print 'N = ', N

for i in xrange(len(gpsdata)):
    gpsdata[i].computeField()
#sys.exit()

for item in manifolds:
    item.info()

# Create the a priori functions
inv.Priors = []
inv.Sampled = []
inv.Fixed = []
inv.mu, inv.tau = [], []
inv.mmin, inv.mmax = [], []

for name, lower, upper, initial, sigma, dist in zip(m_name, m_min, m_max, m_init, sigmam, pdist):
    # print dist
    sigma = sigma/2 # 1-sigma: sigma/2 ... !!! bidouille !!! 
    if upper!= lower:
        # builing model
        if dist is 'Gaus':
            # sys.exit(0) 
            sigma = sigma/2
            p = pymc.Normal(name, mu = initial, tau = 1./sigma**2, value = initial)
        elif dist is 'Logn':
            # colors = ["#348ABD", "#A60628", "#7A68A6"]
            # m1 = [pymc.rnormal(initial,1./sigma**2) for i in range(20000)]
            # m2 = [pymc.rnormal(initial,1./sigma) for i in range(20000)]
            # #m3 = [pymc.rnormal(initial,sigma) for i in range(20000)]
            # plt.hist(m,bins=50,histtype='stepfilled',normed=True,label='{},sigma={}'.format(name,sigma),color="#7A68A6")
            # plt.legend()
            # plt.show()
            # sys.exit()
            frac = (sigma*2)/5
            p = pymc.Lognormal(name, mu = math.log(initial), tau = 1./(0.125*frac)**2)
        else:
            p = pymc.Uniform(name, lower, upper, value = initial)
        
        inv.Sampled.append(name)
        inv.Priors.append(p)
        # save parameters for plots
        inv.mu.append(initial)
        inv.tau.append(sigma)
        inv.mmin.append(lower)
        inv.mmax.append(upper)
    
    # Not sampled if sigma is null
    elif upper == lower:
        inv.Fixed.append(name)  
    
    else:
        print('Probleme dans la definition du prior pour le parametre {}'.format(name))
        sys.exit(1)

@pymc.deterministic(plot = False)
def forward(theta = inv.Priors):
    msample = theta
    g = np.zeros((N))
    start = 0
    M = 0
    # Rebuild the full m vector
    m = []
    uu = 0
    for name, initial in zip(m_name, m_init):
        if name in inv.Sampled:
            m.append(msample[uu])
            uu +=  1
        elif name in inv.Fixed:
            m.append(initial)
    for i in xrange(len(manifolds)):
        # define baseline network
        if 3==manifolds[i].dim:
                manifolds[i].a,manifolds[i].b,manifolds[i].c = m[inv.Minv+M],m[inv.Minv+M+1],m[inv.Minv+M+2]
                M+= 3
        else:
                manifolds[i].a,manifolds[i].b = m[inv.Minv+M],m[inv.Minv+M+1]
                M+= 2
        
        # build g
        g[start:start+manifolds[i].N] = manifolds[i].g(m)
        start+= manifolds[i].N
    return g

## full covariance matrix
if fullcov is 'yes':
    def Cov():
        Cov = np.zeros((N,N))
        start = 0
        for i in xrange(len(manifolds)):
            if manifolds[i].dim == 1:
                Cd = manifolds[i].computeFullCovariance(xbounds = [None, None],
                                                    plot = True,
                                                    distmax = 50.,
                                                    every = 1.,
                                                    outdir = inv.outdir,
                                                    ramp=rampcov,
                                                    maskcov=maskcov
                                                    )
                # weight insar data
                Cd *= manifolds[i].wd**2
                # change sigmad for plot
                manifolds[i].sigmad = np.diag(Cd)
            else:
                Cd = np.diag(manifolds[i].sigmad**2,k = 0)
            manifolds[i].Cd = Cd
            manifolds[i].invCd = np.linalg.inv(Cd)
            Cov[start:start+manifolds[i].N,start:start+manifolds[i].N] = Cd
            start+= manifolds[i].N

        fig = plt.figure(2)
        ax = fig.add_subplot(111)
        cax = plt.imshow(Cov)
        fig.colorbar(cax)
        fig.savefig(outdir+'/insar/'+'COV_mat.eps', format='EPS')
        plt.show()
        return Cov
    d = pymc.MvNormalCov('Data', mu = forward, C = Cov(), value = data(), observed = True)

# autocovariance only
else:
    def Cov():
        Cov = np.zeros((N))
        start = 0
        for i in xrange(len(manifolds)):
            Cd = np.diag(manifolds[i].sigmad**2,k = 0)
            # And the diagonal of its inverse
            Cov[start:start+manifolds[i].N] = np.diag(np.linalg.inv(Cd))
            # Save Covariance matrix for each data set
            manifolds[i].Cd = np.diag(Cd)
            manifolds[i].invCd = np.diag(np.linalg.inv(Cd))
            start+= manifolds[i].N
        return Cov
    d = pymc.Normal('Data', mu = forward, tau = Cov(), value = data(), observed = True) # observed = True for constant stochastic values


print
print "---------------------------------------------------------------------------"
print "    BEST-FIT solution .....    "
print "---------------------------------------------------------------------------"
print


Parameters = pymc.Model( inv.Priors + [d] )
# Create a sampler
map_ = pymc.MAP(Parameters)
# map_.fit() # best-fit solution
# print 'Akaike information criterion for the model: ', map_.AIC
# print 'The Bayesian information criterion for the model: ', map_.BIC
# for variable in map_.variables:
#     print {str(variable): variable.value}

#map_.coef_
model = pymc.MCMC(Parameters)
for p,sigma in zip(inv.Priors,inv.tau):
    model.use_step_method(pymc.Metropolis, p)
    ### tests optimisation parameters...pymc black box
    #model.use_step_method(pymc.AdaptiveMetropolis, p, shrink_if_necessary=False, interval=5000)
    #model.use_step_method(pymc.Metropolis, p, proposal_distribution='Prior')
    #model.use_step_method(pymc.Metropolis, p, proposal_sd=sigma/np.sqrt(N), proposal_distribution='Normal')
    #model.use_step_method(pymc.AdaptiveMetropolis, p, delay=nburn, greedy=True, shrink_if_necessary=True)

print
print "---------------------------------------------------------------------------"
print "    SAMPLING .....    "
print "---------------------------------------------------------------------------"
print

model.sample(iter = niter, burn = nburn, thin=1)

# check acceptance rates
#print
#print 'Checking accpetance rates:'
#for p,name in zip(inv.Priors,inv.Sampled):
#    print '{} acceptance rates: {}'.format(name, model.step_method_dict[p][0].ratio)

print
print "---------------------------------------------------------------------------"
print ' POSTERIOR MODEL   '
print "---------------------------------------------------------------------------"
print


print
print "# Average model: "

# non-graphical summary of the posterior 
print 'Writing output files...'
#model.stats()

# Compute 95% HDI (Highest Density Interval) is the smallest width interval to contain 95% of the posterior probability. 
def hdi(trace, cred_mass = 0.95):
    hdi_min, hdi_max = pymc.utils.calc_min_interval(np.sort(trace), 1.0-cred_mass)
    return hdi_min, hdi_max

opts = {'c':'green', 'linestyle':'--'}
nfigures = 100

# attention a bien respecter ordre remplissage mf...
mf = []
inv.traces = []
inv.labels = []

if '{} Strike Slip'.format(inv.fmodel[0].name) in inv.Sampled:
    m = model.trace('{} Strike Slip'.format(inv.fmodel[0].name))[:]
    mf.append(np.mean(m))
    inv.fmodel[0].ss = np.mean(m)
    inv.fmodel[0].sigmass = 2*np.std(m) 
    inv.traces.append(m)
    inv.labels.append('{} Strike Slip'.format(inv.fmodel[0].name))

    fid = open(outstat+'{}_ss.txt'.format(inv.fmodel[0].name), 'w')
    np.savetxt(fid, m)
    fid.close()
else:
    inv.fmodel[0].ss = m_init[0]
   
if '{} Shortening'.format(inv.fmodel[0].name) in inv.Sampled:
    m = model.trace('{} Shortening'.format(inv.fmodel[0].name))[:]
    mf.append(np.mean(m))
    inv.fmodel[0].ds = np.mean(m)
    inv.fmodel[0].sigmads = 2*np.std(m) 
    inv.traces.append(m)
    inv.labels.append('{} Shortening'.format(inv.fmodel[0].name))

    fid = open(outstat+'{}_ds.txt'.format(inv.fmodel[0].name), 'w')
    m = model.trace('{} Shortening'.format(inv.fmodel[0].name))[:]
    np.savetxt(fid, m)
    fid.close()

if '{} H'.format(inv.fmodel[0].name) in inv.Sampled:
    m = model.trace('{} H'.format(inv.fmodel[0].name))[:]
    mf.append(np.mean(m))
    inv.fmodel[0].H,inv.fmodel[0].w = np.mean(m),np.mean(m)
    inv.fmodel[0].sigmaw,inv.fmodel[0].sigmaH = 2*np.std(m), 2*np.std(m) 
    inv.fmodel[0].tracew = np.asarray(m)
    inv.traces.append(m)
    inv.labels.append('{} H'.format(inv.fmodel[0].name))

    fid = open(outstat+'{}_H.txt'.format(inv.fmodel[0].name), 'w')
    np.savetxt(fid, m)
    fid.close()

for j in xrange(1,inv.Mseg):
    if '{} Strike Slip'.format(inv.fmodel[j].name) in inv.Sampled:
        m = model.trace('{} Strike Slip'.format(inv.fmodel[j].name))[:]
        mf.append(np.mean(m))
        inv.fmodel[j].ss = np.mean(m)
        inv.fmodel[j].sigmass = 2*np.std(m)
        inv.traces.append(m)
        inv.labels.append('{} Strike Slip'.format(inv.fmodel[j].name))

        fid = open(outstat+'{}_ss.txt'.format(inv.fmodel[j].name), 'w')
        np.savetxt(fid, m)
        fid.close()
    
    if '{} D'.format(inv.fmodel[j].name) in inv.Sampled:
        m = model.trace('{} D'.format(inv.fmodel[j].name))[:]
        mf.append(np.mean(m))
        inv.fmodel[j].D = np.mean(m)
        inv.fmodel[j].sigmaD = 2*np.std(m) 
        inv.fmodel[j].traceD = np.asarray(m)
        inv.traces.append(m)
        inv.labels.append('{} D'.format(inv.fmodel[j].name))

        fid = open(outstat+'{}_L.txt'.format(inv.fmodel[j].name), 'w')
        np.savetxt(fid, m)
 
    if '{} H'.format(inv.fmodel[j].name) in inv.Sampled:
        m = model.trace('{} H'.format(inv.fmodel[j].name))[:]
        mf.append(np.mean(m))
        inv.fmodel[j].H=np.mean(m)
        inv.fmodel[j].sigmaw,inv.fmodel[j].sigmaH = 2*np.std(m),2*np.std(m) 
        inv.fmodel[j].traceH = np.asarray(m)
        inv.traces.append(m)
        inv.labels.append('{} H'.format(inv.fmodel[j].name))

        fid = open(outstat+'{}_H.txt'.format(inv.fmodel[j].name), 'w')
        np.savetxt(fid, m)
        fid.close()
        fid.close()

# Compute ss main fault
#tot_ss = sum([inv.fmodel[j].ss for j in xrange(1,inv.Mseg)])
tot_ss = 0
for j in xrange(1,inv.Mseg):
        tot_ss += inv.fmodel[j].ss
inv.fmodel[0].sst = inv.fmodel[0].ss + tot_ss 
# print inv.fmodel[0].sst, tot_ss, inv.fmodel[0].ss
# sys.exit()

# save all plaussible models for plot
inv.fmodel[0].traceF = inv.fmodel[0].D*np.ones((inv.nsample))
if inv.structures[0].Mseg >1:
    inv.fmodel[1].tracew,inv.fmodel[2].tracew = inv.fmodel[0].tracew - inv.fmodel[1].traceH, inv.fmodel[0].tracew - inv.fmodel[2].traceH
    inv.fmodel[1].traceF,inv.fmodel[2].traceF = inv.fmodel[1].traceD,inv.fmodel[2].traceD
Mtemp = inv.structures[0].Mseg
for j in xrange(1,inv.Mstruc):
    for k in xrange(inv.structures[j].Mseg):
        if abs(inv.fmodel[Mtemp+k].fperp) > 1:
            inv.fmodel[Mtemp+k].tracew = inv.fmodel[Mtemp-1].tracew - inv.fmodel[Mtemp+k].H
            inv.fmodel[Mtemp+k].traceF = inv.fmodel[Mtemp-1].traceF + inv.fmodel[Mtemp+k].D
        else:
            inv.fmodel[Mtemp+k].tracew = inv.fmodel[0].tracew - inv.fmodel[Mtemp+k].H
            inv.fmodel[Mtemp+k].traceF = inv.fmodel[Mtemp+k].D*np.ones((inv.nsample))
    Mtemp += inv.structures[j].Mseg

for j in xrange(0,inv.Mvol):
    if '{} Shortening'.format(inv.volum[j].name) in inv.Sampled:
        m = model.trace('{} Shortening'.format(inv.volum[j].name))[:]
        mf.append(np.mean(m))
        inv.volum[j].ds = np.mean(m)
        inv.volum[j].sigmads = 2*np.std(m)
    
    if '{} D'.format(inv.volum[j].name) in inv.Sampled:
        m = model.trace('{} D'.format(inv.volum[j].name))[:]
        mf.append(np.mean(m))
        inv.volum[j].D = np.mean(m)
        inv.volum[j].sigmaD = 2*np.std(m)

uu = 1
for i in xrange(len(manifolds)):
    if 3 == manifolds[i].dim:
        m = model.trace('{} Baseline {}'.format(profile.name,uu))[:]
        manifolds[i].a = np.mean(m)
        manifolds[i].sigmaa=np.std(m)
        mf.append(manifolds[i].a)

        m = model.trace('{} Baseline {}'.format(profile.name,uu+1))[:]       
        manifolds[i].b = np.mean(m)
        manifolds[i].sigmab=np.std(m)
        mf.append(manifolds[i].b)

        m = model.trace('{} Baseline {}'.format(profile.name,uu+2))[:]
        manifolds[i].c = np.mean(m)
        manifolds[i].sigmac=np.std(m)
        mf.append(manifolds[i].c)
        
        print '%s: %f Upar + %f Uperp %f Uv '%(manifolds[i].reduction, manifolds[i].a,manifolds[i].b,manifolds[i].c)
        fidf.write('{}\n'.format(manifolds[i].reduction))
        fidf.write('a:par\tb:perp\tc:vert\n')
        np.savetxt(fidf, np.vstack([manifolds[i].a,manifolds[i].b,manifolds[i].c]).T, fmt='%.6f',delimiter='\t', newline='\n')
        np.savetxt(fidf, np.vstack([manifolds[i].sigmaa,manifolds[i].sigmab,manifolds[i].sigmac]).T, fmt='%.6f',delimiter='\t', newline='\n')

        uu+= 3
    else:
        m = model.trace('{} Baseline {}'.format(profile.name,uu))[:]
        manifolds[i].a = np.mean(m)
        mf.append(manifolds[i].a)
        manifolds[i].sigmaa=np.std(m)
       
        m = model.trace('{} Baseline {}'.format(profile.name,uu+1))[:]
        manifolds[i].b = np.mean(m)
        mf.append(manifolds[i].b)
        manifolds[i].sigmab=np.std(m)
        
        if 2 == manifolds[i].dim:
            print '%s: %f Upar + %f Uperp'%(manifolds[i].reduction, manifolds[i].a,manifolds[i].b)

            fidf.write('{}\n'.format(manifolds[i].reduction))
            fidf.write('a:par\tb:perp\n')
            np.savetxt(fidf, np.vstack([manifolds[i].a,manifolds[i].b]).T,fmt='%.6f',delimiter='\t', newline='\n' )
            np.savetxt(fidf, np.vstack([manifolds[i].sigmaa,manifolds[i].sigmab]).T,fmt='%.6f',delimiter='\t', newline='\n' )

        else:
            print 'a + b yperp %s: %f + %f yperp '%(manifolds[i].reduction, manifolds[i].a,manifolds[i].b)
            fidf.write('{}\n'.format(manifolds[i].reduction))
            fidf.write('a:cst\tb:lin\n')
            np.savetxt(fidf, np.vstack([manifolds[i].a,manifolds[i].b]).T, fmt='%.6f' , delimiter='\t',newline='\n')
            np.savetxt(fidf, np.vstack([manifolds[i].sigmaa,manifolds[i].sigmab]).T, fmt='%.6f' , delimiter='\t',newline='\n')

        uu+= 2

# mean forward model 
uu = 0
m = []
for name, initial in zip(m_name, m_init):
  if name in inv.Sampled:
    m.append(mf[uu])
    uu += 1
  elif name in inv.Fixed:
    m.append(initial)
print m

start = 0
if fullcov is 'yes':
    g=np.zeros((N))
    r=np.zeros((N,N))
    for i in xrange(len(manifolds)):
        g[start:start+manifolds[i].N]=manifolds[i].g(m)
        r[start:start+manifolds[i].N,start:start+manifolds[i].N]=manifolds[i].residual(m,np.linalg.det(manifolds[i].Cd),manifolds[i].invCd) 
        start+=manifolds[i].N
else:    
    g=np.zeros((N))
    r=np.zeros((N))
    for i in xrange(len(manifolds)):
        g[start:start+manifolds[i].N]=manifolds[i].g(m)
        r[start:start+manifolds[i].N]=manifolds[i].residual(m,manifolds[i].invCd,manifolds[i].invCd) 
        start+=manifolds[i].N

print
print '-----------------------'
print 'likelihoods of the residus: ', np.sum(pow(r,2))
print '-----------------------'
print

print
print '--------------------------------------------'
print 'Checking horizontal conservation of motion of mean model:'
if inv.structures[0].Mseg >1:
    print 'flower structure: ', 1
    print '{} = {} - {}'.format(inv.fmodel[0].name, inv.fmodel[1].name, inv.fmodel[2].name)  
    print '{} = {} - {}'.format(inv.fmodel[0].vh, inv.fmodel[1].vh, inv.fmodel[2].vh)
else:
    print inv.fmodel[0].name, inv.fmodel[0].vh
Mtemp = inv.structures[0].Mseg 
for j in xrange(1,inv.Mstruc):
    if inv.structures[j].Mseg > 1:
        print 'flower structure: ', j
        print '{} = {} - {}'.format(inv.fmodel[Mtemp-1].name,inv.fmodel[Mtemp].name, inv.fmodel[Mtemp+1].name) 
        print '{} = {} - {}'.format(inv.fmodel[Mtemp-1].vh,inv.fmodel[Mtemp].vh, inv.fmodel[Mtemp+1].vh)
    else:
        print inv.fmodel[Mtemp].name, inv.fmodel[Mtemp].vh

    Mtemp+=inv.structures[j].Mseg
print '--------------------------------------------'
print

fig2 = plt.figure(1000)
ax = fig2.add_subplot(111)
colors=['blue','m','orange','yellow','green','black']
for i in xrange(len(manifolds)):
    if 3==manifolds[i].dim:
        ax.plot(manifolds[i].d[::manifolds[i].dim]-manifolds[i].a,manifolds[i].dp[::manifolds[i].dim]-manifolds[i].a,'+',color='blue',label='GPS profile-perpendicular')
        ax.errorbar(manifolds[i].d[::manifolds[i].dim]-manifolds[i].a,manifolds[i].dp[::manifolds[i].dim]-manifolds[i].a,xerr=2*manifolds[i].sigmad[::manifolds[i].dim],ecolor='blue',barsabove='True',fmt=None)
        ax.plot(manifolds[i].d[1::manifolds[i].dim]-manifolds[i].b,manifolds[i].dp[1::manifolds[i].dim]-manifolds[i].b,'+',color='green',label='GPS profile-parallel')
        ax.errorbar(manifolds[i].d[1::manifolds[i].dim]-manifolds[i].b,manifolds[i].dp[1::manifolds[i].dim]-manifolds[i].b,xerr=2*manifolds[i].sigmad[1::manifolds[i].dim],ecolor='green',barsabove='True',fmt=None)
        ax.plot(manifolds[i].d[2::manifolds[i].dim]-manifolds[i].c,manifolds[i].dp[2::manifolds[i].dim]-manifolds[i].c,'+',color='red',label='GPS vertical')
        ax.errorbar(manifolds[i].d[2::manifolds[i].dim]-manifolds[i].c,manifolds[i].dp[2::manifolds[i].dim]-manifolds[i].c,xerr=2*manifolds[i].sigmad[2::manifolds[i].dim],ecolor='red',barsabove='True',fmt=None)
    elif 2==manifolds[i].dim:
        ax.plot(manifolds[i].d[::manifolds[i].dim]-manifolds[i].a,manifolds[i].dp[::manifolds[i].dim]-manifolds[i].a,'x',color='blue',label='GPS profile-perpendicular')
        ax.errorbar(manifolds[i].d[::manifolds[i].dim]-manifolds[i].a,manifolds[i].dp[::manifolds[i].dim]-manifolds[i].a,xerr=2*manifolds[i].sigmad[::manifolds[i].dim],ecolor='blue',barsabove='True',fmt=None)
        ax.plot(manifolds[i].d[1::manifolds[i].dim]-manifolds[i].b,manifolds[i].dp[1::manifolds[i].dim]-manifolds[i].b,'x',color='green',label='GPS profile-parallel')
        ax.errorbar(manifolds[i].d[1::manifolds[i].dim]-manifolds[i].b,manifolds[i].dp[1::manifolds[i].dim]-manifolds[i].b,xerr=2*manifolds[i].sigmad[1::manifolds[i].dim],ecolor='green',barsabove='True',fmt=None)
    else:
        ax.plot(manifolds[i].d-(manifolds[i].a+manifolds[i].b*manifolds[i].yp),manifolds[i].dp[::manifolds[i].dim]-(manifolds[i].a+manifolds[i].b*manifolds[i].yp),'.',color=manifolds[i].color,label=manifolds[i].reduction)
        #ax.errorbar(manifolds[i].d-(manifolds[i].a+manifolds[i].b*manifolds[i].yp),manifolds[i].dp[::manifolds[i].dim]-(manifolds[i].a+manifolds[i].b*manifolds[i].yp),xerr=manifolds[i].sigmad,ecolor='red',barsabove='True',fmt=None)

plt.xlabel=('d')
plt.ylabel=('g(mf)')
title('d=f(g(mf))')
ax.legend(loc=2)
ax=plt.gca()
ax.grid(True)
xmin,xmax=ax.get_xlim()
plt.plot([xmin,xmax],[xmin,xmax])
fig2.savefig(outdir+profile.name+'foward.eps', format='EPS')

plotLOS(inv,nfigures)
nfigures +=  1
plotMap(inv,x1,y1,nfigures)
nfigures +=  1
plotHist(inv,model,nfigures)
nfigures +=  1
# pymc plot function
#pymc.Matplot.plot(model,format = 'eps',path = outstat)

# display amd save results 
for j in xrange(inv.structures[0].Mseg):
    inv.fmodel[j].info()
    inv.fmodel[j].write(fidf)
Mtemp = inv.structures[0].Mseg
for j in xrange(1,inv.Mstruc):
    for k in xrange(inv.structures[j].Mseg):
        # calcul L, dip, ds components
        inv.fmodel[Mtemp+k].info()
        inv.fmodel[Mtemp+k].write(fidf)
    Mtemp += inv.structures[j].Mseg
for j in xrange(inv.Mvol):
    inv.volum[j].info()
fidf.close()

plt.show()
sys.exit()

