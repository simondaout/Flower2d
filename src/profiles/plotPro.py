#!/usr/bin/env python2.7 

from __future__ import print_function
import numpy as np
from readgmt import *
from matplotlib import pyplot as plt
#from matplotlib import mpl
import matplotlib
import matplotlib.cm as cm

from sys import argv,exit,stdin,stdout
import getopt
import os
from os import path

import logging

def usage():
  print('plotPro.py infile.py [-v] [-h]')
  print('-v Verbose mode. Show more information about the processing')
  print('-h Show this screen')

#load input file 
try:
    opts,args = getopt.getopt(sys.argv[1:], "h", ["help"])
except getopt.GetoptError, err:
    print(str(err))
    print("for help use --help")
    sys.exit()

level = 'basic'
for o in sys.argv:
    if o in ("-h","--help"):
       usage()
       sys.exit()
    if o in ("-v","--verbose"):
      level = 'debug'

# init logger 
if level == 'debug':
    logging.basicConfig(level=logging.DEBUG,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
else:
    logging.basicConfig(level=logging.INFO,\
        format='%(asctime)s -- %(levelname)s -- %(message)s')
logger = logging.getLogger('plotPro.log')

if 1==len(sys.argv):
  usage()
  assert False, "no input file"
  logger.critical('No input file')
  sys.exit()

if len(sys.argv)>1:
  fname=sys.argv[1]
  logger.info('Read input file {0}'.format(fname))
  try:
    sys.path.append(path.dirname(path.abspath(fname)))
    exec ("from "+path.basename(fname)+" import *")
  except:
    execfile(path.abspath(fname))

if not os.path.exists(outdir):
    logger.info('Creating output directory {0}'.format(outdir))
    os.makedirs(outdir)

# distance between fault and center of profile in fault azimuth 
# Fault model
Mfault=len(fmodel)
fperp=np.zeros(Mfault)

# Load data
for i in xrange(len(topodata)):
    plot=topodata[i]
    logger.info('Load data {0}'.format(plot.name))
    plot.load()

for i in xrange(len(gpsdata)):
    gps = gpsdata[i]
    logger.info('Load data {0}'.format(gps.network))
    gps.loadgps()

for i in xrange(len(insardata)):
    insar = insardata[i]
    logger.info('Load data {0}'.format(insar.network))
    insar.loadinsar()
    if insar.theta is True:
      logger.debug('Read incidence angle...')
      insar.losm = np.mean(insar.los)
      # insar.ulos = insar.ulos * \
      #   (np.sin(np.deg2rad(insar.losm))/np.sin(np.deg2rad(insar.los)))

# MAP
fig=plt.figure(0,figsize = (9,8))
ax = fig.add_subplot(1,1,1)
ax.axis('equal')

logger.info('Plot Map ....') 

for i in xrange(len(insardata)):
  
  insar=insardata[i]
  samp = insar.samp

  logger.info('Plot data {0} between {1} and {2}'.format(insar.network, insar.lmin, insar.lmax))
  logger.info('Subsample data every {0} point'.format(insar.samp))
  norm = matplotlib.colors.Normalize(vmin=insar.lmin, vmax=insar.lmax)
  m = cm.ScalarMappable(norm = norm, cmap = 'rainbow')
  m.set_array(insar.ulos[::samp])
  facelos = m.to_rgba(insar.ulos[::samp])
  ax.scatter(insar.x[::samp],insar.y[::samp], s=1, marker = 'o',color = facelos, label = 'LOS Velocity %s'%(insar.reduction))

  # save flatten map
  # np.savetxt('{}_flat'.format(insardata[i].network), np.vstack([insar.x,insar.y,insar.ulos]).T, fmt='%.6f')

for i in xrange(len(gpsdata)):
  gps=gpsdata[i]
  logger.info('Plot GPS data {0}'.format(gps.network))
  ax.quiver(gps.x,gps.y,gps.ux,gps.uy,scale = 100, width = 0.005, color = 'red')

# plot faults
for kk in xrange(Mfault):
  xf,yf = np.zeros((2)),np.zeros((2))
  strike=fmodel[kk].strike
  str=(strike*math.pi)/180
  s=[math.sin(str),math.cos(str),0]
  n=[math.cos(str),-math.sin(str),0]
  xf[0] = fmodel[kk].x+2*-150*s[0]
  xf[1] = fmodel[kk].x+2*150*s[0]
  yf[0] = fmodel[kk].y+2*-150*s[1]
  yf[1] = fmodel[kk].y+2*150*s[1]
  # plot fault
  ax.plot(xf[:],yf[:],'--',color = 'black',lw = 1.)

for ii in xrange(len(gmtfiles)):
  name = gmtfiles[ii].name
  wdir = gmtfiles[ii].wdir
  filename = gmtfiles[ii].filename
  color = gmtfiles[ii].color
  width = gmtfiles[ii].width
  fx,fy = gmtfiles[ii].load()
  for i in xrange(len(fx)):
    ax.plot(fx[i],fy[i],color = color,lw = width)

if 'xmin' in locals(): 
  logger.info('Found boundaries map plot {0}-{1} and {2}-{3} in locals'.format(xmin,xmax,ymin,ymax))
  ax.set_xlim(xmin,xmax)
  ax.set_ylim(ymin,ymax)

# add colorbar los
if 'm' in locals():
  fig.colorbar(m,shrink = 0.5, aspect = 5)

# fig pro topo
if len(profiles) > 2:
  fig1=plt.figure(1,figsize=(7,8))
else:
  fig1=plt.figure(1,figsize=(10,3))
fig1.subplots_adjust(hspace=0.0001)

# fig pro insar
if len(profiles) > 2:
  fig2=plt.figure(4,figsize=(7,8))
else:
  fig2=plt.figure(4,figsize=(10,3))
fig2.subplots_adjust(hspace=0.0001)

# fig pro gps
if len(profiles) > 2:
  fig3=plt.figure(5,figsize=(7,8))
else:
  fig3=plt.figure(5,figsize=(10,3))
fig3.subplots_adjust(hspace=0.0001)

logger.info('Plot Profiles ....')

# Plot profile
for k in xrange(len(profiles)): 

  l=profiles[k].l
  w=profiles[k].w
  x0=profiles[k].x
  y0=profiles[k].y
  strike = profiles[k].strike
  name=profiles[k].name
  typ=profiles[k].typ

  logger.info('Plot profile {0}. length: {1}, width :{2}, strike: {3}'.format(name, l, w, strike)) 

  # lim profile
  ypmax,ypmin=l/2,-l/2
  xpmax,xpmin=w/2,-w/2

  # profile azimuth
  profiles[k].str=(profiles[k].strike*math.pi)/180
  profiles[k].s=[math.sin(profiles[k].str),math.cos(profiles[k].str),0]
  profiles[k].n=[math.cos(profiles[k].str),-math.sin(profiles[k].str),0]

  for j in xrange(Mfault):
    fperp[j]=(fmodel[j].x-profiles[k].x)*profiles[k].n[0]+(fmodel[j].y-profiles[k].y)*profiles[k].n[1]

  ax1=fig1.add_subplot(len(profiles),1,k+1)
  ax1.set_xlim([-l/2,l/2])
  for i in xrange(len(topodata)):
        plot=topodata[i]

        # perp and par composante ref to the profile 
        plot.ypp=(plot.x-profiles[k].x)*profiles[k].n[0]+(plot.y-profiles[k].y)*profiles[k].n[1]
        plot.xpp=(plot.x-profiles[k].x)*profiles[k].s[0]+(plot.y-profiles[k].y)*profiles[k].s[1]

        index=np.nonzero((plot.xpp>xpmax)|(plot.xpp<xpmin)|(plot.ypp>ypmax)|(plot.ypp<ypmin))
        plotxpp,plotypp,plotz=np.delete(plot.xpp,index),np.delete(plot.ypp,index),np.delete(plot.z,index)

            
        nb = np.float(l/(len(plotz)/20.))
        logger.debug('Load {0}. Create bins every {1:.3f} km'.format(plot.name, nb)) 
        bins = np.arange(-l/2,l/2, nb)
        inds = np.digitize(plotypp,bins)
        distance = []
        moy_topo = []
        std_topo = []
        for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            if len(uu)>0:
                distance.append(bins[j] + (bins[j+1] - bins[j])/2.)
                std_topo.append(np.std(plotz[uu]))
                moy_topo.append(np.median(plotz[uu]))
        
        distance = np.array(distance)
        std_topo = np.array(std_topo)
        moy_topo = np.array(moy_topo)

        ax1.plot(distance,-moy_topo,label=plot.name,color='black',lw=1)
        if plot.plotminmax == True:
          logger.debug('plotminmax set to True')
          ax1.plot(distance,-moy_topo-std_topo,color='black',lw=1)
          ax1.plot(distance,-moy_topo+std_topo,color='black',lw=1)
        
        # for kk in xrange(Mfault):    
            # ax1.plot([fperp[kk],fperp[kk]],[8,-8],color='red')
            # ax1.text(fperp[kk],0.5,fmodel[kk].name,color='red')
        
        if (plot.topomin is not None) and (plot.topomax is not None) :
            logger.info('Set ylim to {} and {}'.format(plot.topomin,plot.topomax))
            ax1.set_ylim([plot.topomin,plot.topomax])
  
  # LOS profile/map
  ax2=fig2.add_subplot(len(profiles),1,k+1)
  ax2.set_xlim([-l/2,l/2])

  # LOS profile/map
  ax3=fig3.add_subplot(len(profiles),1,k+1)
  ax3.set_xlim([-l/2,l/2])
  
  # plot profiles
  xp,yp = np.zeros((7)),np.zeros((7))
  xp[:] = x0-w/2*profiles[k].s[0]-l/2*profiles[k].n[0],x0+w/2*\
  profiles[k].s[0]-l/2*profiles[k].n[0],x0+w/2*profiles[k].s[0]+l/2*profiles[k].n[0],x0-w/2*profiles[k].s[0]+l/2*profiles[k].n[0],x0-w/2*profiles[k].s[0]-l/2*profiles[k].n[0],x0-l/2*profiles[k].n[0],x0+l/2*profiles[k].n[0]
  yp[:] = y0-w/2*profiles[k].s[1]-l/2*profiles[k].n[1],y0+w/2*\
  profiles[k].s[1]-l/2*profiles[k].n[1],y0+w/2*profiles[k].s[1]+l/2*profiles[k].n[1],y0-w/2*profiles[k].s[1]+l/2*profiles[k].n[1],y0-w/2*profiles[k].s[1]-l/2*profiles[k].n[1],y0-l/2*profiles[k].n[1],y0+l/2*profiles[k].n[1]

  ax.plot(xp[:],yp[:],color = 'black',lw = 1.)

  # GPS plot
  markers = ['+','d','x','v']
  for i in xrange(len(gpsdata)):
      gps=gpsdata[i]
      gpsmin = gps.lmin
      gpsmax = gps.lmax
      logger.info('Load GPS {0}'.format(gps.network)) 

      # perp and par composante ref to the profile 
      gps.ypp=(gps.x-profiles[k].x)*profiles[k].n[0]+(gps.y-profiles[k].y)*profiles[k].n[1]
      gps.xpp=(gps.x-profiles[k].x)*profiles[k].s[0]+(gps.y-profiles[k].y)*profiles[k].s[1]

      # select data within profile
      index=np.nonzero((gps.xpp>xpmax)|(gps.xpp<xpmin)|(gps.ypp>ypmax)|(gps.ypp<ypmin))
      gpsux,gpsuy,gpssigmax,gpssigmay,gpsx,gpsy,gpsxp,gpsyp=np.delete(gps.ux,index),np.delete(gps.uy,index)\
      ,np.delete(gps.sigmax,index),np.delete(gps.sigmay,index),np.delete(gps.x,index),np.delete(gps.y,index),\
      np.delete(gps.xpp,index),np.delete(gps.ypp,index)

      # co;pute fault parallel and perpendicular for each profiles
      gpsupar = gpsux*profiles[k].s[0]+gpsuy*profiles[k].s[1]
      gpsuperp = gpsux*profiles[k].n[0]+gpsuy*profiles[k].n[1]
      gpssigmaperp=((gpssigmax*np.cos(profiles[k].str))**2 + (gpssigmay*np.sin(profiles[k].str))**2)**0.5
      gpssigmapar=((gpssigmax*np.sin(profiles[k].str))**2 + (gpssigmay*np.cos(profiles[k].str))**2)**0.5

      ax3.plot(gpsyp,gpsupar,markers[i],color = 'blue',mew = 1.,label =\
       '%s fault-parallel velocities'%gpsdata[i].reduction )
      ax3.errorbar(gpsyp,gpsupar,yerr = gpssigmapar,ecolor = 'blue',barsabove = 'True',fmt = "none")
      ax3.plot(gpsyp,gpsuperp,markers[i],color = 'green',mew = 1.,\
        label = '%s fault-perpendicular velocities'%gpsdata[i].reduction)
      ax3.errorbar(gpsyp,gpsuperp,yerr = gpssigmaperp,ecolor = 'green',fmt = "none")

      logger.debug('Number of GPS left within profile {0}'.format(len(gpsyp))) 

      if 3 == gps.dim:
          gpsuv,gpssigmav = np.delete(gps.uv,index), np.delete(gps.sigmav,index)

          # plot gps los
          ax2.plot(gpsyp,gpsulos,'+',color='black',mew=5.,label='%s GPS LOS'%gpsdata[i].reduction)
          
          ax3.plot(gpsyp,gpsuv,markers[i],color = 'red',mew = 1.,label = '%s vertical velocities'%gpsdata[i].reduction)
          ax3.errorbar(gpsyp,gpsuv,yerr = gpssigmav,ecolor = 'red',fmt = "none")          

      # set born profile equal to map
      logger.debug('Set ylim GPS profile to {0}-{1}'.format(gpsmin,gpsmax))
      ax3.set_ylim([gpsmin,gpsmax])

  colors = ['blue','red','orange','magenta']
  cst=0
  for i in xrange(len(insardata)):
      insar=insardata[i]
      losmin=insar.lmin
      losmax=insar.lmax

      logger.info('Load InSAR {0}'.format(insar.network)) 

      # perp and par composante ref to the profile 
      insar.ypp=(insar.x-profiles[k].x)*profiles[k].n[0]+(insar.y-profiles[k].y)*profiles[k].n[1]
      insar.xpp=(insar.x-profiles[k].x)*profiles[k].s[0]+(insar.y-profiles[k].y)*profiles[k].s[1]

      # select data within profile
      index=np.nonzero((insar.xpp>xpmax)|(insar.xpp<xpmin)|(insar.ypp>ypmax)|(insar.ypp<ypmin))
      insar.uu,insar.xx,insar.yy,insar.xxpp,insar.yypp=np.delete(insar.ulos,index),np.delete(insar.x,index),\
      np.delete(insar.y,index),np.delete(insar.xpp,index),np.delete(insar.ypp,index)

      logger.debug('Number of InSAR point left within profile {0}'.format(len(insar.uu))) 
     
      if len(insar.uu) > 50:

        nb = np.float(l/(len(insar.uu)/20.))
        logger.debug('Create bins every {0:.3f} km'.format(nb)) 

        bins = np.arange(-l/2-1,l/2+1,nb)
        inds = np.digitize(insar.yypp,bins)
        insar.distance = []
        insar.moy_los = []
        insar.std_los = []
        insar.xperp = []
        insar.yperp = []
        insar.uulos =  []
 
        for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            # remove NaN
            kk = np.flatnonzero(~np.isnan(insar.uu[uu]))
            _los = np.copy(insar.uu[uu][kk])
            _xperp = np.copy(insar.xxpp[uu][kk])
            _yperp = np.copy(insar.yypp[uu][kk])

            
            if len(kk)>10:
                insar.distance.append(bins[j] + (bins[j+1] - bins[j])/2.)

                indice = np.flatnonzero(np.logical_and(_los>np.percentile(\
                  _los,100-insar.perc),_los<np.percentile(_los,insar.perc)))

                insar.std_los.append(np.std(_los[indice]))
                insar.moy_los.append(np.median(_los[indice]))
                insar.xperp.append(_xperp[indice])
                insar.yperp.append(_yperp[indice])
                insar.uulos.append(_los[indice])
            else:
                logger.debug('{} points within the bin'.format(len(kk)))
                logger.debug('Less than 10 points within the bin. Nothing to be plot')

        del _los; del _xperp; del _yperp
        insar.distance = np.array(insar.distance)
        insar.std_los = np.array(insar.std_los)
        insar.moy_los = np.array(insar.moy_los)

        try:
          insar.xperp = np.concatenate(np.array(insar.xperp))
          insar.yperp = np.concatenate(np.array(insar.yperp))
          insar.uulos = np.concatenate(np.array(insar.uulos))
        except:
          insar.xperp = np.array(insar.xperp)
          insar.yperp = np.array(insar.yperp)
          insar.uulos = np.array(insar.uulos)

        # PLOT
        if typ is 'distscale':
          logger.info('Plot InSAR with distscale option')
          # colorscale fct of the parrallel distance to the profile
          norm = matplotlib.colors.Normalize(vmin=xpmin, vmax=xpmax)
          m1 = cm.ScalarMappable(norm=norm,cmap='cubehelix_r')
          m1.set_array(insar.xperp)
          facelos=m1.to_rgba(insar.xperp)
          ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,\
             label=insardata[i].reduction,color=facelos, rasterized=True)
        
        elif typ is 'std':
          logger.info('Plot InSAR with std option')
          # plot mean and standard deviation
          # ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,\
          #    label=insardata[i].reduction,color=colors[i])
          ax2.plot(insar.distance,insar.moy_los,color=insar.color,lw=3.,label=insardata[i].reduction)
          ax2.plot(insar.distance,insar.moy_los-insar.std_los,color=insar.color,lw=.5)
          ax2.plot(insar.distance,insar.moy_los+insar.std_los,color=insar.color,lw=.5)

        elif typ is 'stdscat':
          logger.info('Plot InSAR with stdscat option')
          # plot mean and standard deviation
          # ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,\
          #    label=insardata[i].reduction,color=colors[i])
          ax2.plot(insar.distance,insar.moy_los,color='black',lw=3.,label=insardata[i].reduction)
          ax2.plot(insar.distance,insar.moy_los-insar.std_los,color='black',lw=.5)
          ax2.plot(insar.distance,insar.moy_los+insar.std_los,color='black',lw=.5)
          ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,color=insar.color,rasterized=True)

        else:
          # plot scattering plot
          logger.info('No type profile give. Plot InSAR scatter point')
          ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,color=insar.color,rasterized=True)

        cst+=1.
      
        # set born profile equal to map
        logger.debug('Set ylim InSAR profile to {0}-{1}'.format(losmin,losmax))
        ax2.set_ylim([losmin,losmax])

        # for j in xrange(Mfault):
        # ax2.plot([fperp[j],fperp[j]],[losmax,losmin],color='red')
        # ax2.plot([fperp[j],fperp[j]],[losmax+cst,losmin-cst],color='red')
        # ax2.plot([fperp[j],fperp[j]],[6,-4],color='red')

      else:
          logger.critical('Number of InSAR points inferior to 50. Exit plot profile!') 
          
  if k is not len(profiles)-1:
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)

  if len(insardata) > 0:  
    if typ is 'distscale':
      fig2.colorbar(m1,shrink=0.5, aspect=5)
    else:
      ax2.legend(loc='best')

ax1.set_xlabel('Distance (km)')
ax1.set_ylabel('Elevation (km)')

ax2.set_xlabel('Distance (km)')
ax2.set_ylabel('LOS velocity (mm/yr)')

logger.debug('Save {0} output file'.format(outdir+profiles[k].name+'protopo.eps'))
fig1.savefig(outdir+profiles[k].name+'protopo.eps', format='EPS', dpi=150)

logger.debug('Save {0} output file'.format(outdir+profiles[k].name+'prolos.pdf'))
fig2.savefig(outdir+profiles[k].name+'prolos.pdf', format='PDF',dpi=150)

logger.debug('Save {0} output file'.format(outdir+profiles[k].name+'progps.eps'))
fig3.savefig(outdir+profiles[k].name+'progps.eps', format='EPS',)

logger.debug('Save {0} output file'.format(outdir+profiles[k].name+'promap.eps'))
fig.savefig(outdir+profiles[k].name+'promap.eps', format='EPS',)

plt.show()


