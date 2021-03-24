#!/usr/bin/env python3 

#from __future__ import print_function
import numpy as np
import scipy.optimize as opt
import scipy.linalg as lst

from osgeo import gdal
import pandas
import geopandas as gpd
import shapely.speedups

from matplotlib import pyplot as plt
#from matplotlib import mpl
import matplotlib
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from readgmt import *
from network2d import *
from model2d import *
from readgmt import *

from sys import argv,exit,stdin,stdout
import getopt
import os, math
from os import path
import logging

def hdi(trace, cred_mass=0.95):
    hdi_min, hdi_max = np.nanpercentile(trace,2.), np.nanpercentile(trace,98.)  
    return hdi_min, hdi_max

def usage():
  print('plotPro.py infile.py [-v] [-h]')
  print('-v Verbose mode. Show more information about the processing')
  print('-h Show this screen')

#load input file 
try:
    opts,args = getopt.getopt(sys.argv[1:], "h", ["help"])
except:
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
        format='%(lineno)s -- %(levelname)s -- %(message)s')
    logger = logging.getLogger('plotPro.log')
    logger.info('Initialise log file {0} in DEBUG mode'.format('plotPro.log'))

else:
    logging.basicConfig(level=logging.INFO,\
        format='%(lineno)s -- %(levelname)s -- %(message)s')
    logger = logging.getLogger('plotPro.log')
    logger.info('Initialise log file {0} in INFO mode. Use option -v for a DEBUG mode'.format('plotPro.log'))

if 1==len(sys.argv):
  usage()
  assert False, "no input file"
  logger.critical('No input file')
  sys.exit()

fname=sys.argv[1]
exec(open(fname).read())
if len(sys.argv)>1:
  try:
    fname=sys.argv[1]
    logger.info('Read input file {0}'.format(fname))
    try:
      sys.path.append(path.dirname(path.abspath(fname)))
      exec ("from "+path.basename(fname)+" import *")
    except:
      exec(open(fname).read())
      #execfile(path.abspath(fname))
  
  except Exception as e: 
    logger.critical('Problem in input file')
    logger.critical(e)
    print(network.__doc__)
    print(topo.__doc__)
    print(profile.__doc__)
    print(gmt.__doc__)
    print(fault2d.__doc__)
    sys.exit()

if 'xmin' in locals():
  xmin =xmin*1e3; xmax=xmax*1e3
  ymin = ymin*1e3; ymax = ymax*1e3
  logger.info('Found boundaries map plot {0:.2e} - {1:.2e} and {2:.2e} - {3:.2e} in locals'.format(xmin,xmax,ymin,ymax))
  xlim=[xmin,xmax];ylim=[ymin,ymax]
  extent = (xmin, xmax, ymin, ymax)
else:
  logger.info('Did not found  boundaries map in input file, set xmin,xmax,ymin,ymax')
  xlim=None;ylim=None
  extent = None

if not os.path.exists(outdir):
    logger.info('Creating output directory {0}'.format(outdir))
    os.makedirs(outdir)

# distance between fault and center of profile in fault azimuth 
# Fault model
try:
  Mfault=len(fmodel)
except:
  Mfault=0
  fmodel=[]
fperp=np.zeros(Mfault)

# Info basemap
if 'plot_basemap' not in locals():
    plot_basemap = False

# Load data
if 'topodata' not in globals():
    topodata = []
    logger.warning('No topodata list defined')
Mtopo = len(topodata)
for i in range(Mtopo):
  plot=topodata[i]
  logger.debug('Load data {0}'.format(plot.name))
  plot.load(xlim=xlim,ylim=ylim)
  if len(plot.z) < 1:
        logger.debug('Empty data file...')
if Mtopo == 0: 
  logger.warning('No topodata defined')
  Mtopo = 0

if 'gpsdata' not in globals():
    gpsdata = []
    logger.warning('No gpsdata list defined')
Mgps = len(gpsdata)
for i in range(Mgps):
      gps = gpsdata[i]
      logger.debug('Load data {0}'.format(gps.network))
      gps.loadgps()
      crs = gps.utm_proj       

if 'gmtfiles' not in globals():
    gmtfiles = []
    logger.warning('No gmtfiles list defined')

if 'seismifiles' not in globals():
    seismifiles = []
    logger.warning('No seismifiles list defined')
Mseismi=len(seismifiles)
for i in range(Mseismi):
    seismi = seismifiles[i]
    logger.debug('Load data {0}'.format(seismi.filename))
    seismi.load(xlim=xlim,ylim=ylim)

if 'shapefiles' not in globals():
    shapefiles = []
    logger.warning('No shapefiles list defined')

if 'insardata' not in globals():
    insardata = []
    logger.warning('No insardata list defined')
Minsar = len(insardata)
for i in range(Minsar):
    insar = insardata[i]
    logger.debug('Load data {0}'.format(insar.network))
    insar.loadinsar()
    if insar.theta == True:
      logger.warning('Convert LOS displacements to mean LOS angle assuming \
        horizontal displacements...')
      logger.warning('Use option theta=False in network class to avoid that')
      insar.losm = np.mean(insar.los)
      insar.uloscor = insar.ulos * \
        (np.sin(np.deg2rad(insar.losm))/np.sin(np.deg2rad(insar.los)))
    else:
      insar.uloscor = insar.ulos
    crs = insar.utm_proj

# except:
#   logger.warning('No insardata defined')
#   Minsar = 0

# MAP
# check if vertical for GPS
vertical_map = False
for i in range(Mgps):
  if gpsdata[i].dim == 3:
    vertical_map = True
if vertical_map:
  fig=plt.figure(0,figsize = (14,6))
  ax = fig.add_subplot(1,2,1)
  ax.axis('equal')
  ax12 = fig.add_subplot(1,2,2)
  ax.axis('equal')
  ax12.axis('equal')
else:
  fig=plt.figure(0,figsize = (9,8))
  ax = fig.add_subplot(1,1,1)
  ax.axis('equal')

logger.info('Plot Map ....') 
if 'xmin' in locals(): 
  ax.set_xlim(xmin,xmax)
  ax.set_ylim(ymin,ymax)
  if vertical_map:
    ax12.set_xlim(xmin,xmax)
    ax12.set_ylim(ymin,ymax)

if plot_basemap == True: 
    import contextily as ctx
    if extent == None:
       ax.axis(extent)
       if vertical_map:
         ax12.axis(extent)
    ctx.add_basemap(ax,crs="EPSG:{}".format(crs), source=ctx.providers.OpenTopoMap,alpha=0.5,zorder=0)
    if vertical_map:
      ctx.add_basemap(ax12,crs="EPSG:{}".format(crs), source=ctx.providers.OpenTopoMap,alpha=0.5,zorder=0)
else:
    print('plot_basemap variable is not defined or is not True. Skip backgroup topography plot')

for ii in range(len(gmtfiles)):
  name = gmtfiles[ii].name
  wdir = gmtfiles[ii].wdir
  filename = gmtfiles[ii].filename
  color = gmtfiles[ii].color
  width = gmtfiles[ii].width
  fx,fy = gmtfiles[ii].load(xlim=xlim,ylim=ylim)
  for i in range(len(fx)):
    ax.plot(fx[i],fy[i],color = color,lw = width,zorder=1)
  if vertical_map:
    ax12.plot(fx[i],fy[i],color = color,lw = width,zorder=1)

for ii in range(len(shapefiles)):
  name = shapefiles[ii].name
  fname = shapefiles[ii].filename
  wdir = shapefiles[ii].wdir
  color = shapefiles[ii].color
  edgecolor = shapefiles[ii].edgecolor
  linewidth = shapefiles[ii].linewidth
  crs = shapefiles[ii].crs
  shape = gpd.read_file(wdir + fname)
  if crs != None:
    shape = shape.to_crs("EPSG:{}".format(crs))
  shape.plot(ax=ax,facecolor='none', color=color,edgecolor=edgecolor,linewidth=linewidth,label=name,zorder=1)
  if vertical_map:
    shape.plot(ax=ax12,facecolor='none', color=color,edgecolor=edgecolor,linewidth=linewidth,label=name,zorder=1)

for ii in range(len(seismifiles)):
  name = seismifiles[ii].name
  x,y = seismifiles[ii].x, seismifiles[ii].y
  wdir = seismifiles[ii].wdir
  color = seismifiles[ii].color
  smin = np.nanmin(seismifiles[ii].mag)
  width = (seismifiles[ii].mag - smin)*seismifiles[ii].width*5
  ax.scatter(x,y,c=color,marker='o',s=width,linewidths=1, edgecolor='black',alpha=0.5,label=seismifiles[ii].name,zorder=2) 

try:
  from matplotlib.colors import LinearSegmentedColormap
  cm_locs = os.environ["Flower2d"] + '/contrib/colormaps/'
  cmap = LinearSegmentedColormap.from_list('roma', np.loadtxt(cm_locs+"roma.txt"))
  cmap = cmap.reversed()
except:
  cmap = cm.rainbow

for i in range(Minsar):
  insar=insardata[i]
  samp = insar.samp*4

  if (insar.lmin or insar.lmax) == None:
       vmin = np.nanpercentile(insar.ulos, 1)    
       vmax = np.nanpercentile(insar.ulos, 99)
  else:
       vmin, vmax = insar.lmin, insar.lmax

  logger.info('Plot data in map view {0} between {1} and {2}'.format(insar.network, vmin, vmax))
  logger.info('Subsample data every {0} point (samp option)'.format(insar.samp))
  norm = matplotlib.colors.Normalize(vmin=insar.lmin, vmax=insar.lmax)
  m = cm.ScalarMappable(norm = norm, cmap = cmap)
  m.set_array(insar.ulos[::samp])
  masked_array = np.ma.array(insar.ulos[::samp], mask=np.isnan(insar.ulos[::samp]))
  facelos = m.to_rgba(masked_array)
  ax.scatter(insar.x[::samp],insar.y[::samp], s=1, marker = 'o',color = facelos, rasterized=True, label = 'LOS Velocity {}'.format(insar.reduction),zorder=3)

gpscolor = ['black','coral','red','darkorange']
for i in range(Mgps):
  gps=gpsdata[i]
  logger.info('Plot GPS data {0}'.format(gps.network))
  ax.quiver(gps.x,gps.y,gps.ux,gps.uy,scale = 150, width = 0.003, color = gpscolor[i%4],zorder=4)
  if vertical_map:
    norm = matplotlib.colors.Normalize(vmin=np.nanpercentile(gps.uv,5), vmax=np.nanpercentile(gps.uv,95))
    mv = cm.ScalarMappable(norm = norm, cmap = cmap)
    mv.set_array(gps.uv)
    facev = mv.to_rgba(gps.uv)
    ax12.scatter(gps.x,gps.y,c=facev,marker='o',s=40,linewidths=1, edgecolor='black',alpha=0.8,label='Vertical velocities network {}'.format(gps.reduction),zorder=2)    

  if gps.plotName == True:
      for kk in range(len(gps.name)):
            ax.text(gps.x[kk], gps.y[kk], gps.name[kk], color ='black')

# add colorbar los
if 'facelos' in locals():
  divider = make_axes_locatable(ax)
  c = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(m, cax=c) 

# add colorbar vertical
if 'facev' in locals():
  divider = make_axes_locatable(ax)
  c = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(mv, cax=c) 
  divider = make_axes_locatable(ax12)
  c = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(mv, cax=c) 

# plot legend
ax.legend(loc = 'upper right',fontsize='x-small')
if vertical_map:
  ax12.legend(loc = 'upper right',fontsize='x-small')
plt.tight_layout()  

# fig pro topo
if len(profiles) > 1:
  fig1=plt.figure(1,figsize=(10,8))
else:
  fig1=plt.figure(1,figsize=(10,3))
fig1.subplots_adjust(hspace=0.0001)

# fig pro insar
if Minsar>0:
  if len(profiles) > 1:
    fig2=plt.figure(4,figsize=(10,8))
  else:
    fig2=plt.figure(4,figsize=(10,3))
    fig2.subplots_adjust(hspace=0.0001)

# fig pro gps
if Mgps>0:
  if len(profiles) > 1:
    fig3=plt.figure(5,figsize=(10,8))
  else:
    fig3=plt.figure(5,figsize=(10,3))
    fig3.subplots_adjust(hspace=0.0001)

if len(seismifiles)>0:
  if len(profiles) > 1:
    fig4=plt.figure(6,figsize=(10,8))
  else:
    fig4=plt.figure(6,figsize=(10,3))
    fig4.subplots_adjust(hspace=0.0001)

logger.info('Plot Profiles ....')

# Plot profile
for k in range(len(profiles)): 

  l=profiles[k].l
  w=profiles[k].w
  x0=profiles[k].x
  y0=profiles[k].y
  strike = profiles[k].strike
  name=profiles[k].name
  typ=profiles[k].typ
  flat=profiles[k].flat
  nb = profiles[k].lbins
  loc_ramp = profiles[k].loc_ramp

  logger.info('Plot profile {0}. length: {1}, width :{2}, strike: {3}'.format(name, l, w, strike)) 

  # lim profile
  ypmax,ypmin=l/2,-l/2
  xpmax,xpmin=w/2,-w/2

  # profile azimuth
  profiles[k].str=(profiles[k].strike*math.pi)/180
  profiles[k].s=[math.sin(profiles[k].str),math.cos(profiles[k].str),0]
  profiles[k].n=[math.cos(profiles[k].str),-math.sin(profiles[k].str),0]

  for j in range(Mfault):
    fperp[j]=(fmodel[j].x-profiles[k].x)*profiles[k].n[0]+(fmodel[j].y-profiles[k].y)*profiles[k].n[1]

  ax1=fig1.add_subplot(len(profiles),1,k+1)
  ax1.set_xlim([-l/2,l/2])

  for i in range(Mtopo):
        plot=topodata[i]

        # perp and par composante ref to the profile 
        plot.ypp=(plot.x-profiles[k].x)*profiles[k].n[0]+(plot.y-profiles[k].y)*profiles[k].n[1]
        plot.xpp=(plot.x-profiles[k].x)*profiles[k].s[0]+(plot.y-profiles[k].y)*profiles[k].s[1]

        index=np.nonzero((plot.xpp>xpmax)|(plot.xpp<xpmin)|(plot.ypp>ypmax)|(plot.ypp<ypmin))
        plotxpp,plotypp,plotz=np.delete(plot.xpp,index),np.delete(plot.ypp,index),np.delete(plot.z,index)
        if nb == None:
          nb = np.float(l/(len(plotz)/100.))
          logger.info('Create bins every {0:.3f} km'.format(nb)) 
        else:
          logger.info('Set nbins to {}, defined in profile class'.format(nb))

        bins = np.arange(-l/2,l/2, nb/2.)
        inds = np.digitize(plotypp,bins)
        distance = []
        moy_topo = []
        std_topo = []
        for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            if len(uu)>0:
                distance.append(bins[j] + (bins[j+1] - bins[j])/2.)
                std_topo.append(np.nanstd(plotz[uu]))
                moy_topo.append(np.nanmedian(plotz[uu]))
        
        distance = np.array(distance)
        std_topo = np.array(std_topo)
        moy_topo = np.array(moy_topo)

        ax1.plot(distance,moy_topo,label=plot.name,color=plot.color,lw=plot.width)
        if plot.plotminmax == True:
          logger.debug('plotminmax set to True')
          ax1.plot(distance,moy_topo-std_topo,color=plot.color,lw=plot.width/2)
          ax1.plot(distance,moy_topo+std_topo,color=plot.color,lw=plot.width/2)
        
        if plot.axis == 'equal':
          ax1.axis('equal')
        else:
          if (plot.topomin != None) and (plot.topomax != None) :
              logger.info('Set ylim to {} and {}'.format(plot.topomin,plot.topomax))
              ax1.set_ylim([plot.topomin,plot.topomax])
              for kk in range(Mfault):    
                ax1.plot([fperp[kk],fperp[kk]],[plot.topomin,plot.topomax],color='red')
                ax1.text(fperp[kk],0.5,fmodel[kk].name,color='red')
          else:
              topomin,topomax= ax1.get_ylim()
              for kk in range(Mfault):    
                ax1.plot([fperp[kk],fperp[kk]],[topomin,topomax],color='red')
                ax1.text(fperp[kk],0.5,fmodel[kk].name,color='red')

  # LOS profile/map
  if Minsar>0:
    ax2=fig2.add_subplot(len(profiles),1,k+1)
    ax2.set_xlim([-l/2,l/2])

  # LOS profile/map
  if Mgps>0:
    ax3=fig3.add_subplot(len(profiles),1,k+1)
    ax3.set_xlim([-l/2,l/2])

  # depth profile
  if len(seismifiles)>0:
     ax4=fig4.add_subplot(len(profiles),1,k+1)
     ax4.axis('equal')
     ax4.set_xlim([-l/2,l/2])
      
  for ii in range(len(seismifiles)):
    name = seismifiles[ii].name
    x,y = seismifiles[ii].x, seismifiles[ii].y
    wdir = seismifiles[ii].wdir
    color = seismifiles[ii].color

    # project in profile
    seismi.ypp=(x-profiles[k].x)*profiles[k].n[0]+(y-profiles[k].y)*profiles[k].n[1]
    seismi.xpp=(x-profiles[k].x)*profiles[k].s[0]+(y-profiles[k].y)*profiles[k].s[1]
    # select data within profile
    index=np.nonzero((seismi.xpp>xpmax)|(seismi.xpp<xpmin)|(seismi.ypp>ypmax)|(seismi.ypp<ypmin))
    seismi.xp,seismi.yp,depth,size=np.delete(seismi.xpp,index),np.delete(seismi.ypp,index),np.delete(seismifiles[ii].depth,index),np.delete(seismifiles[ii].mag,index)
    smin = np.nanmin(size)
    size = (size - smin)* np.float(seismifiles[ii].width)*5
    # plot
    ax4.scatter(seismi.yp,-depth,s=size,c=color,marker='o',linewidths=1, edgecolor='black',alpha=0.5,label=seismifiles[ii].name,zorder=2) 

  # plot profiles
  xp,yp = np.zeros((7)),np.zeros((7))
  xp[:] = x0-w/2*profiles[k].s[0]-l/2*profiles[k].n[0],x0+w/2*\
  profiles[k].s[0]-l/2*profiles[k].n[0],x0+w/2*profiles[k].s[0]+l/2*profiles[k].n[0],x0-w/2*profiles[k].s[0]+l/2*profiles[k].n[0],x0-w/2*profiles[k].s[0]-l/2*profiles[k].n[0],x0-l/2*profiles[k].n[0],x0+l/2*profiles[k].n[0]
  yp[:] = y0-w/2*profiles[k].s[1]-l/2*profiles[k].n[1],y0+w/2*\
  profiles[k].s[1]-l/2*profiles[k].n[1],y0+w/2*profiles[k].s[1]+l/2*profiles[k].n[1],y0-w/2*profiles[k].s[1]+l/2*profiles[k].n[1],y0-w/2*profiles[k].s[1]-l/2*profiles[k].n[1],y0-l/2*profiles[k].n[1],y0+l/2*profiles[k].n[1]

  # plot in map view  
  ax.plot(xp[:],yp[:],color = 'black',lw = 1.)
  if vertical_map:
    ax12.plot(xp[:],yp[:],color = 'black',lw = 1.)

  # GPS plot
  markers = ['+','d','x','v']
  for i in range(Mgps):
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

      # compute fault parallel and perpendicular for each profiles
      gpsupar = gpsux*profiles[k].s[0]+gpsuy*profiles[k].s[1]
      gpsuperp = gpsux*profiles[k].n[0]+gpsuy*profiles[k].n[1]
      gpssigmaperp=((gpssigmax*np.cos(profiles[k].str))**2 + (gpssigmay*np.sin(profiles[k].str))**2)**0.5
      gpssigmapar=((gpssigmax*np.sin(profiles[k].str))**2 + (gpssigmay*np.cos(profiles[k].str))**2)**0.5

      ax3.plot(gpsyp,gpsupar,markers[i],color = 'blue',mew = 1.5,label =\
       '%s fault-parallel velocities'%gpsdata[i].reduction )
      ax3.errorbar(gpsyp,gpsupar,yerr = gpssigmapar,ecolor = 'blue',barsabove = 'True',fmt = "none",alpha=.5)
      ax3.plot(gpsyp,gpsuperp,markers[i],color = 'green',mew = 1.5,\
        label = '%s fault-perpendicular velocities'%gpsdata[i].reduction)
      ax3.errorbar(gpsyp,gpsuperp,yerr = gpssigmaperp,ecolor = 'green',fmt = "none",alpha=.5)

      logger.debug('Number of GPS left within profile {0}'.format(len(gpsyp))) 

      if 3 == gps.dim:
          gpsuv,gpssigmav,gpsulos,gpssigmalos = np.delete(gps.uv,index), np.delete(gps.sigmav,index),np.delete(gps.ulos,index),np.delete(gps.sigmalos,index)

          ax3.plot(gpsyp,gpsuv,markers[i],color = 'red',mew = 1.5,label = '%s vertical velocities'%gpsdata[i].reduction)
          ax3.errorbar(gpsyp,gpsuv,yerr = gpssigmav,ecolor = 'red',fmt = "none",alpha=.5)          
          
          if gps.proj != None:
            # plot gps los
            ax2.plot(gpsyp,gpsulos,'+',color='red',mew=5.,label='%s GPS LOS'%gpsdata[i].reduction)
            ax2.errorbar(gpsyp,gpsulos,yerr = gpssigmalos,ecolor = 'red',fmt = "none")          

      for j in range(Mfault):
          ax3.plot([fperp[j],fperp[j]],[gpsmax,gpsmin],color='red')
     
      # plot vertical lines
      ax3.hlines(np.linspace(gpsmin,gpsmax,6),xmin=-l/2,xmax=l/2,ls='--',color='black',lw=.5)

      # set born profile equal to map
      logger.debug('Set ylim GPS profile to {0}-{1}'.format(gpsmin,gpsmax))
      ax3.set_ylim([gpsmin,gpsmax])

  if Mgps>0:
    ax3.legend(loc = 'best',fontsize='x-small')
  if len(seismifiles)>0:
    ax4.legend(loc = 'best',fontsize='x-small')
          
  cst=0
  for i in range(Minsar):
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
       
      # Initialise for plot in case no data for this profile
      insar.distance = []
      insar.moy_los = []
      insar.std_los = []
      insar.xperp = []
      insar.yperp = []
      insar.uulos =  []    

      if len(insar.uu) > 50:

        if nb == None:
          nb = np.float(l/(len(insar.uu)/100.))
          logger.info('Create bins every {0:.3f} km'.format(nb)) 
        else:
          logger.info('Set nbins to {} defined in profile class'.format(nb)) 

        bins = np.arange(-l/2-1,l/2+1,nb)
        inds = np.digitize(insar.yypp,bins)
 
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

                insar.std_los.append(np.nanstd(_los[indice]))
                insar.moy_los.append(np.nanmedian(_los[indice]))
                insar.xperp.append(_xperp[indice])
                insar.yperp.append(_yperp[indice])
                insar.uulos.append(_los[indice])
            else:
                logger.debug('{} points within the bin'.format(len(kk)))
                logger.debug('Less than 10 points within the bin. Nothing to be plot')

        del _los; del _xperp; del _yperp

        try:
          insar.xperp = np.concatenate(np.array(insar.xperp))
          insar.yperp = np.concatenate(np.array(insar.yperp))
          insar.uulos = np.concatenate(np.array(insar.uulos))
          insar.distance = np.concatenate(np.array(insar.distance))
          insar.std_los = np.concatenate(np.array(insar.std_los))
          insar.moy_los = np.concatenate(np.array(insar.moy_los))
        except:
          insar.xperp = np.array(insar.xperp)
          insar.yperp = np.array(insar.yperp)
          insar.uulos = np.array(insar.uulos)
          insar.distance = np.array(insar.distance)
          insar.std_los = np.array(insar.std_los)
          insar.moy_los = np.array(insar.moy_los)

      else:
          logger.critical('Number of InSAR points inferior to 50 for track {}. Exit plot profile!'.format(insar.reduction)) 

  # FLATEN
  if (flat != None):

    # remove ramp between two profiles
    if len(insardata)==2:

      logger.info('Flat is not None and 2 InSAR network defined: flattening based on the differences in the overlapping areas')
      insar1,insar2=insardata[0],insardata[1]
      kk2 = np.in1d(insar2.distance, insar1.distance)
      kk2 = np.flatnonzero(kk2==True)
      kk1 = np.in1d(insar1.distance, insar2.distance)
      kk1 = np.flatnonzero(kk1==True)

      temp_los = (insar1.moy_los[kk1] - insar2.moy_los[kk2])[::]
      temp_yp = insar1.distance[kk1][::]
      temp_std = np.sqrt((insar1.std_los[kk1]**2 + insar2.std_los[kk2]**2))[::]

      # # cut longueurs tracks
      kmax1,kmax2=np.max(insar1.distance), np.max(insar2.distance)
      kmin1,kmin2=np.min(insar1.distance), np.min(insar2.distance)
      kmax,kmin = np.min([kmax1,kmax2]), np.max([kmin1,kmin2])

    # remove ramp along profile on one LOS
    else:

      logger.info('Flat is not None and 1 InSAR network defined: flattening along the profile')
      insar2 = insardata[0]
      kmax,kmin = np.max(insar2.distance), np.min(insar2.distance)

      if loc_ramp=="positive":
        logger.info('Estimate ramp in the postive distances of the profile')
        kk2 = np.nonzero((insar2.distance>0))
      elif loc_ramp=="negative":
        logger.info('Estimate ramp in the negative distances of the profile')
        kk2 = np.nonzero((insar2.distance<0))
      else:
        logger.info('Estimate ramp within the whole profile')
        kk2 = np.arange(len(insar2.distance))

      temp_los = insar2.moy_los[kk2]
      temp_yp = insar2.distance[kk2]
      temp_std = insar2.std_los[kk2]

    # temp_std = np.ones(len(temp_yp))
    # Cd = np.diag(temp_std**2,k=0)
    # remove residuals NaN
    kk = np.flatnonzero(~np.isnan(temp_los))
    temp_los,temp_yp,temp_std = temp_los[kk],temp_yp[kk],temp_std[kk]
   
    if flat == 'quad': 
        G = np.zeros((len(temp_los),3))
        G[:,0] = temp_yp**2
        G[:,1] = temp_yp
        G[:,2] = 1
    elif flat == 'cub':
        G = np.zeros((len(temp_los),4))
        G[:,0] = temp_yp**3
        G[:,1] = temp_yp**2
        G[:,2] = temp_yp
        G[:,3] = 1
    else:
        G = np.zeros((len(temp_los),2))
        G[:,0] = temp_yp
        G[:,1] = 1

    try:
      x0 = lst.lstsq(G,temp_los)[0]
    except Exception as e:
      logger.warning(e)
      x0 = np.zeros(np.shape(G)[1])

    # print x0
    _func = lambda x: np.sum(((np.dot(G,x)-temp_los)/temp_std)**2)
    _fprime = lambda x: 2*np.dot(G.T/temp_std, (np.dot(G,x)-temp_los[::])/temp_std)
    pars = opt.fmin_slsqp(_func,x0,fprime=_fprime,iter=2000,full_output=True,iprint=0)[0]

    # if flat == 'quad': 
    #   a = pars[0]; b = pars[1]; c = pars[2]
    #   blos = a*insar2.distance[kk2]**2 + b*insar2.distance[kk2] + c
    # elif flat == 'cub':
    #   a = pars[0]; b = pars[1]; c = pars[2]; d =pars[3]
    #   blos = a*insar2.distance[kk2]**3 + b*insar2.distance[kk2]**2 + c*insar2.distance[kk2] + d
    # else:
    #   a = pars[0]; b = pars[1]
    #   blos = a*insar2.distance[kk2] + b
    
    if flat == 'quad':
        a = pars[0]; b = pars[1]; c = pars[2]
        logger.info('Remove ramp: {0} yperp**2  + {1} yperp  + {2}'.format(a,b,c))

        blos = a*insar2.distance**2 + b*insar2.distance + c
        insar2.moy_los = insar2.moy_los + blos
        diff = temp_los - blos[kk]

        blos = a*insar2.yperp**2 + b*insar2.yperp + c
        insar2.uulos = insar2.uulos + blos

        blos = a*insar2.yypp**2 + b*insar2.yypp + c
        insar2.uu = insar2.uu + blos

        blos = a*insar2.ypp**2 + b*insar2.ypp + c
        insar2.ulos = insar2.ulos + blos
        insar2.uloscor = insar2.uloscor + blos

        x = np.arange(kmin,kmax,1)
        ysp = a*x**2 + b*x + c

    elif flat == 'cub':
        a = pars[0]; b = pars[1]; c = pars[2]; d =pars[3]
        logger.info('Remove ramp: {0} yperp**3 + {1} yperp**2  + {2} yperp  + {3}'.format(a,b,c,d))
    
        blos = a*insar2.distance**3 + b*insar2.distance**2 + c*insar2.distance + d
        insar2.moy_los = insar2.moy_los + blos
        diff = temp_los - blos[kk]

        blos = a*insar2.yperp**3 + b*insar2.yperp**2 + c*insar2.yperp + d
        insar2.uulos = insar2.uulos + blos

        blos = a*insar2.yypp**3 + b*insar2.yypp**2 + c*insar2.yypp + d
        insar2.uu = insar2.uu + blos

        blos = a*insar2.ypp**3 + b*insar2.ypp**2 + c*insar2.ypp + d
        insar2.ulos = insar2.ulos + blos
        insar2.uloscor = insar2.uloscor + blos

        x = np.arange(kmin,kmax,1)
        ysp = a*x**3 + b*x**2 + c*x +d

    else:
        a = pars[0]; b = pars[1]
        logger.info('Remove ramp: {0} yperp  + {1}'.format(a,b))
    
        blos = a*insar2.distance + b
        insar2.moy_los = insar2.moy_los + blos
        diff = temp_los - blos[kk]
        
        blos = a*insar2.yperp + b
        insar2.uulos = insar2.uulos + blos

        blos = a*insar2.yypp + b
        insar2.uu = insar2.uu + blos

        blos = a*insar2.ypp + b
        insar2.ulos = insar2.ulos + blos
        insar2.uloscor = insar2.uloscor + blos

        x = np.arange(kmin,kmax,1)
        ysp = a*x + b

    if len(insardata)==2:

        # diff = temp_los - blos[kk]
        # coor1 = np.array([insardata[0].xx, insardata[0].yy])
        # coor2 = np.array([insardata[1].xx, insardata[1].yy]) 
        # indices = coor1 == coor2 

        # Plot histogram
        fig4=plt.figure(5,figsize=(9,6))
        ax4 = fig4.add_subplot(1,1,1)
        ax4.hist(diff,bins=40,density=True,histtype='stepfilled', \
          color='grey',alpha=0.4,label='{}-{}'.format(insardata[0].reduction,insardata[1].reduction))
        hdi_min, hdi_max = hdi(diff)
        opts = {'c':'green', 'linestyle':'--'}
        ax4.axvline(x=hdi_min, **opts)
        ax4.axvline(x=hdi_max, **opts)
        ax4.set_xlabel("Mean: {:0.3f}\n95% HDI: {:0.3f} - {:0.3f}".format(\
          diff.mean(), hdi_min, hdi_max))
        print("Compute Difference")
        print("Mean: {:0.3f}95% HDI: {:0.3f} - {:0.3f}".format(diff.mean(), hdi_min, hdi_max))
        print("Std: {:0.3f}".format(diff.std()))
        ax4.legend(loc='best')
        ax4.set_xlim(math.floor(np.nanmin(diff)),math.ceil(np.nanmax(diff)))
        logger.debug('Save {0} output file'.format(outdir+profiles[0].name+'_'+flat+'_histo.eps'))
        fig4.savefig(outdir+profiles[0].name+'_'+flat+'_histo.eps', format='EPS',dpi=150)
    
    # plot ramp
    ax2.plot(x,ysp,color='red',lw=1.,label='Estimated ramp')

  for i in range(Minsar):
        insar=insardata[i]
        losmin=insar.lmin
        losmax=insar.lmax

        if (flat != None) and len(insardata)==2:
            logger.info('Plot InSAR with std option')
            # plot mean and standard deviation
            ax2.plot(insar.distance,insar.moy_los,color=insar.color,lw=2.,label=insardata[i].reduction)
            ax2.plot(insar.distance,insar.moy_los-insar.std_los,color=insar.color,lw=.5)
            ax2.plot(insar.distance,insar.moy_los+insar.std_los,color=insar.color,lw=.5)
            ax2.scatter(insar.yperp,insar.uulos,s = .01, marker='o',alpha=0.1,color=insar.color,rasterized=True)

        else:
          if len(insar.distance) >0:
            # PLOT
            if typ == 'distscale':
              logger.info('Plot InSAR with distscale option')
              # colorscale fct of the parrallel distance to the profile
              norm = matplotlib.colors.Normalize(vmin=xpmin, vmax=xpmax)
              m1 = cm.ScalarMappable(norm=norm,cmap='cubehelix_r')
              m1.set_array(insar.xperp)
              facelos=m1.to_rgba(insar.xperp)
              ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,\
                 label=insardata[i].reduction,color=facelos, rasterized=True)
            
            elif typ == 'std':
              logger.info('Plot InSAR with std option')
              # plot mean and standard deviation
              ax2.plot(insar.distance,insar.moy_los,color=insar.color,lw=2.,label=insardata[i].reduction)
              ax2.plot(insar.distance,insar.moy_los-insar.std_los,color=insar.color,lw=.5)
              ax2.plot(insar.distance,insar.moy_los+insar.std_los,color=insar.color,lw=.5)

            elif typ == 'stdscat':
              logger.info('Plot InSAR with stdscat option')
              # plot mean and standard deviation
              ax2.plot(insar.distance,insar.moy_los,color=insar.color,lw=2.,label=insardata[i].reduction)
              ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.1,color=insar.color,rasterized=True)
              ax2.plot(insar.distance,insar.moy_los-insar.std_los,color='black',lw=.5)
              ax2.plot(insar.distance,insar.moy_los+insar.std_los,color='black',lw=.5)

            else:
              # plot scattering plot
              logger.info('No type profile give. Plot InSAR scatter point')
              ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,color=insar.color,rasterized=True)

            cst+=1.
          
        # set born profile equal to map
        if (losmin != None) and (losmax != None):
          logger.debug('Set ylim InSAR profile to {0}-{1}'.format(losmin,losmax))
          ax2.set_ylim([losmin,losmax])

  if Minsar>0:
    for j in range(Mfault):
      ax2.plot([fperp[j],fperp[j]],[losmax,losmin],color='red')
    if k != len(profiles)-1:
      plt.setp(ax2.get_xticklabels(), visible=False)
      plt.setp(ax1.get_xticklabels(), visible=False)
    if typ == 'distscale':
      fig2.colorbar(m1,shrink=0.5, aspect=5)
    else:
      ax2.legend(loc='best')

if (flat != None) and len(insardata)==2:
  logger.info('Plot fatten Maps...')
  # MAP
  fig=plt.figure(6,figsize = (9,7))
  ax = fig.add_subplot(1,1,1)
  ax.axis('equal')
  for i in range(len(insardata)):
    insar=insardata[i]
    samp = insar.samp*4

    logger.info('Plot data in map view {0} between {1} and {2}'.format(insar.network, vmin, vmax))
    logger.info('Subsample data every {0} point (samp option)'.format(insar.samp))
    norm = matplotlib.colors.Normalize(vmin=insar.lmin, vmax=insar.lmax)
    m = cm.ScalarMappable(norm = norm, cmap = 'rainbow')
    m.set_array(insar.ulos[::samp])
    masked_array = np.ma.array(insar.ulos[::samp], mask=np.isnan(insar.ulos[::samp]))
    facelos = m.to_rgba(masked_array)
    ax.scatter(insar.x[::samp],insar.y[::samp],s = 2,marker = 'o',color = facelos,\
      label = 'LOS Velocity %s'%(insar.reduction),rasterized=True)

    # save flatten map
    if i==1:
      np.savetxt('{}_flat'.format(insardata[i].network), np.vstack([insar.x,insar.y,insar.ulos]).T, fmt='%.6f')

    # plot faults
    for kk in range(Mfault):
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

    if 'xmin' in locals(): 
      ax.set_xlim(xmin,xmax)
      ax.set_ylim(ymin,ymax)

    for ii in range(len(gmtfiles)):
      name = gmtfiles[ii].name
      wdir = gmtfiles[ii].wdir
      filename = gmtfiles[ii].filename
      color = gmtfiles[ii].color
      width = gmtfiles[ii].width
      fx,fy = gmtfiles[ii].load()
      for i in range(len(fx)):
        ax.plot(fx[i],fy[i],color = color,lw = width)

  # add colorbar los
  if len(insardata) > 0:
    fig.colorbar(m,shrink = 0.5, aspect = 5)

  # plot profile
  ax.plot(xp[:],yp[:],color = 'black',lw = 1.)
  ax.set_title('Flatten LOS')

ax1.set_xlabel('Distance (km)')
ax1.set_ylabel('Elevation (km)')

if len(seismifiles)>0:
  ax4.set_xlabel('Distance (km)')
  ax4.set_ylabel('Elevation (m)')

plt.show()

if Minsar>0:
  ax2.set_xlabel('Distance (km)')
  ax2.set_ylabel('LOS velocity (mm)')
  logger.debug('Save {0} output file'.format(outdir+profiles[k].name+'pro-los.pdf'))
  fig2.savefig(outdir+profiles[k].name+'pro-los.pdf', format='PDF',dpi=150)

logger.debug('Save {0} output file'.format(outdir+profiles[k].name+'protopo.eps'))
fig1.savefig(outdir+profiles[k].name+'pro-topo.pdf', format='PDF', dpi=150)

if Mgps>0:
  logger.debug('Save {0} output file'.format(outdir+profiles[k].name+'progps.eps'))
  fig3.savefig(outdir+profiles[k].name+'pro-gps.pdf', format='PDF',dpi=75)

if len(seismifiles)>0 : 
  logger.debug('Save {0} output file'.format(outdir+profiles[k].name+'pro-depth.eps'))
  fig4.savefig(outdir+profiles[k].name+'pro-depth.pdf', format='PDF', dpi=150)

logger.debug('Save {0} output file'.format(outdir+profiles[k].name+'promap.eps'))
fig.savefig(outdir+profiles[k].name+'pro-map.pdf', format='PDF', dpi=150)



