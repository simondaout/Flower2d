#!/usr/bin/env python2.7 

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

#load input file 
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

if 1==len(sys.argv):
  usage()
  assert False, "no input file"

if 2==len(sys.argv):
  fname=sys.argv[1]
  print
  print 'input file :', fname
  sys.path.append(path.dirname(path.abspath(fname)))
  exec ("from "+path.basename(fname)+" import *")
  #exec ("from " + fname+ " import *")
else:
  assert False, "too many arguments"

if not os.path.exists(outdir):
    os.makedirs(outdir)

# distance between fault and center of profile in fault azimuth 
# Fault model
Mfault=len(fmodel)
fperp=np.zeros(Mfault)

# Load data
for i in xrange(len(topodata)):
    plot=topodata[i]
    plot.load()

for i in xrange(len(gpsdata)):
    gps = gpsdata[i]
    gps.loadgps()

for i in xrange(len(insardata)):
    insar = insardata[i]
    insar.loadinsar()
    if insar.theta is True:
		  sys.exit()
		  insar.losm = np.mean(insar.los)
		  insar.ulos = insar.ulos * \
        (np.sin(np.deg2rad(insar.losm))/np.sin(np.deg2rad(insar.los)))


# MAP
fig=plt.figure(0,figsize = (9,8))
ax = fig.add_subplot(1,1,1)
ax.axis('equal')

for i in xrange(len(insardata)):
  insar=insardata[i]
  samp = 1

  norm = matplotlib.colors.Normalize(vmin = losmin, vmax = losmax)
  m = cm.ScalarMappable(norm = norm, cmap = 'rainbow')
  m.set_array(insar.ulos[::samp])
  facelos = m.to_rgba(insar.ulos[::samp])
  #facelos = m.to_rgba(insarulos)
  ax.scatter(insar.x[::samp],insar.y[::samp],s = 2,marker = 'o',color = facelos,label = 'LOS Velocity %s'%(insar.reduction))
  #ax.scatter(insar.xx,insar.yy,s = 2,marker = 'o',color = facelos,alpha=0.5,label = 'LOS Velocity %s'%(insar.reduction))

  # save flatten map
  # np.savetxt('{}_flat'.format(insardata[i].network), np.vstack([insar.x,insar.y,insar.ulos]).T, fmt='%.6f')

for i in xrange(len(gpsdata)):
  gps=gpsdata[i]
  ax.quiver(gps.x,gps.y,gps.ux,gps.uy,scale = 100,width = 0.005,color = 'red')


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
    
# Plot profile
for k in xrange(len(profiles)): 

  l=profiles[k].l
  w=profiles[k].w
  x0=profiles[k].x
  y0=profiles[k].y
  name=profiles[k].name
  typ=profiles[k].typ

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
        topomax, topomin = np.max(-plotz), np.min(-plotz) 
        
        bins = np.arange(-l/2,l/2,1)
        inds = np.digitize(plotypp,bins)
        distance = []
        moy_topo = []
        std_topo = []
        for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            if len(uu)>0:
                distance.append(bins[j] + (bins[j+1] - bins[j])/2.)
                std_topo.append(np.std(plotz[uu]))
                moy_topo.append(np.mean(plotz[uu]))
        
        distance = np.array(distance)
        std_topo = np.array(std_topo)
        moy_topo = np.array(moy_topo)

        ax1.scatter(plotypp,-plotz,s=plot.width, marker='o',label=plot.name,color=plot.color)
        ax1.plot(distance,-moy_topo,label=plot.name,color='black',lw=1)
        #ax1.plot(distance,-moy_topo-std_topo,color='black',lw=1)
        #ax1.plot(distance,-moy_topo+std_topo,color='black',lw=1)
        
        # for kk in xrange(Mfault):    
            # ax1.plot([fperp[kk],fperp[kk]],[8,-8],color='red')
            # ax1.text(fperp[kk],0.5,fmodel[kk].name,color='red')
        
        #plt.ylim([-30,5])
        #plt.ylim([topomin,topomax])
        # plt.ylim([0.5,6.5])
        # plt.title('Profile %s'%(profiles[k].name))
  
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
      ax3.errorbar(gpsyp,gpsupar,yerr = gpssigmapar,ecolor = 'blue',barsabove = 'True',fmt = None)
      ax3.plot(gpsyp,gpsuperp,markers[i],color = 'green',mew = 1.,\
        label = '%s fault-perpendicular velocities'%gpsdata[i].reduction)
      ax3.errorbar(gpsyp,gpsuperp,yerr = gpssigmaperp,ecolor = 'green',fmt = None)

      if 3 == gps.dim:
          gpsuv,gpssigmav = np.delete(gps.uv,index), np.delete(gps.sigmav,index)

          # plot gps los
          ax2.plot(gpsyp,gpsulos,'+',color='black',mew=5.,label='%s GPS LOS'%gpsdata[i].reduction)
          
          ax3.plot(gpsyp,gpsuv,markers[i],color = 'red',mew = 1.,label = '%s vertical velocities'%gpsdata[i].reduction)
          ax3.errorbar(gpsyp,gpsuv,yerr = gpssigmav,ecolor = 'red',fmt = None)          

      # set born profile equal to map
      if 'gpsmin' in locals():
          ax3.set_ylim([gpsmin,gpsmax])

  colors = ['blue','red','orange','magenta']
  cst=0
  for i in xrange(len(insardata)):
      insar=insardata[i]
      # perp and par composante ref to the profile 
      insar.ypp=(insar.x-profiles[k].x)*profiles[k].n[0]+(insar.y-profiles[k].y)*profiles[k].n[1]
      insar.xpp=(insar.x-profiles[k].x)*profiles[k].s[0]+(insar.y-profiles[k].y)*profiles[k].s[1]

      # select data within profile
      index=np.nonzero((insar.xpp>xpmax)|(insar.xpp<xpmin)|(insar.ypp>ypmax)|(insar.ypp<ypmin))
      insar.uu,insar.xx,insar.yy,insar.xxpp,insar.yypp=np.delete(insar.ulos,index),np.delete(insar.x,index),\
      np.delete(insar.y,index),np.delete(insar.xpp,index),np.delete(insar.ypp,index)
      if len(insar.uu) > 1000:
        #binsarlos=np.mean(insarulos)

        # clean data
        # index = np.flatnonzero(np.logical_and(insar.uu>np.mean(insar.uu)+np.percentile(insar.uu,95),insar.uu<np.mean(insar.uu)-np.percentile(insar.uu,95)))
        # print np.percentile(insar.uu,35), np.percentile(insar.uu,65)
        # index = np.flatnonzero(np.logical_or(insar.uu<np.percentile(insar.uu,35),insar.uu>np.percentile(insar.uu,65)))
        # insar.uu,insar.xx,insar.yy,insar.xxpp,insar.yypp=np.delete(insar.uu,index),np.delete(insar.xx,index),\
        # np.delete(insar.yy,index),np.delete(insar.xxpp,index),np.delete(insar.yypp,index)
        
        # print insar.xxpp.max(),insar.xxpp.min()
        # print xpmax, xpmin

        bins = np.arange(-l/2,l/2,2)
        inds = np.digitize(insar.yypp,bins)
        insar.distance = []
        insar.moy_los = []
        insar.std_los = []
        insar.xperp = []
        insar.yperp = []
        insar.uulos =  []
 
        for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            # at leat 500 insar points
            # print len(uu)
            if len(uu)>50:
                insar.distance.append(bins[j] + (bins[j+1] - bins[j])/2.)

                # indice = np.flatnonzero(np.logical_and(insar.uu[uu]>np.percentile(\
                # 	insar.uu[uu],5.),insar.uu[uu]<np.percentile(insar.uu[uu],95.)))
                indice = np.flatnonzero(np.logical_and(insar.uu[uu]>np.percentile(\
                  insar.uu[uu],100-insar.perc),insar.uu[uu]<np.percentile(insar.uu[uu],insar.perc)))

                insar.std_los.append(np.std(insar.uu[uu][indice]))
                insar.moy_los.append(np.mean(insar.uu[uu][indice]))
                insar.xperp.append(insar.xxpp[uu][indice])
                insar.yperp.append(insar.yypp[uu][indice])
                insar.uulos.append(insar.uu[uu][indice])

                # insar.std_los.append(np.std(insar.uu[uu])) 
                # insar.moy_los.append(np.mean(insar.uu[uu]))
                # print len(uu), insar.moy_los[-1], insar.std_los[-1]
        insar.distance = np.array(insar.distance)
        insar.std_los = np.array(insar.std_los)
        insar.moy_los = np.array(insar.moy_los)
       	insar.xperp = np.concatenate(np.array(insar.xperp))
        insar.yperp = np.concatenate(np.array(insar.yperp))
        insar.uulos =  np.concatenate(np.array(insar.uulos))

        # PLOT
        if typ is 'distscale':
          # colorscale fct of the parrallel distance to the profile
          norm = matplotlib.colors.Normalize(vmin=xpmin, vmax=xpmax)
          m1 = cm.ScalarMappable(norm=norm,cmap='cubehelix_r')
          m1.set_array(insar.xperp)
          facelos=m1.to_rgba(insar.xperp)
          ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,\
          	 label=insardata[i].reduction,color=facelos, rasterized=True)
        
        elif typ is 'std':
          # plot mean and standard deviation
          # ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,\
          #    label=insardata[i].reduction,color=colors[i])
          ax2.plot(insar.distance,insar.moy_los,color=insar.color,lw=1.,label=insardata[i].reduction)
          ax2.plot(insar.distance,insar.moy_los-insar.std_los,color=insar.color,lw=0.5)
          ax2.plot(insar.distance,insar.moy_los+insar.std_los,color=insar.color,lw=0.5)

        elif typ is 'stdscat':
          # plot mean and standard deviation
          # ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,\
          #    label=insardata[i].reduction,color=colors[i])
          ax2.plot(insar.distance,insar.moy_los,color='black',lw=1.,label=insardata[i].reduction)
          ax2.plot(insar.distance,insar.moy_los-insar.std_los,color='black',lw=0.5)
          ax2.plot(insar.distance,insar.moy_los+insar.std_los,color='black',lw=0.5)
          ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,color=insar.color,rasterized=True)


        else:
          # plot scattering plot
          ax2.scatter(insar.yperp,insar.uulos,s = .1, marker='o',alpha=0.4,color=insar.color,rasterized=True)

        cst+=1.
        
        # set born profile equal to map
        if 'losmin' in locals():
          ax2.set_ylim([losmin,losmax])

        # for j in xrange(Mfault):
          # ax2.plot([fperp[j],fperp[j]],[losmax,losmin],color='red')
          # ax2.plot([fperp[j],fperp[j]],[losmax+cst,losmin-cst],color='red')
          # ax2.plot([fperp[j],fperp[j]],[6,-4],color='red')

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

fig1.savefig(outdir+profiles[k].name+'protopo.eps', format='EPS', dpi=150)
fig2.savefig(outdir+profiles[k].name+'prolos.pdf', format='PDF',dpi=150)
fig.savefig(outdir+profiles[k].name+'promap.eps', format='EPS',)
fig3.savefig(outdir+profiles[k].name+'progps.eps', format='EPS',)

plt.show()


