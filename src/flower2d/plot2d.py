from os import path, environ
import numpy as np
from numpy.lib.stride_tricks import as_strided
import matplotlib
if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as tic
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
from readgmt import *
from flatten import *
#import plot2Ddist
import math,sys,getopt
from os import path

def plotLOS(flt,nfigure):
    
    logger = flt.logger

    gpsdata, insardata = flt.gpsdata, flt.insardata
    fmodel = flt.fmodel
    profile = flt.profile
    Mseg = flt.Mseg
    Mstruc = flt.Mstruc
    Mvol = flt.Mvol
    volum = flt.volum
    topodata = flt.topodata
    plotdata = flt.plotdata
    seismifiles = flt.seismifiles
    outdir = flt.outdir

    l = profile.l
    w = profile.w
    x0 = profile.x
    y0 = profile.y
    name = profile.name
    xpmin,xpmax = profile.xpmin,profile.xpmax
    ypmin,ypmax = profile.ypmin,profile.ypmax
    strike = flt.strike
    s = flt.s
    n = flt.n

    proj = profile.proj
    
    fig = plt.figure(nfigure,figsize = (10,6))
    fig.subplots_adjust(hspace = 0.1)
    fig.subplots_adjust(hspace = 0.7)
    ax1 = fig.add_subplot(2,1,1)
    
    ax1.set_xlim([-l/2,l/2])
    tmin, tmax = 0, 0
    
    for j in xrange(len(topodata)):
        plot = topodata[j] 
        
        nb = np.float(l/(len(plot.z)/100.))
        logger.debug('Load {0}. Create bins every {1:.3f} km'.format(plot.name, nb)) 
        bins = np.arange(min(plot.yp),max(plot.yp),nb)
        inds = np.digitize(plot.yp,bins)
        distance = []
        moy_topo = []
        std_topo = []
        for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            if len(uu)>0:
                distance.append(bins[j] + (bins[j+1] - bins[j])/2.)
                std_topo.append(np.std(plot.z[uu]))
                moy_topo.append(np.mean(plot.z[uu]))

        distance = np.array(distance)
        std_topo = np.array(std_topo)
        moy_topo = np.array(moy_topo)
        max_topo,min_topo = np.array(moy_topo)+2*std_topo, np.array(moy_topo)-2*std_topo

        ax1.plot(distance,moy_topo,color=plot.color,lw=plot.width,label='topography')
        if plot.plotminmax == True:
          logger.info('plotminmax topo set to True')
          ax1.plot(distance,moy_topo-std_topo,color=plot.color,lw=int(plot.width/2))
          ax1.plot(distance,moy_topo+std_topo,color=plot.color,lw=int(plot.width/2))

        tmin, tmax = np.min(np.append(min_topo,tmin)), np.max(np.append(max_topo,tmax))

    if len(topodata) > 0:
        if (plot.topomin is not None) and (plot.topomax is not None) :
            logger.info('Set ylim topography to {} and {}'.format(plot.topomin,plot.topomax))
            ax1.set_ylim([plot.topomin,plot.topomax])
        else:
            ax1.set_ylim([np.nanmin(tmin)-0.5,np.nanmax(tmax)+0.5])

        plt.setp( ax1.get_xticklabels(), visible = False)
        ax1.set_ylabel('Elevation (km)')
        ax1.yaxis.set_major_locator(tic.MaxNLocator(3))
        ax1.legend(loc = 'best',fontsize='x-small')
        ax1.set_title('Profile %s, azimuth: %s'%(profile.name, strike))

    # number of plot data
    nb = int(flt.nsample/500) + 1
    
    #ax2 = fig.add_subplot(3+len(insardata),1,2)
    ax2 = fig.add_subplot(2,1,2)
    #ax2.axis('equal')
    ax2.set_xlim([-l/2,l/2])
    
    # plot decollement
    ax2.text(fmodel[0].fperp+1,-(fmodel[0].w+1),fmodel[0].name,color = 'black',style='italic',size='xx-small')
    ax2.text(fmodel[0].fperp+1,-(fmodel[0].w+2),'SS: %4.1f mm'%(fmodel[0].ss),style='italic',size='xx-small') 
    ax2.text(fmodel[0].fperp+1,-(fmodel[0].w+3),'DS: %4.1f mm'%(fmodel[0].ds),style='italic',size='xx-small') 
    
    ax2.scatter(fmodel[0].fperp,-fmodel[0].w, s = 30,marker = 'x',color = 'blue')
    
    if (fmodel[0].dip == 0 or fmodel[0].dip == 180):
        for i in xrange(0,len(fmodel[0].tracew),nb):
            ax2.plot([ fmodel[0].traceF[i], fmodel[0].traceF[i]],[ -fmodel[0].tracew[i], -100], '-' ,lw=.01, color='blue',zorder=20)
    else:
        for i in xrange(0,len(fmodel[0].tracew),nb):
            ax2.plot([ fmodel[0].traceF[i], fmodel[0].traceF[i] + fmodel[0].traceL[i]*np.cos(np.deg2rad(fmodel[0].tracedip[i]))], 
               [ -fmodel[0].tracew[i], -fmodel[0].tracew[i] - fmodel[0].traceL[i]*np.sin(np.deg2rad(fmodel[0].tracedip[i])) ], '-' ,lw=.01, color='blue',zorder=20)

    # plot ramp and kink 
    for k in xrange(1,flt.structures[0].Mseg):
        ax2.scatter(fmodel[k].fperp,-fmodel[k].w, s = 30,marker = 'x',color = 'royalblue')
        ax2.text(fmodel[k].fperp+1,-fmodel[k].w+1,fmodel[k].name,color = 'black',style='italic',size='xx-small')
        ax2.text(fmodel[k].fperp+1,-(fmodel[k].w+2),'SS: %4.1f mm'%((fmodel[k].ss)),style='italic',size='xx-small') 
        ax2.text(fmodel[k].fperp+1,-(fmodel[k].w+3),'DS: %4.1f mm'%((fmodel[k].ds)),style='italic',size='xx-small') 
        for i in xrange(0,len(fmodel[0].tracew),nb):
            ax2.plot([ fmodel[0].traceF[i], fmodel[k].traceF[i] ],[ -fmodel[0].tracew[i], -fmodel[k].tracew[i] ], '-' ,lw= .02, color='royalblue',zorder=20)
   
    Mtemp = flt.structures[0].Mseg
    for j in xrange(1,Mstruc):
        for k in xrange(flt.structures[j].Mseg):
            if fmodel[Mtemp+k].type != "creep" and fmodel[Mtemp+k].type != "vertical":
                ax2.scatter(fmodel[Mtemp+k].fperp,-fmodel[Mtemp+k].w, s = 30,marker = 'x',color = 'royalblue')
            else:
                ax2.scatter(fmodel[Mtemp+k].fperp,-fmodel[Mtemp+k].w, s = 30,marker = 'x',color = 'cyan')
            ax2.text(fmodel[Mtemp+k].fperp+1,-fmodel[Mtemp+k].w+1,fmodel[Mtemp+k].name,color = 'black',style='italic',size='xx-small')
            ax2.text(fmodel[Mtemp+k].fperp+1,-(fmodel[Mtemp+k].w+2),'SS: %4.1f mm'%(fmodel[Mtemp+k].ss),style='italic',size='xx-small') 
            ax2.text(fmodel[Mtemp+k].fperp+1,-(fmodel[Mtemp+k].w+3),'DS: %4.1f mm'%(fmodel[Mtemp+k].ds),style='italic',size='xx-small') 
            for i in xrange(0,len(fmodel[0].tracew),nb):
                if fmodel[Mtemp+k].type ==  "creep":
                  ax2.plot([ fmodel[0].traceF[i], fmodel[Mtemp+k].traceF[i] ],[ -fmodel[0].tracew[i], -fmodel[Mtemp+k].traceH[i] ], '-' ,lw= .02, color='cyan',zorder=20)  
                elif fmodel[Mtemp+k].type ==  "vertical":
                    ax2.plot([ fmodel[Mtemp+k].traceF[i] , fmodel[Mtemp+k].traceF[i] ],[ -fmodel[Mtemp+k].tracew[i], -100 ], '-' ,lw= .02, color='cyan',zorder=20)
                else:
                  ax2.plot([ fmodel[Mtemp+k].traceF[i] - fmodel[Mtemp+k].traceD[i] , fmodel[Mtemp+k].traceF[i] ],[ -(fmodel[Mtemp+k].tracew[i] + fmodel[Mtemp+k].traceH[i]), -fmodel[Mtemp+k].tracew[i] ], '-' ,lw= .02, color='royalblue',zorder=20)
        Mtemp += flt.structures[j].Mseg

    for i in xrange(len(plotdata)):
        plot = plotdata[i]
        ax2.scatter(plot.yp, -plot.z, s = plot.width, marker = 'o', linewidths=1, edgecolor='black', color = plot.color, label = plot.name, alpha=0.6,zorder=10)
   
    for ii in range(len(seismifiles)):
        color = seismifiles[ii].color
        smin = np.nanmin(seismifiles[ii].mag)
        width = (seismifiles[ii].mag - smin)*seismifiles[ii].width*5
        size = (seismifiles[ii].mag - smin)* float(seismifiles[ii].width)*5
        ax2.scatter(seismifiles[ii].yp,-seismifiles[ii].z,s=size,c=color,marker='o',linewidths=1, edgecolor='black',alpha=0.5,label=seismifiles[ii].name,zorder=10)
 
    l_arrow = float(profile.l)/20.
    w_arrow= l_arrow/5.
    if abs(fmodel[0].ds) > 0:
      # compressional arrow
      if (fmodel[0].ds > 0) and (fmodel[0].dip <= 90) or (fmodel[0].ds < 0) and (fmodel[0].dip > 90):
        x2_arrow = patches.FancyArrow(
        x=fmodel[0].fperp-(profile.l/4), y=-fmodel[0].w,
        dx=l_arrow, dy=0,
        width=w_arrow,
        head_length=w_arrow,
        head_width=l_arrow/2,
        length_includes_head=True,
        alpha=.8, fc='red')

        x1_arrow = patches.FancyArrow(
        x=fmodel[0].fperp+(profile.l/4), y=-fmodel[0].w,
        dx=-l_arrow, dy=0,
        width=w_arrow,
        head_length=w_arrow,
        head_width=l_arrow/2,
        length_includes_head=True,
        alpha=.8, fc='red')

        ax2.add_artist(x1_arrow)
        ax2.text(fmodel[0].fperp-(profile.l/4),-fmodel[0].w-2,'Short/Ext: %4.1f mm'%(fmodel[0].vh),style='italic',size='xx-small')
        ax2.add_artist(x2_arrow)

      # extensional arrow
      if (fmodel[0].ds > 0) and (fmodel[0].dip > 90) or (fmodel[0].ds < 0) and (fmodel[0].dip <= 90):
        x2_arrow = patches.FancyArrow(
        x=fmodel[0].fperp-(profile.l/4), y=-fmodel[0].w,
        dx=-l_arrow, dy=0,
        width=w_arrow,
        head_length=w_arrow,
        head_width=l_arrow/2,
        length_includes_head=True,
        alpha=.8, fc='red')

        x1_arrow = patches.FancyArrow(
        x=fmodel[0].fperp+(profile.l/4), y=-fmodel[0].w,
        dx=l_arrow, dy=0,
        width=w_arrow,
        head_length=w_arrow,
        head_width=l_arrow/2,
        length_includes_head=True,
        alpha=.8, fc='red')

        ax2.add_artist(x1_arrow)
        ax2.text(fmodel[0].fperp-(profile.l/4),-fmodel[0].w-2,'Short/Ext: %4.1f mm'%(fmodel[0].vh),style='italic',size='xx-small')
        ax2.add_artist(x2_arrow)


    plt.setp( ax2.get_xticklabels(), visible = False)
     
    if flt.depthmax == None:
        wmax = fmodel[0].w + 15
    else:
        wmax = flt.depthmax
    ax2.set_ylim([-wmax,5])
    ax2.set_ylabel('Depth (km)')
    ax2.yaxis.set_major_locator(tic.MaxNLocator(5))
    ax2.legend(loc = 'best',fontsize='x-small')
    ymin,ymax = ax2.get_ylim()
        
    # 1) Create synt data
    ysp = np.arange(-l/2,l/2,0.5)
    u = fmodel[0].displacement(ysp)
    for j in xrange(1,Mseg):
        u = u+fmodel[j].displacement(ysp)

    for j in xrange(Mvol):
        u = u+volum[j].displacement(ysp)

    # Sum of each faults
    uspar = u[:,0]
    usperp = u[:,1]
    usz = u[:,2]

    # Calculs x,y,z,los component for synt data
    usx = u[:,0]*s[0]+u[:,1]*n[0] 
    usy = u[:,0]*s[1]+u[:,1]*n[1]
    usz = u[:,2]

    # 2) Plot 
    markers = ['+','d','x','v']

    plt.tight_layout()
    logger.info('Save profile topographic figure in {}'.format(outdir+'/profile/'+name+'_topo.pdf'))
    np.savetxt(outdir+'/profile/'+name+'_topo.txt', np.vstack([distance,moy_topo]).T, header = '# distance (km)      z(km)', fmt='%.6f')
    fig.savefig(outdir+'/profile/'+name+'_topo.pdf', format = 'PDF')

    if len(insardata) < 3:
      fig = plt.figure(nfigure+1,figsize = (10,6))
    else:
      fig = plt.figure(nfigure+1,figsize = (10,9))
    fig.subplots_adjust(hspace = 0.1) 
    
    # fault parrallele
    ax4 = fig.add_subplot(1+len(insardata),1,1+len(insardata))

    ax4.set_xlim([-l/2,l/2])
    ax4.plot(ysp,uspar,'-',lw = 2.,color = 'blue',label='modeled fault-parallel displacements')
    # fault perpendicular
    ax4.plot(ysp,usperp,'-',lw = 2.,color = 'green',label='modeled fault-perpendicular displacements')
    ax4.plot(ysp,usz,'-',lw = 2.,color = 'red',label='modeled vertical displacements')
    ax4.set_ylabel('Displacements')
    ax4.yaxis.set_major_locator(tic.MaxNLocator(5))
    ax4.set_xlabel('Distance (km)')
    
    # 3) InSAR data and model
    colors = ['red','darkred','r','red','darkred','r','red','darkred','r']
    for i in xrange(len(insardata)):
        ax3 = fig.add_subplot(1+len(insardata),1,1+i)
        ax3.set_xlim([-l/2,l/2])
        plt.setp( ax3.get_xticklabels(), visible = False)
        ax3.set_ylabel('Velocities/Displacements')
        ax3.yaxis.set_major_locator(tic.MaxNLocator(3)) 
        
        insar = insardata[i]

        # Create Baseline component
        binsarlos = insar.a+insar.b*insar.yp
        ilos = insar.ulos-binsarlos

        ax3.scatter(insar.yp,insar.ulos-binsarlos,s = insar.width ,marker = 'o',label = insardata[i].reduction,color = insardata[i].color, alpha=.3)
        # problem if several insar with differents weight...
        

        if (insar.lmin != None) and (insar.lmax != None):
            ax3.set_ylim([insar.lmin,insar.lmax])

        uslos = usx*insar.projm[0]+usy*insar.projm[1]+usz*insar.projm[2]
        ax3.plot(ysp,uslos,'-',color='red', label = 'modeled {} displacements'.format(insar.reduction),lw = 2.)
        # sigmalos =  insar.sigmad[0]*np.ones(len(uslos))
        
        # plot legend
        ax3.legend(bbox_to_anchor = (0.,1.02,1.,0.102),loc = 3,ncol = 2,mode = 'expand',borderaxespad = 0.,fontsize = 'x-small')
  
    # 4) GPS data and model
    for i in xrange(len(gpsdata)):
        gps = gpsdata[i]
        markers = ['^','v','+','x']

        bpar = gps.a
        bperp = gps.b
        bv = gps.c
        
        # remove wieght in gps uncertainties
        #wd = gps.wd
        #gps.sigmapar, gps.sigmaperp = (gps.sigmapar/wd)*2 , (gps.sigmaperp/wd)*2 
        
        ax4.plot(gps.yp,gps.upar-bpar,markers[i],color = 'blue',mew = 1.,label = '%s fault-parallel displacements'%gpsdata[i].reduction )
        ax4.errorbar(gps.yp,gps.upar-bpar,yerr = gps.sigmapar,ecolor = 'blue',barsabove = 'True', fmt = 'none')

        ax4.plot(gps.yp,gps.uperp-bperp,markers[i],color = 'green',mew = 1.,label = '%s fault-perpendicular displacements'%gpsdata[i].reduction)
        ax4.errorbar(gps.yp,gps.uperp-bperp,yerr = gps.sigmaperp,ecolor = 'green',barsabove = 'True', fmt = 'none')
        #ymax,ymin = np.max(np.hstack([gps.upar-bpar+gps.sigmapar,gps.uperp-bperp+gps.sigmaperp]))+5,np.min(np.hstack([gps.upar-bpar-gps.sigmapar,gps.uperp-bperp-gps.sigmaperp]))-5
        ymax,ymin = np.max(np.hstack([gps.upar-bpar,gps.uperp-bperp]))+4,np.min(np.hstack([gps.upar-bpar,gps.uperp-bperp]))-4
        # ax4.set_ylim([ymin,ymax])
        
        if 3 == gps.dim:
            #gps.sigmav = (gps.sigmav/wd)*2
            bx = bpar*s[0]+bperp*n[0] 
            by = bpar*s[1]+bperp*n[1]
            # add gps in los if dim is 3
            if proj is not None:
                blos = bx*proj[0]+by*proj[1]+bv*proj[2]
                ax3.plot(gps.yp,gps.ulos-blos,'+',color='black',mew=1.,label='%s GPS LOS'%gpsdata[i].reduction)
            ax4.plot(gps.yp,gps.uv-bv,markers[i],color = 'red',mew = 1.,label = '%s vertical displacements'%gpsdata[i].reduction)
            ax4.errorbar(gps.yp,gps.uv-bv,yerr = gps.sigmav,ecolor = 'red',fmt = 'none')
            
            ymax,ymin = np.max(np.hstack([gps.upar-bpar,gps.uperp-bperp,gps.uv-bv]))+2,np.min(np.hstack([gps.upar-bpar,gps.uperp-bperp,gps.uv-bv]))-2
            ax4.set_ylim([ymin,ymax])

        if (gps.lmin != None) and (gps.lmax != None):
            ax4.set_ylim([gps.lmin,gps.lmax])
            ymax,ymin = gps.lmax,gps.lmin
    
    if len(insardata) > 0:
      # Plot fault and legend    
      ymin,ymax = ax3.get_ylim()
      for j in xrange(Mseg):
        ax3.plot([fmodel[j].fperp,fmodel[j].fperp],[ymin,ymax],'--',lw=1,color = 'black')
    
    ax4.legend(bbox_to_anchor = (0.,1.02,1.,0.102),loc = 3,ncol = 2,mode = 'expand',borderaxespad = 0.,fontsize ='x-small')
    ymin,ymax = ax4.get_ylim()
    for j in xrange(Mseg):
        ax4.plot([fmodel[j].fperp,fmodel[j].fperp],[ymin,ymax],'--',lw=1,color = 'black')
   
    plt.tight_layout()
    logger.info('Save profile LOS figure in {}'.format(outdir+'/profile/'+name+'_LOS.pdf'))
    fig.savefig(outdir+'/profile/'+name+'_LOS.pdf', format = 'PDF')
    
def plotMap(flt,nfigure):
    
    profile = flt.profile
    fmodel = flt.fmodel
    Mseg = flt.Mseg
    insardata = flt.insardata
    gpsdata = flt.gpsdata
    outdir = flt.outdir
    gmtfiles = flt.gmtfiles
    shapefiles = flt.shapefiles

    logger = flt.logger

    for i in xrange(len(insardata)):
        nn = i
        fig = plt.figure(nfigure,figsize = (12,5)) #width, height
        nfigure = nfigure + 1
        fig.subplots_adjust()
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2) 
        ax3 = fig.add_subplot(1,3,3) 
        markers = ['^','v','/^']
        colors = ['orange','m','yellow','red','blue','orange','m','yellow','red','blue']

        # bundary profile
        xp,yp = np.zeros((7)),np.zeros((7))
        x0 = profile.x
        y0 = profile.y
        l = profile.l
        w = profile.w
        name = profile.name
        xpmin,xpmax = profile.xpmin,profile.xpmax
        ypmin,ypmax = profile.ypmin,profile.ypmax
        
        strike = flt.strike
        proj = insardata[i].projm
        outdir = outdir

        s = flt.s
        n = flt.n

        xp[:] = x0-w/2*s[0]-l/2*n[0],x0+w/2*s[0]-l/2*n[0],x0+w/2*s[0]+l/2*n[0],x0-w/2*s[0]+l/2*n[0],x0-w/2*s[0]-l/2*n[0],x0-l/2*n[0],x0+l/2*n[0]
        yp[:] = y0-w/2*s[1]-l/2*n[1],y0+w/2*s[1]-l/2*n[1],y0+w/2*s[1]+l/2*n[1],y0-w/2*s[1]+l/2*n[1],y0-w/2*s[1]-l/2*n[1],y0-l/2*n[1],y0+l/2*n[1]
     
        logger.info('Save profile coordiantes in {}'.format(outdir+'/profile/'+name+'_coord.xy'))
        fid = open(outdir+'/profile/'+name+'_coord.xy','w')
        np.savetxt(fid, np.vstack([xp[:],yp[:]]).T ,header = 'x(km)     y(km) ',comments = '# ')
        fid.write('\n')
        fid.close

        ##InSAR DATA
        insar = insardata[i]
        imax = np.percentile(insar.ulos,99.5)
        kk = np.flatnonzero(insar.ulos>imax)
        insar.ulos[kk] = 0
        
        name = insar.reduction
        orb = insar.a+insar.b*insar.yp
        
        vmax = np.mean((insar.ulos)) + np.percentile(abs(insar.ulos-orb),98)
        vmin = np.mean((insar.ulos)) - np.percentile(abs(insar.ulos-orb),98)
        
        # cannot have different color scales
        norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
        m = cm.ScalarMappable(norm = norm, cmap = cm.jet)
        
        index = np.nonzero((insar.xp>xpmax)|(insar.xp<xpmin)|(insar.yp>ypmax)|(insar.yp<ypmin))
        insarx,insary,los,model = np.delete(insar.x,index),np.delete(insar.y,index),np.delete(insar.ulos,index),np.delete(insar.mlos+orb,index)

        logger.info('Write InSAR data in {} text file'.format(outdir+'/insar/'+name+'.xylos'))
        fid = open(outdir+'/insar/'+name+'.xylos','w')
        np.savetxt(fid, np.vstack([insarx, insary, los]).T ,header = 'x(km)     y(km)    los  ',comments = '# ')
        fid.write('\n')
        fid.close
        logger.info('Write InSAR model in {} text file'.format(outdir+'/insar/'+name+'model.xylos'))
        fid = open(outdir+'/insar/'+name+'model.xylos','w')
        np.savetxt(fid, np.vstack([insarx, insary, model]).T ,header = 'x(km)     y(km)    los  ',comments = '# ')
        fid.write('\n')
        fid.close
        logger.info('Write InSAR residual in {} text file'.format(outdir+'/insar/'+name+'residus.xylos'))
        fid=open(outdir+'/insar/'+name+'residus.xylos','w')
        np.savetxt(fid, np.vstack([insarx, insary, los-model]).T ,header='x(km)     y(km)    los  ',comments='# ')
        fid.write('\n')
        fid.close

        m.set_array(los)
        
        # InSAR model 
        facelos = m.to_rgba(los)
        facemodel = m.to_rgba(model)
        faceres = m.to_rgba(los-model) 

        ax1.scatter(insarx,insary,s = insar.width*3, marker = 'o',color = facelos,label = insar.reduction)
        ax2.scatter(insarx,insary,s = insar.width*3, marker = 'o',color = facemodel)
        ax3.scatter(insarx,insary,s = insar.width*3, marker = 'o',color = faceres)

        for j in xrange(len(gpsdata)):
            gps = gpsdata[j]
            name = gps.reduction

            index = np.nonzero((gps.xp>xpmax)|(gps.xp<xpmin)|(gps.yp>ypmax)|(gps.yp<ypmin))        
            if gps.dim==2:
                gpsname,gpsx,gpsy,gpsux,gpsuy,gpsmx,gpsmy,gpssigmax,gpssigmay = np.delete(gps.name,index),np.delete(gps.x,index),np.delete(gps.y,index),np.delete(gps.ux,index),np.delete(gps.uy,index),np.delete(gps.mx,index),np.delete(gps.my,index),np.delete(gps.sigmax,index),np.delete(gps.sigmay,index)
            else:
                gpsname,gpsx,gpsy,gpsux,gpsuy,gpsuv,gpsmx,gpsmy,gpsmz,gpssigmax,gpssigmay,gpssigmav = np.delete(gps.name,index),np.delete(gps.x,index),np.delete(gps.y,index),np.delete(gps.ux,index),np.delete(gps.uy,index),np.delete(gps.uv,index),np.delete(gps.mx,index),np.delete(gps.my,index),np.delete(gps.mz,index),np.delete(gps.sigmax,index),np.delete(gps.sigmay,index),np.delete(gps.sigmav,index)
            # remove wieght in gps uncertainties
            wd = gps.wd
            gpssigmax,gpssigmay = gpssigmax/wd , gpssigmay/wd
            
            bx = (gps.a*flt.s[0]+gps.b*flt.n[0])
            by = (gps.a*flt.s[1]+gps.b*flt.n[1])
            
            if gps.dim==2:
            
                logger.info('Write GPS data in {} text file'.format(outdir+'/gps/'+name+'_%s.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux,gpsuy,gpssigmax,gpssigmay]).T ,header = 'x(km)  y(km)   East_Vel    North_Vel   East_StdDev North_StdDev',comments = '# ')
                fid.write('\n')
                fid.close
                logger.info('Write GPS model in {} text file'.format(outdir+'/gps/'+name+'model_%s.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'model_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsmx,gpsmy,gpssigmax,gpssigmay]).T ,header = 'x(km)  y(km)   East_Vel    North_Vel   East_StdDev North_StdDev',comments = '# ')
                fid.write('\n')
                fid.close    
                logger.info('Write GPS model without optimised baselines in {} text file'.format(outdir+'/gps/'+name+'residus_%s.psvelo'%(j)))
                fid=open(outdir+'/gps/'+name+'residus_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy,gpsux-gpsmx-bx,gpsuy-gpsmy-by,gpssigmax,gpssigmay]).T ,header='x(km)  y(km)   East_Dev    North_Dev   East_StdDev North_StdDev',comments='# ')
                fid.write('\n')
                fid.close
                logger.info('Write GPS residual in {} text file'.format(outdir+'/gps/'+name+'_%s_rf.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'_%s_rf.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux-bx,gpsuy-by,gpssigmax,gpssigmay]).T ,header = 'x(km)   y(km)   East_Vel    North_Vel   East_StdDev North_StdDev',comments = '# ')
                fid.write('\n')
                fid.close
            
            else:
                
                bv=gps.c
                logger.info('Write GPS data in {} text file'.format(outdir+'/gps/'+name+'_%s.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux,gpsuy,gpsuv,gpssigmax,gpssigmay,gpssigmav]).T ,header = 'x(km)  y(km)   East_Vel    North_Vel Up_Vel  East_StdDev North_StdDev Up_StdDev',comments = '# ')
                fid.write('\n')
                fid.close
                logger.info('Write GPS model in {} text file'.format(outdir+'/gps/'+name+'model_%s.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'model_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsmx,gpsmy,gpsmz,gpssigmax,gpssigmay,gpssigmav]).T ,header = 'x(km)  y(km)   East_Vel    North_Vel Up_Vel  East_StdDev North_StdDev Up_StdDev',comments = '# ')
                fid.write('\n')
                fid.close  
                logger.info('Write GPS model without optimised baselines in {} text file'.format(outdir+'/gps/'+name+'_%s_rf.psvelo'%(j)))  
                fid = open(outdir+'/gps/'+name+'_%s_rf.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux-bx,gpsuy-by,gpsuv-bv,gpssigmax,gpssigmay,gpssigmav]).T ,header = 'x(km)   y(km)   East_Vel    North_Vel  Up_Vel  East_StdDev North_StdDev Up_StdDev',comments = '# ')
                fid.write('\n')
                fid.close
                logger.info('Write GPS residual in {} text file'.format(outdir+'/gps/'+name+'residus_%s.psvelo'%(j)))
                fid=open(outdir+'/gps/'+name+'residus_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy,gpsux-gpsmx-bx,gpsuy-gpsmy-by,gpsuv-gpsmz-bv,gpssigmax,gpssigmay,gpssigmav]).T ,header='x(km)  y(km)   East_Dev    North_Dev  Up_Dev   East_StdDev North_StdDev  Up_StdDev',comments='# ')
                fid.write('\n')
                fid.close
            
            data = ax1.quiver(gpsx,gpsy,gpsux-bx,gpsuy-by,scale = 10,width = 0.01,color = 'black')
            model = ax2.quiver(gpsx,gpsy,gpsmx,gpsmy,scale = 10,width = 0.01,color = 'red')
            residus = plt.quiver(gpsx,gpsy,gpsux-gpsmx-bx,gpsuy-gpsmy-by,scale = 10,width = 0.01,color = 'black')
                    
            ## display name of the station
            if gps.plotName is True:
                for kk in xrange(len(gpsname)):
                    ax1.text(gpsx[kk]-4*kk,gpsy[kk]-10,gpsname[kk],color = 'black')
            ax1.quiverkey(data,0.1,1.015,20.,'GPS displacements',coordinates = 'axes',color = 'black')
            ax2.quiverkey(model,0.1,1.015,20.,'Model',coordinates = 'axes',color = 'red')
            ax3.quiverkey(model,1.3,1.015,20,'Residuals',coordinates = 'axes',color = 'red')

        ylim, xlim = ax1.get_ylim(), ax1.get_xlim()
        
        xf,yf = np.zeros((Mseg,2)),np.zeros((Mseg,2))
        for j in xrange(Mseg):
            xf[j,0] = fmodel[j].x-w*s[0]
            xf[j,1] = fmodel[j].x+w*s[0]
            yf[j,0] = fmodel[j].y-w*s[1]
            yf[j,1] = fmodel[j].y+w*s[1]
        
        axes = [ax1, ax2, ax3]
        titles = ['Data', 'Model', 'Residual']
        for ax, title in zip(axes,titles):
            ax.axis('equal')
            ax.legend(loc='best',fontsize='x-small')
            for ii in xrange(len(gmtfiles)):
                        name = gmtfiles[ii].name
                        wdir = gmtfiles[ii].wdir
                        filename = gmtfiles[ii].filename
                        color = gmtfiles[ii].color
                        width = gmtfiles[ii].width
                        fx,fy = gmtfiles[ii].load()
                        for i in xrange(len(fx)):
                            ax.plot(fx[i],fy[i],color = color,lw = width)

            for ii in xrange(len(shapefiles)):
                try:
                    import geopandas as gpd
                    import shapely.speedups
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
                        shape.plot(ax=ax,facecolor='none', color=color,edgecolor=edgecolor,lw=linewidth)
                except:
                        logger.warning('Cannot import geopandas and shapely. Pass shapefiles plot')    
            
            ax.plot(xp[:],yp[:],color = 'black',lw = 1.)
            for f in xrange(Mseg):
                ax.plot(xf[f,:],yf[f,:],'--',color = 'black',lw = 1.)
            
            # plot LOS
            h = xlim[1] - xlim[0]
            x1, y1 = xlim[1] - int(h/5) , ylim[0] + (ylim[1]-ylim[0])/5

            if proj is not None:
                l_arrow = np.float(w)/2
                w_arrow = np.sqrt((proj[0]*l_arrow)**2+(proj[1]*l_arrow)**2)/8.
                los_arrow = patches.FancyArrow(
                        x=x1, y=y1,
                        dx=proj[0]*l_arrow, dy=proj[1]*l_arrow,
                        width=w_arrow,
                        head_length=w_arrow*2.,
                        head_width=w_arrow*2.,
                        alpha=.8, fc='k',
                        length_includes_head=True)

                ax.add_artist(los_arrow)
                # ax.arrow(x1,y1,30*proj[0],30*proj[1],fill = True,width = 0.5,color = 'black')
                # ax.text(x1-1,y1-1,'LOS',color = 'black',alpha=.8)

            ax.set_xlabel('Distance (km)')
            ax.set_title(title)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            if len(insardata)>0:
                divider = make_axes_locatable(ax)
                c = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(m, cax=c)
                # fig.colorbar(m,shrink = 0.5, aspect = 5)

        plt.suptitle('Geodetic map')
        logger.info('Save map figure in {}'.format(outdir+'/map/'+profile.name+'_map_%s.pdf'%(nn)))
        fig.savefig(outdir+'/map/'+profile.name+'_map_%s.pdf'%(nn),format = 'PDF')

    if (len(insardata) == 0) and (len(gpsdata)>0):
      for j in xrange(len(gpsdata)):
        nn = j
        fig = plt.figure(nfigure,figsize = (12,5)) #width, height
        nfigure = nfigure + 1
        fig.subplots_adjust()
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2) 
        ax3 = fig.add_subplot(1,3,3) 
        markers = ['^','v','/^']
        colors = ['orange','m','yellow','red','blue','orange','m','yellow','red','blue']

        # bundary profile
        xp,yp = np.zeros((7)),np.zeros((7))
        x0 = profile.x
        y0 = profile.y
        l = profile.l
        w = profile.w
        name = profile.name
        xpmin,xpmax = profile.xpmin,profile.xpmax
        ypmin,ypmax = profile.ypmin,profile.ypmax
        
        strike = flt.strike
        proj = profile.proj
        outdir = outdir

        s = flt.s
        n = flt.n

        xp[:] = x0-w/2*s[0]-l/2*n[0],x0+w/2*s[0]-l/2*n[0],x0+w/2*s[0]+l/2*n[0],x0-w/2*s[0]+l/2*n[0],x0-w/2*s[0]-l/2*n[0],x0-l/2*n[0],x0+l/2*n[0]
        yp[:] = y0-w/2*s[1]-l/2*n[1],y0+w/2*s[1]-l/2*n[1],y0+w/2*s[1]+l/2*n[1],y0-w/2*s[1]+l/2*n[1],y0-w/2*s[1]-l/2*n[1],y0-l/2*n[1],y0+l/2*n[1]
     
        logger.info('Save profile coordiantes in {}'.format(outdir+'/profile/'+name+'_coord.xy'))
        fid = open(outdir+'/profile/'+name+'_coord.xy','w')
        np.savetxt(fid, np.vstack([xp[:],yp[:]]).T ,header = 'x(km)     y(km) ',comments = '# ')
        fid.write('\n')
        fid.close

        gps = gpsdata[j]
        name = gps.reduction

        index = np.nonzero((gps.xp>xpmax)|(gps.xp<xpmin)|(gps.yp>ypmax)|(gps.yp<ypmin))        
        if gps.dim==2:
                gpsname,gpsx,gpsy,gpsux,gpsuy,gpsmx,gpsmy,gpssigmax,gpssigmay = np.delete(gps.name,index),np.delete(gps.x,index),np.delete(gps.y,index),np.delete(gps.ux,index),np.delete(gps.uy,index),np.delete(gps.mx,index),np.delete(gps.my,index),np.delete(gps.sigmax,index),np.delete(gps.sigmay,index)
        else:
                gpsname,gpsx,gpsy,gpsux,gpsuy,gpsuv,gpsmx,gpsmy,gpsmz,gpssigmax,gpssigmay,gpssigmav = np.delete(gps.name,index),np.delete(gps.x,index),np.delete(gps.y,index),np.delete(gps.ux,index),np.delete(gps.uy,index),np.delete(gps.uv,index),np.delete(gps.mx,index),np.delete(gps.my,index),np.delete(gps.mz,index),np.delete(gps.sigmax,index),np.delete(gps.sigmay,index),np.delete(gps.sigmav,index)
        # remove wieght in gps uncertainties
        wd = gps.wd
        gpssigmax,gpssigmay = gpssigmax/wd , gpssigmay/wd
            
        bx = (gps.a*flt.s[0]+gps.b*flt.n[0])
        by = (gps.a*flt.s[1]+gps.b*flt.n[1])
            
        if gps.dim==2:
            
                logger.info('Write GPS data in {} text file'.format(outdir+'/gps/'+name+'_%s.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux,gpsuy,gpssigmax,gpssigmay]).T ,header = 'x(km)  y(km)   East_Vel    North_Vel   East_StdDev North_StdDev',comments = '# ')
                fid.write('\n')
                fid.close
                logger.info('Write GPS model in {} text file'.format(outdir+'/gps/'+name+'model_%s.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'model_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsmx,gpsmy,gpssigmax,gpssigmay]).T ,header = 'x(km)  y(km)   East_Vel    North_Vel   East_StdDev North_StdDev',comments = '# ')
                fid.write('\n')
                fid.close    
                logger.info('Write GPS model without optimised baselines in {} text file'.format(outdir+'/gps/'+name+'residus_%s.psvelo'%(j)))
                fid=open(outdir+'/gps/'+name+'residus_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy,gpsux-gpsmx-bx,gpsuy-gpsmy-by,gpssigmax,gpssigmay]).T ,header='x(km)  y(km)   East_Dev    North_Dev   East_StdDev North_StdDev',comments='# ')
                fid.write('\n')
                fid.close
                logger.info('Write GPS residual in {} text file'.format(outdir+'/gps/'+name+'_%s_rf.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'_%s_rf.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux-bx,gpsuy-by,gpssigmax,gpssigmay]).T ,header = 'x(km)   y(km)   East_Vel    North_Vel   East_StdDev North_StdDev',comments = '# ')
                fid.write('\n')
                fid.close
            
        else:
                
                bv=gps.c
                logger.info('Write GPS data in {} text file'.format(outdir+'/gps/'+name+'_%s.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux,gpsuy,gpsuv,gpssigmax,gpssigmay,gpssigmav]).T ,header = 'x(km)  y(km)   East_Vel    North_Vel Up_Vel  East_StdDev North_StdDev Up_StdDev',comments = '# ')
                fid.write('\n')
                fid.close
                logger.info('Write GPS model in {} text file'.format(outdir+'/gps/'+name+'model_%s.psvelo'%(j)))
                fid = open(outdir+'/gps/'+name+'model_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsmx,gpsmy,gpsmz,gpssigmax,gpssigmay,gpssigmav]).T ,header = 'x(km)  y(km)   East_Vel    North_Vel Up_Vel  East_StdDev North_StdDev Up_StdDev',comments = '# ')
                fid.write('\n')
                fid.close  
                logger.info('Write GPS model without optimised baselines in {} text file'.format(outdir+'/gps/'+name+'_%s_rf.psvelo'%(j)))  
                fid = open(outdir+'/gps/'+name+'_%s_rf.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux-bx,gpsuy-by,gpsuv-bv,gpssigmax,gpssigmay,gpssigmav]).T ,header = 'x(km)   y(km)   East_Vel    North_Vel  Up_Vel  East_StdDev North_StdDev Up_StdDev',comments = '# ')
                fid.write('\n')
                fid.close
                logger.info('Write GPS residual in {} text file'.format(outdir+'/gps/'+name+'residus_%s.psvelo'%(j)))
                fid=open(outdir+'/gps/'+name+'residus_%s.psvelo'%(j),'w')
                np.savetxt(fid, np.vstack([gpsx, gpsy,gpsux-gpsmx-bx,gpsuy-gpsmy-by,gpsuv-gpsmz-bv,gpssigmax,gpssigmay,gpssigmav]).T ,header='x(km)  y(km)   East_Dev    North_Dev  Up_Dev   East_StdDev North_StdDev  Up_StdDev',comments='# ')
                fid.write('\n')
                fid.close
            
        data = ax1.quiver(gpsx,gpsy,gpsux-bx,gpsuy-by,scale = 10,width = 0.01,color = 'black')
        model = ax2.quiver(gpsx,gpsy,gpsmx,gpsmy,scale = 10,width = 0.01,color = 'red')
        residus = plt.quiver(gpsx,gpsy,gpsux-gpsmx-bx,gpsuy-gpsmy-by,scale = 10,width = 0.01,color = 'black')
                    
        ## display name of the station
        if gps.plotName is True:
                for kk in xrange(len(gpsname)):
                    ax1.text(gpsx[kk]-4*kk,gpsy[kk]-10,gpsname[kk],color = 'black')
        #ax1.quiverkey(data,0.1,1.015,20.,'GPS displacements',coordinates = 'axes',color = 'black')
        #ax2.quiverkey(model,0.1,1.015,20.,'Model',coordinates = 'axes',color = 'red')
        #ax3.quiverkey(model,1.3,1.015,20,'Residuals',coordinates = 'axes',color = 'red')
        #ylim, xlim = ax1.get_ylim(), ax1.get_xlim()
        
        xf,yf = np.zeros((Mseg,2)),np.zeros((Mseg,2))
        for j in xrange(Mseg):
            xf[j,0] = fmodel[j].x+-w*s[0]
            xf[j,1] = fmodel[j].x+w*s[0]
            yf[j,0] = fmodel[j].y+-w*s[1]
            yf[j,1] = fmodel[j].y+w*s[1]
        
        axes = [ax1, ax2, ax3]
        titles = ['Data', 'Model', 'Residual']
        for ax, title in zip(axes,titles):
            ax.axis('equal')
            ax.legend(loc='best',fontsize='x-small')
            for ii in xrange(len(gmtfiles)):
                        name = gmtfiles[ii].name
                        wdir = gmtfiles[ii].wdir
                        filename = gmtfiles[ii].filename
                        color = gmtfiles[ii].color
                        width = gmtfiles[ii].width
                        fx,fy = gmtfiles[ii].load()
                        for i in xrange(len(fx)):
                            ax.plot(fx[i],fy[i],color = color,lw = width)

            for ii in xrange(len(shapefiles)):
                try:
                    import geopandas as gpd
                    import shapely.speedups
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
                        shape.plot(ax=ax,facecolor='none', color=color,edgecolor=edgecolor,linewidth=linewidth)
                except:
                        logger.warning('Cannot import geopandas and shapely. Pass shapefiles plot')    
            

            ax.plot(xp[:],yp[:],color = 'black',lw = 1.)
            for f in xrange(Mseg):
                ax.plot(xf[f,:],yf[f,:],'--',color = 'black',lw = 1.)
            
            ax.set_xlabel('Distance (km)')
            ax.set_title(title)


def plotHist(flt,model,nfigure):
   
    from scipy.stats import norm

    Mseg = flt.Mseg
    fmodel = flt.fmodel
    Sampled = flt.Sampled
    mu, tau = flt.mu, flt.tau
    mmin, mmax = flt.mmin, flt.mmax
    profile = flt.profile
    outdir = flt.outdir
    logger = flt.logger

    # for histo
    traces = []
    labels = []
   
    if len(Sampled) < 10:
        fig = plt.figure(nfigure, figsize = (10,6))
    else:
        fig = plt.figure(nfigure, figsize = (10,8))
    fig.subplots_adjust(hspace = 0.75,wspace = 0.35)
   
    def histo(traces,labels,Mseg,Mker,l,mu,tau,mmin,mmax):
        ll = 1
        for trace,label,mu,tau,mmin,mmax in zip(traces,labels,mu,tau,mmin,mmax):    
            index = l*Mker + ll
            ax = fig.add_subplot(Mseg,Mker,index)
           
            # plot Posterior model
            histo = ax.hist(trace, bins=50, density = True, histtype = 'stepfilled', color = 'dodgerblue', label = label, alpha=.6)
            mean = trace.mean()
            sig = np.std(trace)
            #plt.axvline(x = mean, c = "red")
            #plt.axvline(x = mean+sig, c = "dodgerblue", linestyle = '--')
            #plt.axvline(x = mean-sig, c = "dodgerblue", linestyle = '--')
            # xmin,xmax = ax.set_xlim(mmin,mmax)
            #ymin,ymax = ax.get_ylim()
            ax.set_title(label)
            ax.set_xlabel("Mean: {0:.2f}+-{1:.2f} (2*std)".format(mean, 2*sig))
            ax.set_ylabel('Norm. PDF')
            ax.set_ylim(bottom=0)

            #f = ax.lines[0]
            #xf = f.get_xydata()[:,0]
            #yf = f.get_xydata()[:,1]
            #ax.fill_between(xf, yf, color="dodgerblue", alpha=0.5, where=(xf>(mean-1*sig)) & (xf<(mean+1*std)))
            #ax1.axvline(xf[np.argmax(yf)], color="dodgerblue",linestyle='dashed')

            ll +=  1
   
    Mker = 0 
    if '{} Strike Slip'.format(fmodel[0].name) in Sampled:
        strike = model.trace('{} Strike Slip'.format(fmodel[0].name))[:]
        traces.append(strike)
        labels.append('{} Strike Slip'.format(fmodel[0].name))
        Mker = Mker + 1
    if '{} Shortening'.format(fmodel[0].name) in Sampled:
        short = model.trace('{} Shortening'.format(fmodel[0].name))[:]
        traces.append(short)
        labels.append('{} Shortening'.format(fmodel[0].name))
        Mker = Mker + 1
    if '{} H'.format(fmodel[0].name) in Sampled:
        H1 = model.trace('{} H'.format(fmodel[0].name))[:]
        traces.append(np.ones((flt.nsample))*flt.fmodel[0].winit - H1)
        labels.append('{} Depth'.format(fmodel[0].name))
        Mker = Mker + 1
    if '{} D'.format(fmodel[0].name) in Sampled:
        D1 = model.trace('{} D'.format(fmodel[0].name))[:]
        traces.append(D1)
        labels.append('{} D'.format(fmodel[0].name))
        Mker = Mker + 1
    if '{} L'.format(fmodel[0].name) in Sampled:
        L1 = model.trace('{} L'.format(fmodel[0].name))[:]
        traces.append(L1)
        labels.append('{} Length'.format(fmodel[0].name))
        Mker = Mker + 1
    if '{} dip'.format(fmodel[0].name) in Sampled:
        dip1 = model.trace('{} dip'.format(fmodel[0].name))[:]
        traces.append(dip1)
        labels.append('{} dip'.format(fmodel[0].name))
        Mker = Mker + 1

    j = 0
    # plot histos main fault
    histo(traces,labels,Mseg,Mker,j,mu[:Mker],tau[:Mker],mmin[:Mker],mmax[:Mker])
   
    M = Mker
    for j in xrange(1,Mseg):
        traces = []
        labels = []
        if '{} Strike Slip'.format(fmodel[j].name) in Sampled:
            m = model.trace('{} Strike Slip'.format(fmodel[j].name))[:]
            traces.append(m)
            labels.append('{} Strike Slip'.format(fmodel[j].name))
            
        if '{} D'.format(fmodel[j].name) in Sampled:
            m = model.trace('{} D'.format(fmodel[j].name))[:]
            traces.append(m)
            labels.append('{} D'.format(fmodel[j].name))

        if '{} H'.format(fmodel[j].name) in Sampled:
            m = model.trace('{} H'.format(fmodel[j].name))[:]
            traces.append(m)
            labels.append('{} H'.format(fmodel[j].name))
        
        histo(traces,labels,Mseg,Mker,j,mu[M:M+fmodel[j].Mker],tau[M:M+fmodel[j].Mker],mmin[M:M+fmodel[j].Mker],mmax[M:M+fmodel[j].Mker])
        M = M + fmodel[j].Mker 
    
    # save
    plt.tight_layout()
    logger.info('Save statistic in {}'.format(outdir+'/stat/'+profile.name+'_histo.pdf'))
    fig.savefig(outdir+'/stat/'+profile.name+'_histo.pdf', format = 'PDF')
    

def plotDist(flt,model,nfigure):

    from scipy.stats import norm
    import plot2Ddist
    scatterstyle={'color':'r', 'alpha':0.1}
    styleargs = {'color':'k', 'scatterstyle':scatterstyle}

    fmodel=flt.fmodel
    profile=flt.profile
    Mseg=len(fmodel)
    outdir = flt.outdir
    Mseg = flt.Mseg
    Mstruc = flt.Mstruc
    Mker = fmodel[0].Mker
    Sampled = flt.Sampled
    nsample = flt.nsample
    
    thin = 1
    rad2deg = 180/np.pi
    
    if '{} Strike Slip'.format(fmodel[0].name) in Sampled:
        Us1 = abs(model.trace('{} Strike Slip'.format(fmodel[0].name))[:])
    else:
        Us1 = abs(fmodel[0].ss)*np.ones((nsample)) 
    if '{} Dip Slip'.format(fmodel[0].name) in Sampled:
        Ushort = abs(model.trace('{} Dip Slip'.format(fmodel[0].name))[:])
    else:
        Ushort = abs(fmodel[0].ds)*np.ones((nsample)) 
    if '{} Depth'.format(fmodel[0].name) in Sampled:
        H1 = abs(model.trace('{} Depth'.format(fmodel[0].name))[:])
    else:
        H1 = abs(fmodel[0].H)*np.ones((nsample)) 

    # joint Pdf SS
    name=outdir+'/stat/'+profile.name+'_jointPDFs-Us1.pdf'
    if (np.std(Us1) > 0) and (np.std(Ushort) > 0) and (np.std(H1) > 0):
        plot2Ddist.plot2DdistsPairs(Us1,[Ushort,H1], mainlabel='Us1', labels = ['Ushort','H1'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)
    if (np.std(Us1) > 0) and (np.std(Ushort) > 0):
        plot2Ddist.plot2DdistsPairs(Us1,[Ushort], mainlabel='Us1', labels = ['Ushort'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)
    if (np.std(Us1) > 0) and (np.std(H1) > 0):
        plot2Ddist.plot2DdistsPairs(Us1,[H1], mainlabel='Us1', labels = ['H1'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)

    # others segments
    for j in xrange(1,Mseg): # probeles if more than 2 structures
        if '{} Strike Slip'.format(fmodel[j].name) in Sampled:
            Us2 = abs(model.trace('{} Strike Slip'.format(fmodel[j].name))[:])
        else:
            Us2 = abs(fmodel[j].ss)*np.ones((nsample))
        if '{} Width'.format(fmodel[j].name) in Sampled:
            D2 = model.trace('{} Width'.format(fmodel[j].name))[:]
        else:
            D2 = fmodel[j].D*np.ones((nsample))
        if '{} Depth'.format(fmodel[j].name) in Sampled:
            H2 = abs(model.trace('{} Depth'.format(fmodel[j].name))[:])
        else:
            H2 = abs(fmodel[j].H)*np.ones((nsample))
            
        if np.mean(D2) > 1:
            L = np.sqrt((H1-H2)**2+D2**2)
            fmodel[j].sigmaL = np.std(L)
            alpha = np.arctan2((H1-H2),-D2)
            fmodel[j].sigmadip = np.std(alpha)*rad2deg
            Ud = (Ushort/np.cos(alpha))
            #dinf = np.percentile(Ud,90)
            #print dinf
            #kk = np.flatnonzero(abs(Ud)>=abs(dinf))
            #Ud[kk]=0
            fmodel[j].sigmads = np.std(Ud)

            # joint Pdf Us2
            name=outdir+'/stat/'+profile.name+'_jointPDFs-Us2.pdf'
            if (np.std(Us2) > 0) and (np.std(H2) > 0 ) and (np.std(Us1) > 0) and  (np.std(D2) > 0):
                plot2Ddist.plot2DdistsPairs(Us2,[H2,Us1,D2], mainlabel='Us2', labels = ['H2','Us1','D2'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)
            elif (np.std(Us2) > 0) and (np.std(H2) > 0 ) and (np.std(Us1) > 0):
                plot2Ddist.plot2DdistsPairs(Us2,[H2,Us1], mainlabel='Us2', labels = ['H2','Us1'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)
            elif (np.std(Us2) > 0) and (np.std(H2) > 0):
                plot2Ddist.plot2DdistsPairs(Us2,[H2], mainlabel='Us2', labels = ['H2'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)
            
            # joint Pdf Ud
            name=outdir+'/stat/'+profile.name+'_jointPDFs-Ud.pdf'
            if (np.std(Ud) > 0) and (np.std(alpha) > 0):
                plot2Ddist.plot2DdistsPairs(Ud,[alpha*rad2deg], mainlabel='Ud', labels = ['Dip angle'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)

        else :
            Ud = np.zeros((nsample))
            fmodel[j].sigmads = 0
            fmodel[j].sigmadip = 0
            
            # joint Pdf Uc
            name=outdir+'/stat/'+profile.name+'_jointPDFs-Uc.pdf'
            if (np.std(Us2) > 0) and (np.std(H2) > 0 ) and (np.std(Us1) > 0):
                plot2Ddist.plot2DdistsPairs(Us2,[H2,Us1], mainlabel='Uc', labels = ['Hc','Us1'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)
            elif (np.std(Us2) > 0) and (np.std(H2) > 0):
                plot2Ddist.plot2DdistsPairs(Us2,[H2], mainlabel='Uc', labels = ['Hc'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)

    uu=1
    for i in xrange(len(manifolds)):
        if 3==manifolds[i].dim:
            a = model.trace('Baseline {}'.format(uu))[:]
            b = model.trace('Baseline {}'.format(uu+1))[:]
            c = model.trace('Baseline {}'.format(uu+2))[:]
            uu+=3
        elif 2==manifolds[i].dim:
            a = model.trace('Baseline {}'.format(uu))[:]
            b = model.trace('Baseline {}'.format(uu+1))[:]
            uu+=2
        else:
            a = model.trace('Baseline {}'.format(uu))[:]
            b = model.trace('Baseline {}'.format(uu+1))[:]
            uu+=1
            # joint Pdf ramp
            name=outdir+'/stat/'+profile.name+'_jointPDFs-ramp.pdf'
            if (np.std(b) > 0) and  (np.std(Us1) > 0) and (np.std(Us2) > 0 ) and (np.std(Ud) > 0):
                plot2Ddist.plot2DdistsPairs(b,[Us1,Us2,Ud], mainlabel='InSAR ramp', labels = ['Us1','Us2','Ud'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)
            if (np.std(b) > 0) and (np.std(Us1) > 0) and (np.std(Ud) > 0):
                plot2Ddist.plot2DdistsPairs(b,[Us1,Ud], mainlabel='InSAR ramp', labels = ['Us1','Ud'],plotcontours=True,thin=thin,contourKDEthin=10,histbinslist=[30,30],out=name,  **styleargs)

