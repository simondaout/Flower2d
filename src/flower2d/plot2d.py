from os import path, environ
import numpy as np
from numpy.lib.stride_tricks import as_strided

import matplotlib
if environ["TERM"].startswith("screen"):
    matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as tic
from matplotlib.backends.backend_pdf import PdfPages

from readgmt import *
from flatten import *
#import plot2Ddist

import math,sys,getopt
from os import path

def plotLOS(flt,nfigure):
    #fig = plt.figure(nfigure,figsize = (14,15))
    fig = plt.figure(nfigure,figsize = (14,10))
    fig.subplots_adjust(hspace = 0.1)

    gpsdata, insardata = flt.gpsdata, flt.insardata
    fmodel = flt.fmodel
    profiles = flt.profiles
    Mseg = flt.Mseg
    Mstruc = flt.Mstruc
    Mvol = flt.Mvol
    volum = flt.volum
    topodata = flt.topodata
    plotdata = flt.plotdata
    outdir = flt.outdir

    l = profiles.l
    w = profiles.w
    x0 = profiles.x
    y0 = profiles.y
    name = profiles.name
    xp0 = profiles.xp0
    yp0 = profiles.yp0
    xpmin,xpmax = profiles.xpmin,profiles.xpmax
    ypmin,ypmax = profiles.ypmin,profiles.ypmax
    strike = flt.strike
    s = flt.s
    n = flt.n

    proj = profiles.proj

    fig.subplots_adjust(hspace = 0.7)
    ax1 = fig.add_subplot(4,1,1)
    
    ax1.set_xlim([-l/2,l/2])
    tmin, tmax = 0, 0
    for j in xrange(len(topodata)):
        plot = topodata[j] 
        
        bins = np.arange(min(plot.yp),max(plot.yp),1)
        inds = np.digitize(plot.yp,bins)
        distance = []
        moy_topo = []
        std_topo = []
        for j in range(len(bins)-1):
            uu = np.flatnonzero(inds == j)
            if len(uu)>0:
                distance.append(bins[j] + (bins[j+1] - bins[j])/2.)
                std_topo.append(np.std(-plot.z[uu]))
                moy_topo.append(np.mean(-plot.z[uu]))

        distance = np.array(distance)
        std_topo = np.array(std_topo)
        moy_topo = np.array(moy_topo)
        max_topo,min_topo = np.array(moy_topo)+std_topo, np.array(moy_topo)-std_topo

        #ax1.scatter(plot.yp,-plot.z,s=plot.width,marker='o',color=plot.color,label=plot.name)
        ax1.plot(distance,moy_topo,color='grey',lw=2.,label='topography')
        ax1.plot(distance,max_topo,color='black',lw=1.)
        ax1.plot(distance,min_topo,color='black',lw=1.)

        tmin, tmax = np.min(np.append(min_topo,tmin)), np.max(np.append(max_topo,tmax))
    
    plt.ylim([np.min(tmin)-0.5,np.max(tmax)+0.5])
    #plt.ylim([-1,6])
    #ax1.set_xlabel('Distance (km)')
    plt.setp( ax1.get_xticklabels(), visible = False)
    ax1.set_ylabel('Elevation (km)')
    ax1.yaxis.set_major_locator(tic.MaxNLocator(3))
    ax1.legend(loc = 'best',fontsize='x-small')
    ax1.set_title('Profile %s, azimuth: %s'%(profiles.name, strike))
    
    ax2 = fig.add_subplot(4,1,2)
    #ax2.axis('equal')
    ax2.set_xlim([-l/2,l/2])
    
    # plot decollement
    ax2.text(fmodel[0].fperp,-(fmodel[0].w+15),fmodel[0].name,color = 'black')
    ax2.text(fmodel[0].fperp,-(fmodel[0].w+20),'SS: %4.1f mm'%(fmodel[0].sst),style='italic',size='smaller') 
    ax2.text(fmodel[0].fperp,-(fmodel[0].w+25),'DS: %4.1f mm'%(fmodel[0].ds),style='italic',size='smaller') 
    # plot ramp and kink 
    for k in xrange(1,flt.structures[0].Mseg):
        ax2.text(fmodel[k].fperp+3,-((fmodel[k].w+fmodel[0].w)/1.3),fmodel[k].name,color = 'black')
        #ax2.text(fmodel[2].fperp+3,-((fmodel[2].w+fmodel[0].w)/1.3),fmodel[2].name,color = 'black')
        # plot ss main flower
        ax2.text(fmodel[k].fperp+3,-((fmodel[k].w+fmodel[0].w)/1.3+5),'SS: %4.1f mm'%((fmodel[k].ss)),style='italic',size='smaller') 
        #ax2.text(fmodel[2].fperp+3,-((fmodel[2].w+fmodel[0].w)/1.3+5),'SS: %4.1f mm'%(abs((fmodel[2].ss))),style='italic',size='smaller') 
        # plot ds main flower
        ax2.text(fmodel[k].fperp+3,-((fmodel[k].w+fmodel[0].w)/1.3+10),'DS: %4.1f mm'%((fmodel[k].ds)),style='italic',size='smaller') 
        #ax2.text(fmodel[2].fperp+3,-((fmodel[2].w+fmodel[0].w)/1.3+10),'DS: %4.1f mm'%((fmodel[2].ds)),style='italic',size='smaller') 
        # plot all plaussible models 
        for i in xrange(0,len(fmodel[0].tracew),10):
            ax2.plot([ fmodel[0].traceF[i], fmodel[k].traceF[i] ],[ -fmodel[0].tracew[i], -fmodel[k].tracew[i] ], '-' ,lw= .01, color='blue')
            #ax2.plot([ fmodel[0].traceF[i], fmodel[2].traceF[i] ],[ -fmodel[0].tracew[i], -fmodel[2].tracew[i] ], '-' ,lw= .01, color='blue')
   
    Mtemp = flt.structures[0].Mseg
    for j in xrange(1,Mstruc):
        for k in xrange(flt.structures[j].Mseg):
            # kink
            ax2.text(fmodel[Mtemp+k].fperp+3,-((fmodel[Mtemp+k].w+fmodel[Mtemp-1].w)/1.3),fmodel[Mtemp+k].name,color = 'black')
            ax2.text(fmodel[Mtemp+k].fperp+3,-((fmodel[Mtemp+k].w+fmodel[Mtemp-1].w)/1.3+5),'SS: %4.1f mm'%(fmodel[Mtemp+k].ss),style='italic',size='smaller') 
            ax2.text(fmodel[Mtemp+k].fperp+3,-((fmodel[Mtemp+k].w+fmodel[Mtemp-1].w)/1.3+10),'DS: %4.1f mm'%(fmodel[Mtemp+k].ds),style='italic',size='smaller') 
            for i in xrange(0,len(fmodel[0].tracew),10):
                ax2.plot([ fmodel[Mtemp-1].traceF[i], fmodel[Mtemp+k].traceF[i] ],[ -fmodel[Mtemp-1].tracew[i], -fmodel[Mtemp+k].tracew[i] ], '-' ,lw= .01, color='blue')
        Mtemp += flt.structures[j].Mseg

    for i in xrange(len(plotdata)):
        plot = plotdata[i]
        ax2.scatter(plot.yp,-plot.z,s = plot.width,marker = 'o',color = plot.color,label = plot.name)
    
    plt.setp( ax2.get_xticklabels(), visible = False)
    wmax = fmodel[0].w+10
    ax2.set_ylim([-wmax,6])
    ax2.set_ylabel('Elevation (km)')
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
    uslos = usx*proj[0]+usy*proj[1]+usz*proj[2]
    print proj[0],proj[1],proj[2]
    # default value if no insar
    sigmalos = np.ones(uslos.shape[0]) 

    for j in xrange(Mseg):
        ax2.plot([fmodel[j].fperp,fmodel[j].fperp],[5,-fmodel[j].w],'--',lw=1.,color = 'black')

    # 2) Plot 
    markers = ['+','d','x','v']
    colors = ['orange','m','yellow','red','blue']
    
    # los component
    ax3 = fig.add_subplot(4,1,3)
    ax3.set_xlim([-l/2,l/2])
    plt.setp( ax3.get_xticklabels(), visible = False)
    ax3.set_ylabel('LOS (mm/yr)')
    ax3.yaxis.set_major_locator(tic.MaxNLocator(3)) 
    # fault parrallele
    ax4 = fig.add_subplot(4,1,4)
    ax4.set_xlim([-l/2,l/2])
    ax4.plot(ysp,uspar,'-',lw = 2.)
    # fault perpendicular
    ax4.plot(ysp,usperp,'-',lw = 2.)
    ax4.plot(ysp,usz,'-',lw = 2.)
    ax4.set_ylabel('Velocity (mm/yr)')
    ax4.yaxis.set_major_locator(tic.MaxNLocator(5))
    ax4.set_xlabel('Distance (km)')
    
    # 3) InSAR data and model
    colors = ['blue','red','aqua','orange']
    ymin,ymax = ax3.get_ylim()
    for i in xrange(len(insardata)):
        insar = insardata[i]

        # Create Baseline component
        binsarlos = insar.a+insar.b*insar.yp
        ilos = insar.ulos-binsarlos

        ax3.scatter(insar.yp,insar.ulos-binsarlos,s = 3.,marker = 'o',label = insardata[i].reduction,color = insardata[i].color)
        # problem if several insar with differents weight...
        sigmalos *=  insar.sigmad[0]

        ymean = np.mean(ilos)
        tempmax = ymean + 2.5*np.std(ilos)
        tempmin = ymean - 2.5*np.std(ilos)
        ymax = np.max([ymax,tempmax])
        ymin = np.min([ymin,tempmin])

    ax3.set_ylim([ymin,ymax])
    
    ax3.plot(ysp,uslos,'-',color = 'red',label = 'model',lw = 2.)
    ax3.plot(ysp,uslos-sigmalos,'-',color = 'red',label = 'uncertainty',lw = .5)
    ax3.plot(ysp,uslos+sigmalos,'-',color = 'red',lw = .5)
  
    # 4) GPS data and model
    for i in xrange(len(gpsdata)):
        gps = gpsdata[i]
        markers = ['^','v','+','x']

        bpar = gps.a
        #print bpar
        bperp = gps.b
        #print bperp
        bv = gps.c
        
        # remove wieght in gps uncertainties
        #wd = gps.wd
        #gps.sigmapar, gps.sigmaperp = (gps.sigmapar/wd)*2 , (gps.sigmaperp/wd)*2 
        
        ax4.plot(gps.yp,gps.upar-bpar,markers[i],color = 'blue',mew = 1.,label = '%s fault-parallel velocities'%gpsdata[i].reduction )
        ax4.errorbar(gps.yp,gps.upar-bpar,yerr = gps.sigmapar,ecolor = 'blue',barsabove = 'True',fmt = None)

        ax4.plot(gps.yp,gps.uperp-bperp,markers[i],color = 'green',mew = 1.,label = '%s fault-perpendicular velocities'%gpsdata[i].reduction)
        ax4.errorbar(gps.yp,gps.uperp-bperp,yerr = gps.sigmaperp,ecolor = 'green',fmt = None)
        #ymax,ymin = np.max(np.hstack([gps.upar-bpar+gps.sigmapar,gps.uperp-bperp+gps.sigmaperp]))+5,np.min(np.hstack([gps.upar-bpar-gps.sigmapar,gps.uperp-bperp-gps.sigmaperp]))-5
        ymax,ymin = np.max(np.hstack([gps.upar-bpar,gps.uperp-bperp]))+4,np.min(np.hstack([gps.upar-bpar,gps.uperp-bperp]))-4
        # ax4.set_ylim([ymin,ymax])
        
        if 3 == gps.dim:
            #gps.sigmav = (gps.sigmav/wd)*2
            bx = bpar*s[0]+bperp*n[0] 
            by = bpar*s[1]+bperp*n[1]
            blos = bx*proj[0]+by*proj[1]+bv*proj[2]
            # add gps in los if dim is 3
            ax3.plot(gps.yp,gps.ulos-blos,'+',color='black',mew=1.,label='%s GPS LOS'%gpsdata[i].reduction)
            ax4.plot(gps.yp,gps.uv-bv,markers[i],color = 'red',mew = 1.,label = '%s vertical velocities'%gpsdata[i].reduction)
            ax4.errorbar(gps.yp,gps.uv-bv,yerr = gps.sigmav,ecolor = 'red',fmt = None)
            
            #ymax,ymin = np.max(np.hstack([gps.upar-bpar+gps.sigmapar,gps.uperp-bperp+gps.sigmaperp,gps.uv-bv+gps.sigmav]))+5,np.min(np.hstack([gps.upar-bpar-gps.sigmapar,gps.uperp-bperp-gps.sigmaperp,gps.uv-bv-gps.sigmav]))-5
            ymax,ymin = np.max(np.hstack([gps.upar-bpar,gps.uperp-bperp,gps.uv-bv]))+2,np.min(np.hstack([gps.upar-bpar,gps.uperp-bperp,gps.uv-bv]))-2
            # ax4.set_ylim([ymin,ymax])

    # Plot fault and legend    
    ax3.legend(bbox_to_anchor = (0.,1.02,1.,0.102),loc = 3,ncol = 2,mode = 'expand',borderaxespad = 0.,fontsize = 'x-small')
    #ax3.legend(loc = 'best')
    ymin,ymax = ax3.get_ylim()
    for j in xrange(Mseg):
        ax3.plot([fmodel[j].fperp,fmodel[j].fperp],[ymin,ymax],'--',lw=1,color = 'black')
    
    ax4.legend(bbox_to_anchor = (0.,1.02,1.,0.102),loc = 3,ncol = 2,mode = 'expand',borderaxespad = 0.,fontsize ='x-small')
    #ax4.legend(loc = 4,fontsize = 'x-small')
    ymin,ymax = ax4.get_ylim()
    for j in xrange(Mseg):
        ax4.plot([fmodel[j].fperp,fmodel[j].fperp],[ymin,ymax],'--',lw=1,color = 'black')
    
    fig.savefig(outdir+'/profiles/'+name+'_LOS.eps', format = 'EPS')
    
def plotMap(flt,x1,y1,nfigure):
    
    profiles = flt.profiles
    fmodel = flt.fmodel
    Mseg = flt.Mseg
    insardata = flt.insardata
    gpsdata = flt.gpsdata
    outdir = flt.outdir
    gmtfiles = flt.gmtfiles

    fig = plt.figure(nfigure,figsize = (18,8))
  
    fig.subplots_adjust()
    ax1 = fig.add_subplot(1,3,1) #(row,column,number)
    ax2 = fig.add_subplot(1,3,2) #(row,column,number)
    ax3 = fig.add_subplot(1,3,3) #(row,column,number)
    markers = ['^','v','/^']
    colors = ['orange','m','yellow','red','blue']
    xmin,xmax = profiles.x-100,profiles.x+100
    # bundary profile
    xp,yp = np.zeros((7)),np.zeros((7))
    x0 = profiles.x
    y0 = profiles.y
    l = profiles.l
    w = profiles.w
    name = profiles.name
    xp0 = profiles.xp0
    yp0 = profiles.yp0
    xpmin,xpmax = profiles.xpmin,profiles.xpmax
    ypmin,ypmax = profiles.ypmin,profiles.ypmax
    
    strike = flt.strike
    proj = profiles.proj
    outdir = outdir
    #var_insar = profiles.var_insar
    #var_gps = profiles.var_gps
    #var = profiles.var
    s = flt.s
    n = flt.n

    xp[:] = x0-w/2*s[0]-l/2*n[0],x0+w/2*s[0]-l/2*n[0],x0+w/2*s[0]+l/2*n[0],x0-w/2*s[0]+l/2*n[0],x0-w/2*s[0]-l/2*n[0],x0-l/2*n[0],x0+l/2*n[0]
    yp[:] = y0-w/2*s[1]-l/2*n[1],y0+w/2*s[1]-l/2*n[1],y0+w/2*s[1]+l/2*n[1],y0-w/2*s[1]+l/2*n[1],y0-w/2*s[1]-l/2*n[1],y0-l/2*n[1],y0+l/2*n[1]
 
    fid = open(outdir+'/profiles/'+name+'_coord.xy','w')
    np.savetxt(fid, np.vstack([xp[:],yp[:]]).T ,header = 'x(km)     y(km) ',comments = '# ')
    fid.write('\n')
    fid.close

    #fault model
    for j in xrange(1,Mseg):
        fmodel[j].x = math.cos(flt.str)*fmodel[j].fperp+fmodel[0].x
        fmodel[j].y = -math.sin(flt.str)*fmodel[j].fperp+fmodel[0].y

    xf,yf = np.zeros((Mseg,2)),np.zeros((Mseg,2))
    for j in xrange(Mseg):
        xf[j,0] = fmodel[j].x+2*-100*s[0]
        xf[j,1] = fmodel[j].x+2*100*s[0]
        yf[j,0] = fmodel[j].y+2*-100*s[1]
        yf[j,1] = fmodel[j].y+2*100*s[1]
    
    ##InSAR DATA
    for i in xrange(len(insardata)):
        insar = insardata[i]
        imax = np.percentile(insar.ulos,99.5)
        kk = np.flatnonzero(insar.ulos>imax)
        insar.ulos[kk] = 0
        
        name = insar.reduction
        orb = insar.a+insar.b*insar.yp
        vmax = np.mean((insar.ulos-orb)) + np.percentile(abs(insar.ulos-orb),98)
        vmin = np.mean((insar.ulos-orb)) - np.percentile(abs(insar.ulos-orb),98)

        norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
        m = cm.ScalarMappable(norm = norm, cmap = cm.jet)
        
        index = np.nonzero((insar.xp>xpmax)|(insar.xp<xpmin)|(insar.yp>ypmax)|(insar.yp<ypmin))
        insarx,insary,los,model = np.delete(insar.x,index),np.delete(insar.y,index),np.delete(insar.ulos-orb,index),np.delete(insar.mlos,index)
       
        fid = open(outdir+'/insar/'+name+'.xylos','w')
        np.savetxt(fid, np.vstack([insarx, insary, los]).T ,header = 'x(km)     y(km)    los(mm/yr)  ',comments = '# ')
        fid.write('\n')
        fid.close

        fid = open(outdir+'/insar/'+name+'model.xylos','w')
        np.savetxt(fid, np.vstack([insarx, insary, model]).T ,header = 'x(km)     y(km)    los(mm/yr)  ',comments = '# ')
        fid.write('\n')
        fid.close

        fid=open(outdir+'/insar/'+name+'residus.xylos','w')
        np.savetxt(fid, np.vstack([insarx, insary, los-model]).T ,header='x(km)     y(km)    los(mm/yr)  ',comments='# ')
        fid.write('\n')
        fid.close

        m.set_array(los)
        
        # InSAR model 
        facelos = m.to_rgba(los)
        facemodel = m.to_rgba(model)
        faceres = m.to_rgba(los-model) 

        ax1.scatter(insarx,insary,s = 30,marker = 'o',color = facelos,label = 'LOS Velocity %s'%(insar.reduction))
        ax2.scatter(insarx,insary,s = 30,marker = 'o',color = facemodel,label = 'LOS Velocity %s'%(insar.reduction))
        ax3.scatter(insarx,insary,s = 30,marker = 'o',color = faceres,label = 'LOS Velocity %s'%(insar.reduction))
        
        ax1.arrow(x1,y1,30*proj[0],30*proj[1],fill = True,width = 0.5,color = 'black')
        ax2.arrow(x1,y1,30*proj[0],30*proj[1],fill = True,width = 0.5,color = 'black')
        ax3.arrow(x1,y1,30*proj[0],30*proj[1],fill = True,width = 0.5,color = 'black')
     
    if len(insardata)>0:
        fig.colorbar(m,shrink = 0.5, aspect = 5)

    for i in xrange(len(gpsdata)):
        gps = gpsdata[i]
        name = gps.reduction

        index = np.nonzero((gps.xp>xpmax)|(gps.xp<xpmin)|(gps.yp>ypmax)|(gps.yp<ypmin))
        #index = np.nonzero((gps.xp>60)|(gps.xp<-60)|(gps.yp>200)|(gps.yp<-200))
        
        if gps.dim==2:
            gpsname,gpsx,gpsy,gpsux,gpsuy,gpsmx,gpsmy,gpssigmax,gpssigmay = np.delete(gps.name,index),np.delete(gps.x,index),np.delete(gps.y,index),np.delete(gps.ux,index),np.delete(gps.uy,index),np.delete(gps.mx,index),np.delete(gps.my,index),np.delete(gps.sigmax,index),np.delete(gps.sigmay,index)
        else:
            gpsname,gpsx,gpsy,gpsux,gpsuy,gpsuv,gpsmx,gpsmy,gpsmz,gpssigmax,gpssigmay,gpssigmav = np.delete(gps.name,index),np.delete(gps.x,index),np.delete(gps.y,index),np.delete(gps.ux,index),np.delete(gps.uy,index),np.delete(gps.uv,index),np.delete(gps.mx,index),np.delete(gps.my,index),np.delete(gps.mz,index),np.delete(gps.sigmax,index),np.delete(gps.sigmay,index),np.delete(gps.sigmav,index)
        # remove wieght in gps uncertainties
        wd = gps.wd
        gpssigmax,gpssigmay = gpssigmax/wd , gpssigmay/wd
        
        bx = (gps.a*flt.s[0]+gps.b*flt.n[0])
        by = (gps.a*flt.s[1]+gps.b*flt.n[1])
        #bx,by = np.mean(gpsux-gpsmx),np.mean(gpsuy-gpsmy)

        
        if gps.dim==2:
        
            fid = open(outdir+'/gps/'+name+'_%s.psvelo'%(i),'w')
            np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux,gpsuy,gpssigmax,gpssigmay]).T ,header = 'x(km)  y(km)   East_Vel_(mm/yr)    North_Vel_(mm/yr)   East_StdDev_(mm/yr) North_StdDev_(mm/yr)',comments = '# ')
            fid.write('\n')
            fid.close

            fid = open(outdir+'/gps/'+name+'model_%s.psvelo'%(i),'w')
            np.savetxt(fid, np.vstack([gpsx, gpsy, gpsmx,gpsmy,gpssigmax,gpssigmay]).T ,header = 'x(km)  y(km)   East_Vel_(mm/yr)    North_Vel_(mm/yr)   East_StdDev_(mm/yr) North_StdDev_(mm/yr)',comments = '# ')
            fid.write('\n')
            fid.close    
            
            fid=open(outdir+'/gps/'+name+'residus_%s.psvelo'%(i),'w')
            np.savetxt(fid, np.vstack([gpsx, gpsy,gpsux-gpsmx-bx,gpsuy-gpsmy-by,gpssigmax,gpssigmay]).T ,header='x(km)  y(km)   East_Dev_(mm/yr)    North_Dev_(mm/yr)   East_StdDev_(mm/yr) North_StdDev_(mm/yr)',comments='# ')
            fid.write('\n')
            fid.close
        
            fid = open(outdir+'/gps/'+name+'_%s_rf.psvelo'%(i),'w')
            np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux-bx,gpsuy-by,gpssigmax,gpssigmay]).T ,header = 'x(km)   y(km)   East_Vel_(mm/yr)    North_Vel_(mm/yr)   East_StdDev_(mm/yr) North_StdDev_(mm/yr)',comments = '# ')
            fid.write('\n')
            fid.close
        
        else:
            
            bv=gps.c
            
            fid = open(outdir+'/gps/'+name+'_%s.psvelo'%(i),'w')
            np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux,gpsuy,gpsuv,gpssigmax,gpssigmay,gpssigmav]).T ,header = 'x(km)  y(km)   East_Vel_(mm/yr)    North_Vel_(mm/yr) Up_Vel_(mm/yr)  East_StdDev_(mm/yr) North_StdDev_(mm/yr) Up_StdDev_(mm/yr)',comments = '# ')
            fid.write('\n')
            fid.close

            fid = open(outdir+'/gps/'+name+'model_%s.psvelo'%(i),'w')
            np.savetxt(fid, np.vstack([gpsx, gpsy, gpsmx,gpsmy,gpsmz,gpssigmax,gpssigmay,gpssigmav]).T ,header = 'x(km)  y(km)   East_Vel_(mm/yr)    North_Vel_(mm/yr) Up_Vel_(mm/yr)  East_StdDev_(mm/yr) North_StdDev_(mm/yr) Up_StdDev_(mm/yr)',comments = '# ')
            fid.write('\n')
            fid.close    
            
            fid = open(outdir+'/gps/'+name+'_%s_rf.psvelo'%(i),'w')
            np.savetxt(fid, np.vstack([gpsx, gpsy, gpsux-bx,gpsuy-by,gpsuv-bv,gpssigmax,gpssigmay,gpssigmav]).T ,header = 'x(km)   y(km)   East_Vel_(mm/yr)    North_Vel_(mm/yr)  Up_Vel_(mm/yr)  East_StdDev_(mm/yr) North_StdDev_(mm/yr) Up_StdDev_(mm/yr)',comments = '# ')
            fid.write('\n')
            fid.close
        
            fid=open(outdir+'/gps/'+name+'residus_%s.psvelo'%(i),'w')
            np.savetxt(fid, np.vstack([gpsx, gpsy,gpsux-gpsmx-bx,gpsuy-gpsmy-by,gpsuv-gpsmz-bv,gpssigmax,gpssigmay,gpssigmav]).T ,header='x(km)  y(km)   East_Dev_(mm/yr)    North_Dev_(mm/yr)  Up_Dev_(mm/yr)   East_StdDev_(mm/yr) North_StdDev_(mm/yr)  Up_StdDev_(mm/yr)',comments='# ')
            fid.write('\n')
            fid.close

        
        #data
        data = ax1.quiver(gpsx,gpsy,gpsux-bx,gpsuy-by,scale = 100,width = 0.005,color = 'black') #saller scale for arrow longer
        #data = ax1.quiver(gpsx,gpsy,gpsux,gpsuy,scale = 100,width = 0.005,color = 'black') #saller scale for arrow longer
        #model
        #model = ax1.quiver(gpsx,gpsy,gpsmx,gpsmy,scale = 100,width = 0.005,color = 'red')
        model = ax2.quiver(gpsx,gpsy,gpsmx,gpsmy,scale = 100,width = 0.005,color = 'red')
        #model = ax2.quiver(gpsx,gpsy,gpsmx+bx,gpsmy+by,scale = 100,width = 0.005,color = 'red')
        #residus
        residus = plt.quiver(gpsx,gpsy,gpsux-gpsmx-bx,gpsuy-gpsmy-by,scale = 100,width = 0.005,color = 'purple')
                
        ax1.scatter(gpsx, gpsy, c = colors[i], s = 30, marker = markers[i], label = gps.reduction)
        ax2.scatter(gpsx, gpsy, c = colors[i], s = 30, marker = markers[i], label = gps.reduction)
        ax3.scatter(gpsx, gpsy, c = colors[i], s = 30, marker = markers[i], label = gps.reduction)
        ## display name of the station
        #for kk in xrange(len(gpsname)):
                #ax1.text(gpsx[kk]-4*kk,gpsy[kk]-10,gpsname[kk],color = 'black')
        ax1.quiverkey(data,0.1,1.015,20.,'GPS velocities',coordinates = 'axes',color = 'black')
        ax2.quiverkey(model,0.1,1.015,20.,'Model',coordinates = 'axes',color = 'red')
        ax3.quiverkey(model,1.3,1.015,20,'Residuals',coordinates = 'axes',color = 'red')

    ####DATA
    # plot gtm files
    for ii in xrange(len(gmtfiles)):
                name = gmtfiles[ii].name
                wdir = gmtfiles[ii].wdir
                filename = gmtfiles[ii].filename
                color = gmtfiles[ii].color
                width = gmtfiles[ii].width

                fx,fy = gmtfiles[ii].load()
                #ax1.plot(fx[0],fy[0],color = color,lw = width,label = name) 
                for i in xrange(len(fx)):
                    ax1.plot(fx[i],fy[i],color = color,lw = width)

    ax1.plot(xp[:],yp[:],color = 'black',lw = 2.)
    #ax1.text(profiles.x+profiles.l/2*n[0],profiles.y+profiles.l/2*n[1],profiles.name) 
    #plot fault model
    for f in xrange(Mseg):
        #ax1.plot(xf[f,:],yf[f,:],'--',color = 'black',lw = 4.,label = '%s fault'%(fmodel[f].name))
        ax1.plot(xf[f,:],yf[f,:],'--',color = 'black',lw = 1.)
    ax1.legend(loc='best',fontsize='x-small')
                
    ax1.text(x1,y1+5,'LOS',color = 'black')
    ax1.set_xlabel('Distance (km)')
    ax1.set_title('Data')
            
    ###MODEL
    # plot gtm files
    for ii in xrange(len(gmtfiles)):
                name = gmtfiles[ii].name
                wdir = gmtfiles[ii].wdir
                filename = gmtfiles[ii].filename
                color = gmtfiles[ii].color
                width = gmtfiles[ii].width

                fx,fy = gmtfiles[ii].load()
                #ax2.plot(fx[0],fy[0],color = color,lw = width,label = name) 
                for i in xrange(len(fx)): 
                    ax2.plot(fx[i],fy[i],color = color,lw = width)
         
    # plot profile
    ax2.plot(xp[:],yp[:],color = 'black',lw = 1.)
    #ax2.text(profiles.x+profiles.l/2*n[0],profiles.y+profiles.l/2*n[1],profiles.name) 
    #plot fault model
    for f in xrange(Mseg):
        ax2.plot(xf[f,:],yf[f,:],'--',color = 'black',lw = 1.)
        #ax2.plot(xf[f,:],yf[f,:],'--',color = 'black',lw = 4.,label = '%s fault'%(fmodel[f].name))
    ax2.legend(loc='best', fontsize='x-small')
    ax2.text(x1,y1+5,'LOS',color = 'black')
    ax2.set_xlabel('Distance(km)')
    ax2.set_title('Model')
            
    ###RESIDUS
    # plot gtm files
    for ii in xrange(len(gmtfiles)):
                name = gmtfiles[ii].name
                wdir = gmtfiles[ii].wdir
                filename = gmtfiles[ii].filename
                color = gmtfiles[ii].color
                width = gmtfiles[ii].width

                fx,fy = gmtfiles[ii].load()
                for i in xrange(len(fx)): 
                    ax3.plot(fx[i],fy[i],color = color,lw = width)
        
    # plot profile
    ax3.plot(xp[:],yp[:],color = 'black',lw = 1.)
    #ax3.text(profiles.x+profiles.l/2*n[0],profiles.y+profiles.l/2*n[1],profiles.name) 
    #plot fault model
    for f in xrange(Mseg):
        ax3.plot(xf[f,:],yf[f,:],'--',color = 'black',lw = 1.)
        #ax3.plot(xf[f,:],yf[f,:],'--',color = 'black',lw = 4.,label = '%s fault'%(fmodel[f].name))
    ax3.legend(loc='best', fontsize='x-small')
    ax3.text(x1,y1+5,'LOS',color = 'black')
    ax3.set_xlabel('Distance (km)')
    ax3.set_title('Residus')
    
    ax1.axis('equal')
    ax1.set_xlim(xmin,xmax)
    ax2.axis('equal')
    ax2.set_xlim(xmin,xmax)
    ax3.axis('equal')
    ax3.set_xlim(xmin,xmax)

    plt.suptitle('InSAR foward model')
    
    #fig.savefig(outdir+'/map/'+'basemap.ps',format = 'PS')
    fig.savefig(outdir+'/map/'+profiles.name+'_map.eps',format = 'EPS')

def plotHist(flt,model,nfigure):
   
    from scipy.stats import norm

    Mseg = flt.Mseg
    fmodel = flt.fmodel
    Mker = fmodel[0].Mker # Moche
    Sampled = flt.Sampled
    mu, tau = flt.mu, flt.tau
    mmin, mmax = flt.mmin, flt.mmax
    profile = flt.profiles
    outdir = flt.outdir

    # for histo
    traces = []
    labels = []
   
    fig = plt.figure(nfigure, figsize = (16,11))
    fig.subplots_adjust(hspace = 0.75,wspace = 0.35)
   
    def histo(traces,labels,Mseg,Mker,l,mu,tau,mmin,mmax):
        ll = 0
        for trace,label,mu,tau,mmin,mmax in zip(traces,labels,mu,tau,mmin,mmax):    
            index = l+ll*Mseg
            ax = fig.add_subplot(Mker,Mseg,index)
           
            # plot Posterior model
            histo = plt.hist(trace,bins = 20,normed = True, histtype = 'stepfilled', color = 'black', label = label)
            mean = trace.mean()
            sig = 2*np.std(trace)
            plt.axvline(x = mean, c = "red")
            plt.axvline(x = mean+sig, c = "green", linestyle = '--')
            plt.axvline(x = mean-sig, c = "green", linestyle = '--')
            
            # plot prior model
            xmin,xmax = ax.set_xlim(mmin,mmax)
            #ymin,ymax = ax.get_ylim()
            #x = np.arange(xmin,xmax,(xmax-xmin)/100)
            #plt.plot(x, ymax*norm.pdf(x,mu,tau))
            #plt.plot((mmin,mmax),(1,1))
            

            #plt.legend(fontsize = 'xx-small')
            plt.title(label)

            #if (l < 3 and ll == 1) or (l < 2): 
                #print mean,sig
                #plt.xlabel("Mean: {:0.3f} (mm/yr) +- \n{:0.3f} (mm/yr)".format(mean, sig))

            #else:
                #print mean,sig
                #plt.xlabel("Mean: {:0.1f} (km) +- \n{:0.1f} (km)".format(mean, sig))

            ll +=  1
    
    count = 0
    if '{} Strike Slip'.format(fmodel[0].name) in Sampled:
        strike = model.trace('{} Strike Slip'.format(fmodel[0].name))[:]
        traces.append(strike)
        labels.append('{} Strike Slip'.format(fmodel[0].name))
        count +=  1 
    if '{} Shortening'.format(fmodel[0].name) in Sampled:
        short = model.trace('{} Shortening'.format(fmodel[0].name))[:]
        traces.append(short)
        labels.append('{} Shortening'.format(fmodel[0].name))
        count +=  1 
    if '{} H'.format(fmodel[0].name) in Sampled:
        H1 = model.trace('{} H'.format(fmodel[0].name))[:]
        traces.append(H1)
        labels.append('{} H'.format(fmodel[0].name))
        count +=  1

    j = 1
    # plot histos main fault
    histo(traces,labels,Mseg,Mker,j,mu[:count],tau[:count],mmin[:count],mmax[:count])
   
    for j in xrange(1,Mseg):
        traces = []
        labels = []
        pcount = count
        if '{} Strike Slip'.format(fmodel[j].name) in Sampled:
            m = model.trace('{} Strike Slip'.format(fmodel[j].name))[:]
            traces.append(m)
            labels.append('{} Strike Slip'.format(fmodel[j].name))
            count +=  1 
            
        if '{} D'.format(fmodel[j].name) in Sampled:
            m = model.trace('{} D'.format(fmodel[j].name))[:]
            traces.append(m)
            labels.append('{} D'.format(fmodel[j].name))
            count +=  1 

        if '{} H'.format(fmodel[j].name) in Sampled:
            m = model.trace('{} H'.format(fmodel[j].name))[:]
            traces.append(m)
            labels.append('{} H'.format(fmodel[j].name))
            count +=  1 
        
        histo(traces,labels,Mseg,Mker,j+1,mu[pcount:count],tau[pcount:count],mmin[pcount:count],mmax[pcount:count])
    
    # save
    fig.savefig(outdir+'/stat/'+profile.name+'_histo.eps', format = 'EPS')
    

def plotDist(flt,model,nfigure):

    from scipy.stats import norm
    import plot2Ddist
    scatterstyle={'color':'r', 'alpha':0.1}
    styleargs = {'color':'k', 'scatterstyle':scatterstyle}

    fmodel=flt.fmodel
    profile=flt.profiles
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

