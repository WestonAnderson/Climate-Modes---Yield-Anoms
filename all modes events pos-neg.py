#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 10:02:28 2018
@author: weston

NOTES: PCs were coded to lag early spring backward. So loading in A-M-M-J-J-A-S-O-N-D 
        in PC-year 1998 will be calendar year 1998, but in J-F-M in 
        PC-year 1998 will be calendar year 1999

EDIT - 1/29/18 - change the data source to be the post-processed data in the data folder
                instead of the data in the project files
               - Include precip, geopotential height, winds and anomalous ascent
"""

import numpy as np
import netCDF4
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #
years = [1981,2011] #from the SVD analysis
lvl = 925
saveNotes = ''#'_ensoLC'
modes=['IOD']
thresh = .65
#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#ENSO life cycle years
#'EN': [1982,1987,1997,2009,2006,1994]#2006,1994, 'LN':[1983,1988,1998,2010,2007,1995]#2007,1995
#EN':[1997-1,1982-1,2009-1,1987-1],'LN':

#strong years
#'EN':[1982,1997,2009,2002,1991,1987], 'LN':[1988,1998,2007,1983,1995,2010]


posModeYrs = {'ENSO':[1982-1,1997-1,2009-1,2002-1],'NAO':[1998-1,1994-1,1983-1,1988-1],
              'IOD':[1992-1,1996-1,1981-1,1984-1],'TAV (AMM)':[2010-1,2005-1,1981-1,1997-1]}

negModeYrs = {'ENSO':[1988-1,1998-1,2007-1,1995-1],'NAO':[2009-1,1997-1,1984-1,1986-1],
              'IOD':[1994-1,1997-1,1982-1,1987-1],'TAV (AMM)':[1991-1,1986-1,1994-1,2002-1]}

#Read in agreement mask. agreeMSK(DF, frac, require agree w/ avg? (Y/N) )
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/agreement.py').read())
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/custom_div_cmap.py').read())
prcp_map = custom_div_cmap(20, mincol='saddlebrown', midcol='white' ,maxcol='CornflowerBlue')

yrList = np.array(range(years[0],years[1]+1))
seasons = ['DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ']
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = 'Jan-DecENSOyr0har'


#read in the netCDF file containing the variables
ncGeo = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/gph.nc',
                        'r',format='NETCDF4') #read in the netCDF file for the variable
geoVar=ncGeo.variables['hgt'][(ncGeo.variables['year'][:]>=years[0])&(ncGeo.variables['year'][:]<years[1]+2),
                              :,np.where(ncGeo.variables['level'][:]==lvl)[0][0],:,:]

ncSM = netCDF4.Dataset('/Volumes/Data_Archive/Data/SoilMoisture/GLDAS/Noah10_nc/detrended/SoilM.nc','r',format='NETCDF4') 
sm = np.nansum(ncSM.variables['SoilM'][(ncSM.variables['year'][:]>=years[0])&(ncSM.variables['year'][:]<years[1]+2),:,:3,:,:],2)

ncP = netCDF4.Dataset('/Volumes/Data_Archive/Data/Precip/GPCP/detrended/GPCP_precip.nc','r',format='NETCDF4')
prcp = ncP.variables['precip'][(ncP.variables['year'][:]>=years[0])&(ncP.variables['year'][:]<years[1]+2),:,:,:]

ncMaxT = netCDF4.Dataset('//Volumes/Data_Archive/Data/Temp/BerkeleyEarth/monthly/Complete_TMAX_LatLong1.nc','r',format='NETCDF4') #read in the netCDF file for the variable
maxT = ncMaxT.variables['temperature'][(ncMaxT.variables['time'][:]>=years[0])&(ncMaxT.variables['time'][:]<years[1]+2),:,:] 

ncUwnd = netCDF4.Dataset('/Volumes/Data_Archive/Data/Wind/detrend/uwnd.nc','r',format='NETCDF4') #read in the netCDF file for the variable
Uwnd=ncUwnd.variables['uwnd'][(ncUwnd.variables['year'][:]>=years[0])&(ncUwnd.variables['year'][:]<years[1]+2),
                              :,np.where(ncUwnd.variables['level'][:]==lvl)[0][0],:,:]
ncVwnd = netCDF4.Dataset('/Volumes/Data_Archive/Data/Wind/detrend/vwnd.nc','r',format='NETCDF4') #read in the netCDF file for the variable
Vwnd=ncVwnd.variables['vwnd'][(ncVwnd.variables['year'][:]>=years[0])&(ncVwnd.variables['year'][:]<years[1]+2),
                              :,np.where(ncVwnd.variables['level'][:]==lvl)[0][0],:,:]
ncOmega = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/omega.nc','r',format='NETCDF4')
omega= ncOmega.variables['omega'][(ncOmega.variables['year'][:]>=years[0])&(ncOmega.variables['year'][:]<years[1]+2),
                              :,np.where(ncOmega.variables['level'][:]==lvl)[0][0],:,:]
#ncP = netCDF4.Dataset('/Volumes/Data_Archive/Data/Precip/GPCC_nc/detrended/sumP.nc','r',format='NETCDF4') #read in the netCDF file for the variable



#Read in the SST anomaly data
sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
#T is months since 1960-01-01, and I want to start at 1950 with one year lag
sstLats=sstData.variables['Y'][:]
sstLons=sstData.variables['X'][:]
sstLons,sstLats=np.meshgrid(sstLons,sstLats)

sstAnom=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-years[0]))&
                        (sstData.variables['T'][:]<=12*(2+years[1]-1960)),...]
sstAnom = np.squeeze(sstAnom);sstAnom=np.ma.filled(sstAnom,0)
sstAnom_M = np.squeeze(sstAnom);sstAnom=np.ma.filled(sstAnom,0)


#Remove the climatology to get anomalies by month
for m in range(0,12):
    geoVar[:,m,...] = geoVar[:,m,...] - np.nanmean(geoVar[:,m,...],axis=0)
    omega[:,m,...] = omega[:,m,...] - np.nanmean(omega[:,m,...],axis=0)
    prcp[:,m,...] = prcp[:,m,...] - np.nanmean(prcp[:,m,...],axis=0)
    sm[:,m,...] = sm[:,m,...] - np.nanmean(sm[:,m,...],axis=0)
    Uwnd[:,m,...] = Uwnd[:,m,...] - np.nanmean(Uwnd[:,m,...],axis=0)
    Vwnd[:,m,...] = Vwnd[:,m,...] - np.nanmean(Vwnd[:,m,...],axis=0)
maxT[:,...] = maxT[:,...] - np.nanmean(maxT[:,...],axis=0)

#Rearrange to have only one time dimension
geoVar_M = geoVar.reshape([geoVar.shape[0]*geoVar.shape[1],geoVar.shape[2],geoVar.shape[3]])
omega_M = omega.reshape([omega.shape[0]*omega.shape[1],omega.shape[2],omega.shape[3]])
prcp_M = prcp.reshape([prcp.shape[0]*prcp.shape[1],prcp.shape[2],prcp.shape[3]])
sm_M = sm.reshape([sm.shape[0]*sm.shape[1],sm.shape[2],sm.shape[3]])
maxT_M = maxT#.reshape([maxT.shape[0]*maxT.shape[1],maxT.shape[2],maxT.shape[3]])
Uwnd_M = Uwnd.reshape([Uwnd.shape[0]*Uwnd.shape[1],Uwnd.shape[2],Uwnd.shape[3]])
Vwnd_M = Vwnd.reshape([Vwnd.shape[0]*Vwnd.shape[1],Vwnd.shape[2],Vwnd.shape[3]])
          
NCyrs =ncGeo.variables['year'][(ncGeo.variables['year'][:]>=years[0])&(ncGeo.variables['year'][:]<years[1]+2)]
gphLats=ncGeo.variables['latitude'][:]
gphLons=ncGeo.variables['longitude'][:]#-180.
gphLons, gphLats = np.meshgrid(gphLons,gphLats)
pLats=ncP.variables['latitude'][:]
pLons=ncP.variables['longitude'][:]#-180.
pLons, pLats = np.meshgrid(pLons,pLats)

mLats=ncMaxT.variables['latitude'][:]
mLons=ncMaxT.variables['longitude'][:]#-180.
mLons, mLats = np.meshgrid(mLons,mLats)
smLats=ncSM.variables['latitude'][:]
smLons=ncSM.variables['longitude'][:]#-180.
smLons, smLats = np.meshgrid(smLons,smLats)


for iS in [7]:
    for mode in modes:
        posCompYrs = posModeYrs[mode]
        negCompYrs = negModeYrs[mode]
        posYrInd = np.array([iY == NCyrs for iY in posCompYrs])
        posYrInd = np.where(posYrInd.sum(0)==1)[0]
        
        negYrInd = np.array([iY == NCyrs for iY in negCompYrs])
        negYrInd = np.where(negYrInd.sum(0)==1)[0]
        yrOffset = 0 + (iS//12)
        posYrInd = np.tile(np.arange((iS%12-1),(iS%12+2)),posYrInd.shape[0])+ \
                  np.repeat(posYrInd+yrOffset+1,3)*12
        negYrInd = np.tile(np.arange((iS%12-1),(iS%12+2)),negYrInd.shape[0])+ \
                  np.repeat(negYrInd+yrOffset+1,3)*12   
        limPos = np.int(posYrInd.shape[0]-np.ceil(np.sum(posYrInd>=sm_M.shape[0])/3)*3)
        limNeg = np.int(negYrInd.shape[0]-np.ceil(np.sum(negYrInd>=sm_M.shape[0])/3)*3)
                  
        #index
        geoVar = geoVar_M[posYrInd,...]
        posGPH = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
        geoVar = geoVar_M[negYrInd,...]
        negGPH = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    
        omega = omega_M[posYrInd,...]
        posOm = np.nanmean(omega.reshape([omega.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
        omega = omega_M[negYrInd,...]
        negOm = np.nanmean(omega.reshape([omega.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    
        Uwnd = Uwnd_M[posYrInd,...]
        posUwnd = np.nanmean(Uwnd.reshape([Uwnd.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
        Uwnd = Uwnd_M[negYrInd,...]
        negUwnd = np.nanmean(Uwnd.reshape([Uwnd.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    
        Vwnd = Vwnd_M[posYrInd,...]
        posVwnd = np.nanmean(Vwnd.reshape([Vwnd.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
        Vwnd = Vwnd_M[negYrInd,...]
        negVwnd = np.nanmean(Vwnd.reshape([Vwnd.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    
        prcp = prcp_M[posYrInd,...]
        oceans = mpl_toolkits.basemap.maskoceans(pLons-180,pLats,prcp[0,...])
        lndMsk = oceans.mask*1
        lndMsk,lons = mpl_toolkits.basemap.shiftgrid(0,lndMsk,pLons[0,:]-180)
        prcp = prcp*lndMsk
        prcp[prcp==0]=np.nan
        posPrcp = np.nanmean(prcp.reshape([prcp.shape[0]//3,3,pLons.shape[0],pLons.shape[1]]),1)
        prcp = prcp_M[negYrInd,...]
        negPrcp = np.nanmean(prcp.reshape([prcp.shape[0]//3,3,pLons.shape[0],pLons.shape[1]]),1)
    
        sm = sm_M[posYrInd[:limPos],...]
        posSM = np.nanmean(sm.reshape([sm.shape[0]//3,3,smLons.shape[0],smLons.shape[1]]),1)
        sm = sm_M[negYrInd[:limNeg],...]
        negSM = np.nanmean(sm.reshape([sm.shape[0]//3,3,smLons.shape[0],smLons.shape[1]]),1)
    
        maxT = maxT_M[posYrInd,...]
        posMaxT = np.nanmean(maxT.reshape([maxT.shape[0]//3,3,mLons.shape[0],mLons.shape[1]]),1)
        maxT = maxT_M[negYrInd,...]
        negMaxT = np.nanmean(maxT.reshape([maxT.shape[0]//3,3,mLons.shape[0],mLons.shape[1]]),1)
    
        sstAnom = sstAnom_M[posYrInd,...]
        posSST = np.nanmean(sstAnom.reshape([sstAnom.shape[0]//3,3,sstLons.shape[0],sstLons.shape[1]]),1)
        sstAnom = sstAnom_M[negYrInd,...]
        negSST = np.nanmean(sstAnom.reshape([sstAnom.shape[0]//3,3,sstLons.shape[0],sstLons.shape[1]]),1)
                  
        posGPH_all = np.nanmean(posGPH,0);posGPH[posGPH==0]=np.nan
        negGPH_all = np.nanmean(negGPH,0);negGPH[negGPH==0]=np.nan
        posGPH = agreeMSK(posGPH,thresh,'N')*np.nanmean(posGPH,0);posGPH[posGPH==0]=np.nan
        negGPH = agreeMSK(negGPH,thresh,'N')*np.nanmean(negGPH,0);negGPH[negGPH==0]=np.nan   

        posOm_all = np.nanmean(posOm,0);
        negOm_all = np.nanmean(negOm,0);
        posOm = agreeMSK(posOm,thresh,'N')*np.nanmean(posOm,0);posOm[posOm==0]=np.nan
        negOm = agreeMSK(negOm,thresh,'N')*np.nanmean(negOm,0);negOm[negOm==0]=np.nan

        posUwnd_all = np.nanmean(posUwnd,0);
        negUwnd_all = np.nanmean(negUwnd,0);
        posUwnd = agreeMSK(posUwnd,thresh,'N')*np.nanmean(posUwnd,0);posUwnd[posUwnd==0]=np.nan
        negUwnd = agreeMSK(negUwnd,thresh,'N')*np.nanmean(negUwnd,0);negUwnd[negUwnd==0]=np.nan

        posVwnd_all = np.nanmean(posVwnd,0);
        negVwnd_all = np.nanmean(negVwnd,0);
        posVwnd = agreeMSK(posVwnd,thresh,'N')*np.nanmean(posVwnd,0);posVwnd[posVwnd==0]=np.nan
        negVwnd = agreeMSK(negVwnd,thresh,'N')*np.nanmean(negVwnd,0);negVwnd[negVwnd==0]=np.nan

        posPrcp = agreeMSK(posPrcp,thresh,'N')*np.nanmean(posPrcp,0);posPrcp[posPrcp==0]=np.nan
        negPrcp = agreeMSK(negPrcp,thresh,'N')*np.nanmean(negPrcp,0);negPrcp[negPrcp==0]=np.nan
        posSST = agreeMSK(posSST,thresh,'N')*np.nanmean(posSST,0);posSST[posSST==0]=np.nan
        negSST = agreeMSK(negSST,thresh,'N')*np.nanmean(negSST,0);negSST[negSST==0]=np.nan        
        posMaxT = agreeMSK(posMaxT,thresh,'N')*np.nanmean(posMaxT,0);posMaxT[posMaxT==0]=np.nan
        negMaxT = agreeMSK(negMaxT,thresh,'N')*np.nanmean(negMaxT,0);negMaxT[negMaxT==0]=np.nan        
        posSM = agreeMSK(posSM,thresh,'N')*np.nanmean(posSM,0);posSM[posSM==0]=np.nan
        negSM = agreeMSK(negSM,thresh,'N')*np.nanmean(negSM,0);negSM[negSM==0]=np.nan
    
        #define color scales 

        if lvl >=500:
            gphMax = 200
            gphLvls = np.arange(-gphMax,gphMax+10,5)
            tScale = 3
            qScale = 10
            if (mode == 'IOD')|(mode == 'TAV (AMM)'):
                tScale = 2
                qScale = 10
                if (mode == 'TAV (AMM)'):
                    gphLvls = np.arange(-gphMax,gphMax+10,5)
                    qScale = 5
        else:
            gphMax = 200
            gphLvls = np.arange(-gphMax,gphMax+10,20) 
            tScale = 3
            qScale = 20
            if (mode == 'IOD')|(mode == 'TAV (AMM)'):
                tScale = 2
                qScale = 10
        
        
        pScale = 5
        smScale = 40
        omegaScale = .025
        
        
        fig = plt.figure();
        ax11=plt.subplot(231);ax12=plt.subplot(232);ax13=plt.subplot(233)
        ax21=plt.subplot(234);ax22=plt.subplot(235);ax23=plt.subplot(236);
        
        m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
        m12 = Basemap(projection='kav7',lon_0=0,ax=ax12); m12.drawcoastlines();m12.drawcountries()
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()
        m22 = Basemap(projection='kav7',lon_0=0,ax=ax22); m22.drawcoastlines();m22.drawcountries()
        m13 = Basemap(projection='kav7',lon_0=0,ax=ax13); m13.drawcoastlines();m13.drawcountries()
        m23 = Basemap(projection='kav7',lon_0=0,ax=ax23); m23.drawcoastlines();m23.drawcountries()
        ax11.set_title('Positive')
        ax12.set_title('Negative')
        ax13.set_title('Pos - Neg')
        
        im11=m11.pcolormesh(pLons, pLats,posPrcp,shading='flat',cmap=prcp_map,latlon=True,vmin=-pScale,vmax=pScale )
        cb12 = m12.colorbar(im11,"bottom", size="5%", pad="2%");
        im11=m11.pcolormesh(smLons, smLats,posSM,shading='flat',cmap='BrBG',latlon=True,vmin=-smScale,vmax=smScale )  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
        im12=m12.pcolormesh(pLons, pLats,negPrcp,shading='flat',cmap=prcp_map,latlon=True,vmin=-pScale,vmax=pScale )  
        im12=m12.pcolormesh(smLons,smLats,negSM,shading='flat',cmap='BrBG',latlon=True,vmin=-smScale,vmax=smScale )  
        im13=m13.pcolormesh(pLons, pLats,posPrcp-negPrcp,shading='flat',cmap=prcp_map,latlon=True,vmin=-pScale,vmax=pScale )
        im13=m13.pcolormesh(smLons, smLats,posSM-negSM,shading='flat',cmap='BrBG',latlon=True,vmin=-smScale,vmax=smScale )  
    
        im11=m11.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2],posVwnd[::2,::2],latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=2.5)
        im12=m12.quiver(gphLons[::2,::2],gphLats[::2,::2],negUwnd[::2,::2],negVwnd[::2,::2],latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=2.5)
        im13=m13.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd_all[::2,::2]-negUwnd_all[::2,::2],posVwnd[::2,::2]-negVwnd[::2,::2],latlon=True, color='grey',pivot='middle', units='inches',scale=qScale,headwidth=2.5)
        im13=m13.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2]-negUwnd[::2,::2],posVwnd[::2,::2]-negVwnd[::2,::2],latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=2.5)
     
    #    im11=m11.contour(gphLons, gphLats,posGPH_all,gphLvls, shading='flat',colors='darkgrey',latlon=True)   
    #    im12= m12.contour(gphLons, gphLats,negGPH_all,gphLvls, shading='flat',colors='darkgrey',latlon=True)
    #    im11=m11.contour(gphLons, gphLats,posGPH,gphLvls, shading='flat',colors='k',latlon=True)   
    #    im12= m12.contour(gphLons, gphLats,negGPH,gphLvls, shading='flat',colors='k',latlon=True)
    #    im13=m13.contour(gphLons, gphLats,posGPH_all-negGPH_all,gphLvls, shading='flat',colors='darkgrey',latlon=True)   
    #    im13=m13.contour(gphLons, gphLats,posGPH-negGPH,gphLvls, shading='flat',colors='k',latlon=True)   
    #   
        im21=m21.pcolormesh(gphLons, gphLats,posOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale )
        im22=m22.pcolormesh(gphLons, gphLats,negOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale )  
        im23=m23.pcolormesh(gphLons, gphLats,posOm-negOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale )  
    
        cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
        im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2],posVwnd[::2,::2],latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=2.5)
        im22=m22.quiver(gphLons[::2,::2],gphLats[::2,::2],negUwnd[::2,::2],negVwnd[::2,::2],latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=2.5)
        im23=m23.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2]-negUwnd[::2,::2],posVwnd[::2,::2]-negVwnd[::2,::2],latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=2.5)
       
        qk = plt.quiverkey(im22, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(50,20);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/event composites/'+mode+'_'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_ENSOyrs_yrOffset'+str(yrOffset)+saveNotes+'.png')
        plt.close()        
    
    
    
    
        fig = plt.figure();
        ax11=plt.subplot(231);ax12=plt.subplot(232);ax13=plt.subplot(233)
        ax21=plt.subplot(234);ax22=plt.subplot(235);ax23=plt.subplot(236);
        
        m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
        m12 = Basemap(projection='kav7',lon_0=0,ax=ax12); m12.drawcoastlines();m12.drawcountries()
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()
        m22 = Basemap(projection='kav7',lon_0=0,ax=ax22); m22.drawcoastlines();m22.drawcountries()
        m13 = Basemap(projection='kav7',lon_0=0,ax=ax13); m13.drawcoastlines();m13.drawcountries()
        m23 = Basemap(projection='kav7',lon_0=0,ax=ax23); m23.drawcoastlines();m23.drawcountries()
        ax11.set_title('Positive')
        ax12.set_title('Negative')
        ax13.set_title('Pos - Neg')
          
        im11=m11.contour(gphLons, gphLats,posGPH_all,gphLvls, shading='flat',colors='darkgrey',latlon=True)   
        im12= m12.contour(gphLons, gphLats,negGPH_all,gphLvls, shading='flat',colors='darkgrey',latlon=True)
        im11=m11.contour(gphLons, gphLats,posGPH,gphLvls, shading='flat',colors='k',latlon=True)   
        im12= m12.contour(gphLons, gphLats,negGPH,gphLvls, shading='flat',colors='k',latlon=True)
        im13=m13.contour(gphLons, gphLats,posGPH_all-negGPH_all,gphLvls, shading='flat',colors='darkgrey',latlon=True)   
        im13=m13.contour(gphLons, gphLats,posGPH-negGPH,gphLvls, shading='flat',colors='k',latlon=True)   
       
        im11=m11.pcolormesh(sstLons, sstLats,posSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-tScale,vmax=tScale )  
        im12=m12.pcolormesh(sstLons, sstLats,negSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-tScale,vmax=tScale )  
        im13=m13.pcolormesh(sstLons, sstLats,posSST-negSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-tScale,vmax=tScale )  
        cb12 = m12.colorbar(im12,"bottom", size="5%", pad="2%");

        im11=m11.pcolormesh(mLons, mLats,posMaxT,shading='flat',cmap='PuOr_r',latlon=True,vmin=-tScale,vmax=tScale )  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
        im12=m12.pcolormesh(mLons, mLats,negMaxT,shading='flat',cmap='PuOr_r',latlon=True,vmin=-tScale,vmax=tScale )  
        im13=m13.pcolormesh(mLons, mLats,posMaxT-negMaxT,shading='flat',cmap='PuOr_r',latlon=True,vmin=-tScale,vmax=tScale )  
        
        im21=m21.pcolormesh(gphLons, gphLats,posOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale )
        cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
        im22=m22.pcolormesh(gphLons, gphLats,negOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale )  
        im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2],posVwnd[::2,::2],latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=2.5)
        im22=m22.quiver(gphLons[::2,::2],gphLats[::2,::2],negUwnd[::2,::2],negVwnd[::2,::2],latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=2.5)
        im23=m23.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2]-negUwnd[::2,::2],posVwnd[::2,::2]-negVwnd[::2,::2],latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=2.5)
         
        qk = plt.quiverkey(im22, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(50,20);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/event composites/'+mode+'_'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_ENSOyrs_yrOffset'+str(yrOffset)+saveNotes+'_temp.png')
        plt.close()   