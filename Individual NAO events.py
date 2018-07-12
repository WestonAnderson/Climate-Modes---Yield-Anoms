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
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #

years = [1981,2011] #from the SVD analysis
#PCyrs = [1983,1987,1992,1993,1998,2010] #PC years to composite on
PCs = [1] #not python indexing
monCutoff = 1 #not python indexing
dayCutoff = '0'
nao_crops = 'wheat'
NAOseason = 'DJF'
numEvents = 4
lvl = 500
yrOffset = -1
timePeriod = '' #'' or ' 1960-2010'
notes = '_janSVD'
saveNotes = ''

#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#Read in agreement mask. agreeMSK(DF, frac, require agree w/ avg? (Y/N) )
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/agreement.py').read())
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/custom_div_cmap.py').read())
prcp_map = custom_div_cmap(20, mincol='saddlebrown', midcol='white' ,maxcol='CornflowerBlue')

yrList = np.array(range(years[0],years[1]+1))
seasons = ['DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ']
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr'+dayCutoff+'har'

#Read in the NAO time series to choose NAO years
NAOdf = pd.read_csv('/Volumes/Data_Archive/Data/NAO/nao_station_seasonal.csv',index_col=0)

#ENSO loading years. Chosen on MAM Nino3.4 / MEI
posCompYrs = NAOdf.loc[years[0]:years[1]-2][NAOseason].sort_values().index[-numEvents:]
negCompYrs = NAOdf.loc[years[0]:years[1]-2][NAOseason].sort_values().index[:numEvents]



#read in the netCDF file containing the variables
ncGeo = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/gph.nc',
                        'r',format='NETCDF4') #read in the netCDF file for the variable
geoVar=ncGeo.variables['hgt'][(ncGeo.variables['year'][:]>=years[0])&(ncGeo.variables['year'][:]<years[1]+2),
                              :,np.where(ncGeo.variables['level'][:]==lvl)[0][0],:,:]

geoVar2=ncGeo.variables['hgt'][(ncGeo.variables['year'][:]>=years[0])&(ncGeo.variables['year'][:]<years[1]+2),
                              :,np.where(ncGeo.variables['level'][:]==925)[0][0],:,:]

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
sstAnom_M = np.squeeze(sstAnom);sstAnom=np.ma.filled(sstAnom,0)


#Remove the climatology to get anomalies by month
for m in range(0,12):
    geoVar[:,m,...] = geoVar[:,m,...] - np.nanmean(geoVar[:,m,...],axis=0)
    geoVar2[:,m,...] = geoVar2[:,m,...] - np.nanmean(geoVar2[:,m,...],axis=0)
    omega[:,m,...] = omega[:,m,...] - np.nanmean(omega[:,m,...],axis=0)
    prcp[:,m,...] = prcp[:,m,...] - np.nanmean(prcp[:,m,...],axis=0)
    sm[:,m,...] = sm[:,m,...] - np.nanmean(sm[:,m,...],axis=0)
    Uwnd[:,m,...] = Uwnd[:,m,...] - np.nanmean(Uwnd[:,m,...],axis=0)
    Vwnd[:,m,...] = Vwnd[:,m,...] - np.nanmean(Vwnd[:,m,...],axis=0)
maxT[:,...] = maxT[:,...] - np.nanmean(maxT[:,...],axis=0)
    
#Rearrange to have only one time dimension
geoVar_M = geoVar.reshape([geoVar.shape[0]*geoVar.shape[1],geoVar.shape[2],geoVar.shape[3]])
geoVar2_M = geoVar2.reshape([geoVar2.shape[0]*geoVar2.shape[1],geoVar2.shape[2],geoVar2.shape[3]])
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

smLats=ncSM.variables['latitude'][:]
smLons=ncSM.variables['longitude'][:]#-180.
smLons, smLats = np.meshgrid(smLons,smLats)
mLats=ncMaxT.variables['latitude'][:]
mLons=ncMaxT.variables['longitude'][:]#-180.
mLons, mLats = np.meshgrid(mLons,mLats)




cPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/crop_'+NAOseason+nao_crops+
              '_'+ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy') 

for PC in PCs:
    cPC = cPCs[:,PC-1]
#    #Generate the list of years, transferring from PC years to NetCDF years
    PCyrs = list(zip(cPC,yrList))
    PCyrs.sort()
    negPCyrs = PCyrs[:numEvents]
    posPCyrs = PCyrs[-numEvents:]
    posYrInd = np.array([iY[1] == NCyrs for iY in posPCyrs])
    posYrInd = np.where(posYrInd.sum(0)==1)[0]
    negYrInd = np.array([iY[1] == NCyrs for iY in negPCyrs])
    negYrInd = np.where(negYrInd.sum(0)==1)[0]
    #composite by season
    for iS in np.arange(0,1,1):#np.arange(monCutoff-2,11+monCutoff,1):
        #create the index to select. Note variables have a year of padding on both sides
        # so the 0 year offset will start in year two of the data
        posYrInd = np.tile(np.arange((iS%12-1),(iS%12+2)),posYrInd.shape[0])+ \
                  np.repeat(posYrInd+yrOffset+1,3)*12
        negYrInd = np.tile(np.arange((iS%12-1),(iS%12+2)),negYrInd.shape[0])+ \
                  np.repeat(negYrInd+yrOffset+1,3)*12   
        limPos = posYrInd.shape[0]-np.ceil(np.sum(posYrInd>=sm_M.shape[0])/3)*3
        limNeg = negYrInd.shape[0]-np.ceil(np.sum(negYrInd>=sm_M.shape[0])/3)*3
        
        #index
        geoVar = geoVar_M[posYrInd,...]
        posGPH = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
        geoVar = geoVar_M[negYrInd,...]
        negGPH = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    
        geoVar2 = geoVar2_M[posYrInd,...]
        posGPH2 = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
        geoVar2 = geoVar2_M[negYrInd,...]
        negGPH2 = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    
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

        posGPH = agreeMSK(posGPH,2/3,'N')*np.nanmean(posGPH,0);posGPH[posGPH==0]=np.nan
        negGPH = agreeMSK(negGPH,2/3,'N')*np.nanmean(negGPH,0);negGPH[negGPH==0]=np.nan
        posGPH2 = agreeMSK(posGPH2,2/3,'N')*np.nanmean(posGPH2,0);posGPH2[posGPH2==0]=np.nan
        negGPH2 = agreeMSK(negGPH2,2/3,'N')*np.nanmean(negGPH2,0);negGPH2[negGPH2==0]=np.nan
        posOm = agreeMSK(posOm,2/3,'N')*np.nanmean(posOm,0);posOm[posOm==0]=np.nan
        negOm = agreeMSK(negOm,2/3,'N')*np.nanmean(negOm,0);negOm[negOm==0]=np.nan
        posUwnd = agreeMSK(posUwnd,2/3,'N')*np.nanmean(posUwnd,0);posUwnd[posUwnd==0]=np.nan
        negUwnd = agreeMSK(negUwnd,2/3,'N')*np.nanmean(negUwnd,0);negUwnd[negUwnd==0]=np.nan
        posVwnd = agreeMSK(posVwnd,2/3,'N')*np.nanmean(posVwnd,0);posVwnd[posVwnd==0]=np.nan
        negVwnd = agreeMSK(negVwnd,2/3,'N')*np.nanmean(negVwnd,0);negVwnd[negVwnd==0]=np.nan
        posPrcp = agreeMSK(posPrcp,2/3,'N')*np.nanmean(posPrcp,0);posPrcp[posPrcp==0]=np.nan
        negPrcp = agreeMSK(negPrcp,2/3,'N')*np.nanmean(negPrcp,0);negPrcp[negPrcp==0]=np.nan
        posSM = agreeMSK(posSM,2/3,'N')*np.nanmean(posSM,0);posSM[posSM==0]=np.nan
        negSM = agreeMSK(negSM,2/3,'N')*np.nanmean(negSM,0);negSM[negSM==0]=np.nan
        posSST = agreeMSK(posSST,2/3,'N')*np.nanmean(posSST,0);posSST[posSST==0]=np.nan
        negSST = agreeMSK(negSST,2/3,'N')*np.nanmean(negSST,0);negSST[negSST==0]=np.nan        
        posMaxT = agreeMSK(posMaxT,2/3,'N')*np.nanmean(posMaxT,0);posMaxT[posMaxT==0]=np.nan
        negMaxT = agreeMSK(negMaxT,2/3,'N')*np.nanmean(negMaxT,0);negMaxT[negMaxT==0]=np.nan        

        #define color scales    
        gphMax = 75
        gphLvls = np.arange(-gphMax,gphMax+10,10)           
        pScale = 5
        smScale = 40
        omegaScale = .01
        qScale = 5*35
        
        fig = plt.figure();
        ax11=plt.subplot(221);ax12=plt.subplot(222);
        ax21=plt.subplot(223);ax22=plt.subplot(224);
        
        m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
        m12 = Basemap(projection='kav7',lon_0=0,ax=ax12); m12.drawcoastlines();m12.drawcountries()
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()
        m22 = Basemap(projection='kav7',lon_0=0,ax=ax22); m22.drawcoastlines();m22.drawcountries()

        #im11=m11.pcolormesh(sstLons, sstLats,posSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-3,vmax=3, alpha=0.75)  
        #im12=m12.pcolormesh(sstLons, sstLats,negSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-3,vmax=3, alpha=0.75,)  
        #cb12 = m12.colorbar(im12,"bottom", size="5%", pad="2%");
             
        im11=m11.pcolormesh(pLons, pLats,posPrcp,shading='flat',cmap=prcp_map,latlon=True,vmin=-pScale,vmax=pScale, alpha=0.75)  
        im11=m11.pcolormesh(smLons, smLats,posSM,shading='flat',cmap='BrBG',latlon=True,vmin=-smScale,vmax=smScale, alpha=0.75)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
        im12=m12.pcolormesh(pLons, pLats,negPrcp,shading='flat',cmap=prcp_map,latlon=True,vmin=-pScale,vmax=pScale, alpha=0.75,)  
        im12=m12.pcolormesh(smLons, smLats,negSM,shading='flat',cmap='BrBG',latlon=True,vmin=-smScale,vmax=smScale, alpha=0.75,)  
        im11=m11.contour(gphLons, gphLats,posGPH,gphLvls, shading='flat',colors='k',latlon=True,alpha=0.75)   
        im12= m12.contour(gphLons, gphLats,negGPH,gphLvls, shading='flat',colors='k',latlon=True,alpha=0.75)
        
        im21=m21.pcolormesh(gphLons, gphLats,posOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale, alpha=0.75)
        cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
        im22=m22.pcolormesh(gphLons, gphLats,negOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale, alpha=0.75)  
        im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2],posVwnd[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1,scale=qScale, alpha=0.75)
        im22=m22.quiver(gphLons[::2,::2],gphLats[::2,::2],negUwnd[::2,::2],negVwnd[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=3,scale=qScale, alpha=0.75)
        
        qk = plt.quiverkey(im22, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(30,20);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/NAO/Crop/seasonal/composites/'+nao_crops+'_'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_PC'+str(PC)+'_yrOffset'+str(yrOffset)+timePeriod+notes+saveNotes+'.png')
        plt.close()        


        fig = plt.figure();
        ax11=plt.subplot(221);ax12=plt.subplot(222);
        ax21=plt.subplot(223);ax22=plt.subplot(224);
        
        m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
        m12 = Basemap(projection='kav7',lon_0=0,ax=ax12); m12.drawcoastlines();m12.drawcountries()
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()
        m22 = Basemap(projection='kav7',lon_0=0,ax=ax22); m22.drawcoastlines();m22.drawcountries()

        im11=m11.pcolormesh(sstLons, sstLats,posSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-3,vmax=3, alpha=0.75)  
        im12=m12.pcolormesh(sstLons, sstLats,negSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-3,vmax=3, alpha=0.75,)  
        cb12 = m12.colorbar(im12,"bottom", size="5%", pad="2%");
             
        im11=m11.pcolormesh(mLons, mLats,posMaxT,shading='flat',cmap='RdBu_r',latlon=True,vmin=-4,vmax=4, alpha=0.75)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
        im12=m12.pcolormesh(mLons, mLats,negMaxT,shading='flat',cmap='RdBu_r',latlon=True,vmin=-4,vmax=4, alpha=0.75,)  
        im11=m11.contour(gphLons, gphLats,posGPH2,gphLvls, shading='flat',colors='k',latlon=True,alpha=0.75)   
        im12= m12.contour(gphLons, gphLats,negGPH2,gphLvls, shading='flat',colors='k',latlon=True,alpha=0.75)
        
        im21=m21.pcolormesh(gphLons, gphLats,posOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale, alpha=0.75)
        cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
        im22=m22.pcolormesh(gphLons, gphLats,negOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale, alpha=0.75)  
        im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2],posVwnd[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1,scale=qScale, alpha=0.75)
        im22=m22.quiver(gphLons[::2,::2],gphLats[::2,::2],negUwnd[::2,::2],negVwnd[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=3,scale=qScale, alpha=0.75)
        
        qk = plt.quiverkey(im22, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(30,20);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/NAO/Crop/seasonal/composites/'+nao_crops+'_'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_PC'+str(PC)+'_yrOffset'+str(yrOffset)+timePeriod+notes+saveNotes+'_temp.png')
        plt.close()                    

#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#
        
           #    Repeat but composite on the NAO years   #
    
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*#

#composite by season
for iS in np.arange(0,1,1):#np.arange(monCutoff-2,11+monCutoff,1):
    #create the index to select. Note variables have a year of padding on both sides
    # so the 0 year offset will start in year two of the data
    posYrInd = np.array([iY == NCyrs for iY in posCompYrs])
    posYrInd = np.where(posYrInd.sum(0)==1)[0]
    negYrInd = np.array([iY == NCyrs for iY in negCompYrs])
    negYrInd = np.where(negYrInd.sum(0)==1)[0]
    posYrInd = np.tile(np.arange((iS%12-1),(iS%12+2)),posYrInd.shape[0])+ \
              np.repeat(posYrInd+yrOffset+1,3)*12
    negYrInd = np.tile(np.arange((iS%12-1),(iS%12+2)),negYrInd.shape[0])+ \
              np.repeat(negYrInd+yrOffset+1,3)*12   
    limPos = posYrInd.shape[0]-np.ceil(np.sum(posYrInd>=sm_M.shape[0])/3)*3
    limNeg = negYrInd.shape[0]-np.ceil(np.sum(negYrInd>=sm_M.shape[0])/3)*3
    
    #index
    geoVar = geoVar_M[posYrInd,...]
    posGPH = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    geoVar = geoVar_M[negYrInd,...]
    negGPH = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)

    geoVar2 = geoVar2_M[posYrInd,...]
    posGPH2 = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    geoVar2 = geoVar2_M[negYrInd,...]
    negGPH2 = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)

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

    posGPH = agreeMSK(posGPH,2/3,'N')*np.nanmean(posGPH,0);posGPH[posGPH==0]=np.nan
    negGPH = agreeMSK(negGPH,2/3,'N')*np.nanmean(negGPH,0);negGPH[negGPH==0]=np.nan
    posGPH2 = agreeMSK(posGPH2,2/3,'N')*np.nanmean(posGPH2,0);posGPH2[posGPH2==0]=np.nan
    negGPH2 = agreeMSK(negGPH2,2/3,'N')*np.nanmean(negGPH2,0);negGPH2[negGPH2==0]=np.nan
    posOm = agreeMSK(posOm,2/3,'N')*np.nanmean(posOm,0);posOm[posOm==0]=np.nan
    negOm = agreeMSK(negOm,2/3,'N')*np.nanmean(negOm,0);negOm[negOm==0]=np.nan
    posUwnd = agreeMSK(posUwnd,2/3,'N')*np.nanmean(posUwnd,0);posUwnd[posUwnd==0]=np.nan
    negUwnd = agreeMSK(negUwnd,2/3,'N')*np.nanmean(negUwnd,0);negUwnd[negUwnd==0]=np.nan
    posVwnd = agreeMSK(posVwnd,2/3,'N')*np.nanmean(posVwnd,0);posVwnd[posVwnd==0]=np.nan
    negVwnd = agreeMSK(negVwnd,2/3,'N')*np.nanmean(negVwnd,0);negVwnd[negVwnd==0]=np.nan
    posPrcp = agreeMSK(posPrcp,2/3,'N')*np.nanmean(posPrcp,0);posPrcp[posPrcp==0]=np.nan
    negPrcp = agreeMSK(negPrcp,2/3,'N')*np.nanmean(negPrcp,0);negPrcp[negPrcp==0]=np.nan
    posSM = agreeMSK(posSM,2/3,'N')*np.nanmean(posSM,0);posSM[posSM==0]=np.nan
    negSM = agreeMSK(negSM,2/3,'N')*np.nanmean(negSM,0);negSM[negSM==0]=np.nan
    posSST = agreeMSK(posSST,2/3,'N')*np.nanmean(posSST,0);posSST[posSST==0]=np.nan
    negSST = agreeMSK(negSST,2/3,'N')*np.nanmean(negSST,0);negSST[negSST==0]=np.nan        
    posMaxT = agreeMSK(posMaxT,2/3,'N')*np.nanmean(posMaxT,0);posMaxT[posMaxT==0]=np.nan
    negMaxT = agreeMSK(negMaxT,2/3,'N')*np.nanmean(negMaxT,0);negMaxT[negMaxT==0]=np.nan       

    #define color scales    
    gphMax = 75
    gphLvls = np.arange(-gphMax,gphMax+10,10)           
    pScale = 5
    smScale = 40
    omegaScale = .01
    qScale = 5*35
    
    fig = plt.figure();
    ax11=plt.subplot(221);ax12=plt.subplot(222);
    ax21=plt.subplot(223);ax22=plt.subplot(224);
    
    m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
    m12 = Basemap(projection='kav7',lon_0=0,ax=ax12); m12.drawcoastlines();m12.drawcountries()
    m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='kav7',lon_0=0,ax=ax22); m22.drawcoastlines();m22.drawcountries()

    #im11=m11.pcolormesh(sstLons, sstLats,posSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-3,vmax=3, alpha=0.75)  
    #im12=m12.pcolormesh(sstLons, sstLats,negSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-3,vmax=3, alpha=0.75,)  
    #cb12 = m12.colorbar(im12,"bottom", size="5%", pad="2%");
         
    im11=m11.pcolormesh(pLons, pLats,posPrcp,shading='flat',cmap=prcp_map,latlon=True,vmin=-pScale,vmax=pScale, alpha=0.75)  
    im11=m11.pcolormesh(smLons, smLats,posSM,shading='flat',cmap='BrBG',latlon=True,vmin=-smScale,vmax=smScale, alpha=0.75)  
    cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
    im12=m12.pcolormesh(pLons, pLats,negPrcp,shading='flat',cmap=prcp_map,latlon=True,vmin=-pScale,vmax=pScale, alpha=0.75,)  
    im12=m12.pcolormesh(smLons,smLats,negSM,shading='flat',cmap='BrBG',latlon=True,vmin=-smScale,vmax=smScale, alpha=0.75,)  
    im11=m11.contour(gphLons, gphLats,posGPH,gphLvls, shading='flat',colors='k',latlon=True,alpha=0.75)   
    im12= m12.contour(gphLons, gphLats,negGPH,gphLvls, shading='flat',colors='k',latlon=True,alpha=0.75)
    
    im21=m21.pcolormesh(gphLons, gphLats,posOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale, alpha=0.75)
    cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
    im22=m22.pcolormesh(gphLons, gphLats,negOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale, alpha=0.75)  
    im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2],posVwnd[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1,scale=qScale, alpha=0.75)
    im22=m22.quiver(gphLons[::2,::2],gphLats[::2,::2],negUwnd[::2,::2],negVwnd[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=3,scale=qScale, alpha=0.75)
    
    qk = plt.quiverkey(im22, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
    fig.set_size_inches(30,20);
    fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/NAO/Crop/seasonal/composites/'+nao_crops+'_'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_NAOyrs_yrOffset'+str(yrOffset)+timePeriod+notes+saveNotes+'.png')
    plt.close()        


    fig = plt.figure();
    ax11=plt.subplot(221);ax12=plt.subplot(222);
    ax21=plt.subplot(223);ax22=plt.subplot(224);
    
    m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
    m12 = Basemap(projection='kav7',lon_0=0,ax=ax12); m12.drawcoastlines();m12.drawcountries()
    m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()
    m22 = Basemap(projection='kav7',lon_0=0,ax=ax22); m22.drawcoastlines();m22.drawcountries()

    im11=m11.pcolormesh(sstLons, sstLats,posSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-3,vmax=3, alpha=0.75)  
    im12=m12.pcolormesh(sstLons, sstLats,negSST,shading='flat',cmap='RdBu_r',latlon=True,vmin=-3,vmax=3, alpha=0.75,)  
    cb12 = m12.colorbar(im12,"bottom", size="5%", pad="2%");
         
    im11=m11.pcolormesh(mLons, mLats,posMaxT,shading='flat',cmap='RdBu_r',latlon=True,vmin=-4,vmax=4, alpha=0.75)  
    cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
    im12=m12.pcolormesh(mLons, mLats,negMaxT,shading='flat',cmap='RdBu_r',latlon=True,vmin=-4,vmax=4, alpha=0.75,)  
    im11=m11.contour(gphLons, gphLats,posGPH2,gphLvls, shading='flat',colors='k',latlon=True,alpha=0.75)   
    im12= m12.contour(gphLons, gphLats,negGPH2,gphLvls, shading='flat',colors='k',latlon=True,alpha=0.75)
    
    im21=m21.pcolormesh(gphLons, gphLats,posOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale, alpha=0.75)
    cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
    im22=m22.pcolormesh(gphLons, gphLats,negOm,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-omegaScale,vmax=omegaScale, alpha=0.75)  
    im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2],posVwnd[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1,scale=qScale, alpha=0.75)
    im22=m22.quiver(gphLons[::2,::2],gphLats[::2,::2],negUwnd[::2,::2],negVwnd[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=3,scale=qScale, alpha=0.75)
    
    qk = plt.quiverkey(im22, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
    fig.set_size_inches(30,20);
    fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/NAO/Crop/seasonal/composites/'+nao_crops+'_'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_NAOyrs_yrOffset'+str(yrOffset)+timePeriod+notes+saveNotes+'_temp.png')
    plt.close()   