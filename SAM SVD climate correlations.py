#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 11:21:57 2018

@author: weston

Read in the EOFs and PCs from crop/SST matrices and correlate the PCs to various
physical quantities.
"""

import numpy as np
import pandas as pd
import time
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
start = time.clock()
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #

crops = ['maize','wheat','soy']#
years = [1981,2011]
PCs = [1] #not python indexing
monCutoff = 1 #not python indexing
lvl = 925 # geopotential height lvl
sigLvl = .1682 #.995, .1682     #vertical velocity
timePeriod =''# ' 1960-2010' #'' or ' 1960-2010'
notes = '_janSVD_allTrops'

#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

## colormap function
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/custom_div_cmap.py').read())
prcp_map = custom_div_cmap(20, mincol='saddlebrown', midcol='white' ,maxcol='CornflowerBlue')
##


seasons = ['DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ']
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr0har'

#~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~#
#read in the netCDF file containing the variables

#read in the netCDF file containing the variables
ncGeo = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/gph.nc',
                        'r',format='NETCDF4') #read in the netCDF file for the variable
geoVar=ncGeo.variables['hgt'][(ncGeo.variables['year'][:]>=years[0]-1)&(ncGeo.variables['year'][:]<years[1]+2),
                              :,np.where(ncGeo.variables['level'][:]==lvl)[0][0],:,:]

ncSM = netCDF4.Dataset('/Volumes/Data_Archive/Data/SoilMoisture/GLDAS/Noah10_nc/detrended/SoilM.nc','r',format='NETCDF4') 
sm = np.nansum(ncSM.variables['SoilM'][(ncSM.variables['year'][:]>=years[0]-1)&(ncSM.variables['year'][:]<years[1]+2),:,:3,:,:],2)

ncMaxT = netCDF4.Dataset('//Volumes/Data_Archive/Data/Temp/BerkeleyEarth/monthly/Complete_TMAX_LatLong1.nc','r',format='NETCDF4') #read in the netCDF file for the variable
maxT = ncMaxT.variables['temperature'][(ncMaxT.variables['time'][:]>=years[0]-1)&(ncMaxT.variables['time'][:]<years[1]+2),:,:] 

ncUwnd = netCDF4.Dataset('/Volumes/Data_Archive/Data/Wind/detrend/uwnd.nc','r',format='NETCDF4') #read in the netCDF file for the variable
Uwnd=ncUwnd.variables['uwnd'][(ncUwnd.variables['year'][:]>=years[0]-1)&(ncUwnd.variables['year'][:]<years[1]+2),
                              :,np.where(ncUwnd.variables['level'][:]==lvl)[0][0],:,:]
ncVwnd = netCDF4.Dataset('/Volumes/Data_Archive/Data/Wind/detrend/vwnd.nc','r',format='NETCDF4') #read in the netCDF file for the variable
Vwnd=ncVwnd.variables['vwnd'][(ncVwnd.variables['year'][:]>=years[0]-1)&(ncVwnd.variables['year'][:]<years[1]+2),
                              :,np.where(ncVwnd.variables['level'][:]==lvl)[0][0],:,:]
ncOmega = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/omega.nc','r',format='NETCDF4')
omega= ncOmega.variables['omega'][(ncOmega.variables['year'][:]>=years[0]-1)&(ncOmega.variables['year'][:]<years[1]+2),
                              :,np.where(ncOmega.variables['level'][:]==lvl)[0][0],:,:]
#ncSM = netCDF4.Dataset('/Volumes/Data_Archive/Data/Precip/GPCC_nc/detrended/sumP.nc','r',format='NETCDF4') #read in the netCDF file for the variable

ncChi = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/chi.nc','r',format='NETCDF4')
chi = ncChi.variables['chi'][(ncChi.variables['year'][:]>=years[0]-1)&(ncChi.variables['year'][:]<years[1]+2),
                              :,np.where(ncChi.variables['level'][:]==sigLvl)[0][0],:,:]
ncP = netCDF4.Dataset('/Volumes/Data_Archive/Data/Precip/GPCP/detrended/GPCP_precip.nc','r',format='NETCDF4')
prcp = ncP.variables['precip'][(ncP.variables['year'][:]>=years[0]-1)&(ncP.variables['year'][:]<years[1]+2),:,:,:]


#Read in the SST anomaly data
sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
#T is months since 1960-01-01, and I want to start at 1950 with one year lag
sstLats=sstData.variables['Y'][:]
sstLons=sstData.variables['X'][:]
sstLons,sstLats=np.meshgrid(sstLons,sstLats)

sstAnom=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-years[0]+1))&
                        (sstData.variables['T'][:]<=12*(2+years[1]-1960)),...]
sstAnom = np.squeeze(sstAnom);sstAnom=np.ma.filled(sstAnom,0)


#Remove the climatology to get anomalies by month
for m in range(0,12):
    geoVar[:,m,...] = geoVar[:,m,...] - np.nanmean(geoVar[:,m,...],axis=0)
    omega[:,m,...] = omega[:,m,...] - np.nanmean(omega[:,m,...],axis=0)
    chi[:,m,...] = chi[:,m,...] - np.nanmean(chi[:,m,...],axis=0)
    prcp[:,m,...] = prcp[:,m,...] - np.nanmean(prcp[:,m,...],axis=0)
    sm[:,m,...] = sm[:,m,...] - np.nanmean(sm[:,m,...],axis=0)
    Uwnd[:,m,...] = Uwnd[:,m,...] - np.nanmean(Uwnd[:,m,...],axis=0)
    Vwnd[:,m,...] = Vwnd[:,m,...] - np.nanmean(Vwnd[:,m,...],axis=0)
maxT[:,...] = maxT[:,...] - np.nanmean(maxT[:,...],axis=0)
    
#Rearrange to have only one time dimension
geoVar = geoVar.reshape([geoVar.shape[0]*geoVar.shape[1],geoVar.shape[2],geoVar.shape[3]])
omega = omega.reshape([omega.shape[0]*omega.shape[1],omega.shape[2],omega.shape[3]])
prcp = prcp.reshape([prcp.shape[0]*prcp.shape[1],prcp.shape[2],prcp.shape[3]])
chi = chi.reshape([chi.shape[0]*chi.shape[1],chi.shape[2],chi.shape[3]])
sm= sm.reshape([sm.shape[0]*sm.shape[1],sm.shape[2],sm.shape[3]])
#maxT = maxT.reshape([maxT.shape[0]*maxT.shape[1],maxT.shape[2],maxT.shape[3]])
Uwnd = Uwnd.reshape([Uwnd.shape[0]*Uwnd.shape[1],Uwnd.shape[2],Uwnd.shape[3]])
Vwnd = Vwnd.reshape([Vwnd.shape[0]*Vwnd.shape[1],Vwnd.shape[2],Vwnd.shape[3]])

  
NCyrs =ncGeo.variables['year'][(ncGeo.variables['year'][:]>=years[0])&(ncGeo.variables['year'][:]<years[1]+2)]
gphLats=ncGeo.variables['latitude'][:]
gphLons=ncGeo.variables['longitude'][:]#-180.
gphLons, gphLats = np.meshgrid(gphLons,gphLats)
smLats=ncSM.variables['latitude'][:]
smLons=ncSM.variables['longitude'][:]#-180.
smLons, smLats = np.meshgrid(smLons,smLats)

mLats=ncMaxT.variables['latitude'][:]
mLons=ncMaxT.variables['longitude'][:]#-180.
mLons, mLats = np.meshgrid(mLons,mLats)

pLats=ncP.variables['latitude'][:]
pLons=ncP.variables['longitude'][:]#-180.
pLons, pLats = np.meshgrid(pLons,pLats)

chiLats=ncChi.variables['latitude'][:]
chiLons=ncChi.variables['longitude'][:]#-180.
chiLons, chiLats = np.meshgrid(chiLons,chiLats)

cPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/SAM/crop_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy') 
sstPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/SAM/sst_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy')             


#~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~#           
for PC in PCs:
    cPC = cPCs[:,PC-1]
    sstPC = sstPCs[:,PC-1]
    #Regressions by season
    for iS in np.arange(10,14,1):#np.arange(-6,21,3):
        yrOffset = iS//12
        #FIRST Correlate the time expansion coefficients with different ENSO indices
        SAM = pd.read_csv('/Volumes/Data_Archive/Data/SAM/monthly.sam.1957.2018.csv',index_col=0)
        ind = (np.arange(0,years[1]-years[0]+1)+1)*12+iS 
        SAM  = SAM.loc[years[0]-1:years[1]+1]
        SAM = np.ravel(SAM)
        if (PC==1):
            fig = plt.figure();ax=plt.subplot(111)
            ax.plot(range(years[0],years[1]+1),SAM[ind]/np.std(SAM[ind]),'b',label='SAM')
            ax.plot(range(years[0],years[1]+1),sstPC/np.std(sstPC),'k',label='SST PC'+str(PC))
            ax.plot(range(years[0],years[1]+1),cPC/np.std(cPC),'--k',label='crop yield PC'+str(PC))
            ax.text(2005,2.7,'c-SAM: '+str(np.round(np.corrcoef(cPC,SAM[ind])[0][1],4)))
            ax.text(2005,2.5,'gph-SAM: '+str(np.round(np.corrcoef(sstPC,SAM[ind])[0][1],4)))
            ax.text(2005,2.3,'c-sst: '+str(np.round(np.corrcoef(cPC,sstPC)[0][1],4)))
            ax.legend()
            fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/SAM/SST/SAM Index/sstPC'+str(PC)+ENSOyr+'_'+seasons[iS%12]+'_yrOffset'+str(yrOffset)+timePeriod+notes+'.png')
            plt.close()     
        continue
        #Now make the full correlation maps       
        #create the index to select. Note variables have a year of padding on both sides
        # so the 0 year offset will start in year two of the data
        ind = np.tile(np.arange((iS%12-1),(iS%12+2)),sstPC.shape[0])+ \
              np.repeat(np.arange(0,years[1]-years[0]+1)+yrOffset+1,3)*12
        
        X = sstPC        
        Y = sstAnom[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1)) #average across the season
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]]) #unravel to a vector
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1)) #correlate
        sst_r = r.reshape([sstAnom.shape[1],sstAnom.shape[2]]) #reshape

        Y = geoVar[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        gph_r = r.reshape([geoVar.shape[1],geoVar.shape[2]])

        Y = sm[ind[:(ind.shape[0]-(3*(yrOffset+1)))],...] #SM only goes to 2010, needs a special index
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X[:Y.shape[0]]) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        sm_r = r.reshape([sm.shape[1],sm.shape[2]])

        Y = maxT[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        maxT_r = r.reshape([maxT.shape[1],maxT.shape[2]])

        Y = omega[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        omega_r = r.reshape([omega.shape[1],omega.shape[2]])

        Y = Uwnd[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        Uwnd_r = r.reshape([Uwnd.shape[1],Uwnd.shape[2]])
        
        Y = Vwnd[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        Vwnd_r = r.reshape([Vwnd.shape[1],Vwnd.shape[2]])        
 
        fig = plt.figure();
        ax11=plt.subplot(211);
        ax21=plt.subplot(212);
        m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()             
        im11=m11.pcolormesh(mLons,mLats,maxT_r,shading='flat',cmap='PuOr_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
        im11=m11.pcolormesh(sstLons, sstLats,sst_r,shading='flat',cmap='RdBu_r',latlon=True,vmin=-1.1,vmax=1.1)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
        im11=m11.contour(gphLons, gphLats,gph_r,np.arange(-1,1.1,.1), shading='flat',colors='k',latlon=True,alpha=0.75)          
        im21=m21.pcolormesh(gphLons, gphLats,omega_r,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)
        cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
        im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],Uwnd_r[::2,::2],Vwnd_r[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1, alpha=0.75)
        qk = plt.quiverkey(im21, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(30,20);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/SAM/SST/seasonal/correlations/'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_PC'+str(PC)+'_yrOffset'+str(yrOffset)+timePeriod+notes+'_temp.png')
        plt.close()

        fig = plt.figure();
        ax11=plt.subplot(211);
        ax21=plt.subplot(212);
        m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()
        im11=m11.pcolormesh(sstLons, sstLats,sst_r,shading='flat',cmap='RdBu_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");             
        im11=m11.pcolormesh(smLons,smLats,sm_r,shading='flat',cmap='RdBu_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
        im11=m11.contour(gphLons, gphLats,gph_r,np.arange(-1,1.1,.1), shading='flat',colors='k',latlon=True,alpha=0.75)          
        im21=m21.pcolormesh(gphLons, gphLats,omega_r,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)
        cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
        im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],Uwnd_r[::2,::2],Vwnd_r[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1, alpha=0.75)
        qk = plt.quiverkey(im21, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(30,20);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/SAM/SST/seasonal/correlations/'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_PC'+str(PC)+'_yrOffset'+str(yrOffset)+timePeriod+notes+'.png')
        plt.close()

        X = cPC
        Y = sstAnom[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1)) #average across the season
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]]) #unravel to a vector
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1)) #correlate
        sst_r = r.reshape([sstAnom.shape[1],sstAnom.shape[2]]) #reshape
        #sst_r[(sst_r<.3)&(sst_r>-.3)]=0
        
        Y = geoVar[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        gph_r = r.reshape([geoVar.shape[1],geoVar.shape[2]])
        #gph_r[(gph_r<.3)&(gph_r>-.3)]=0
        
        lim = ind.shape[0]-np.ceil(np.sum(ind>=sm.shape[0])/3)*3
        Y = sm[ind[:lim],...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X[:Y.shape[0]]) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        sm_r = r.reshape([sm.shape[1],sm.shape[2]])
        #sm_r[(sm_r<.3)&(sm_r>-.3)]=0
        
        Y = maxT[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        maxT_r = r.reshape([maxT.shape[1],maxT.shape[2]])
        #maxT_r[(maxT_r<.3)&(maxT_r>-.3)]=0
        
        Y = omega[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        omega_r = r.reshape([omega.shape[1],omega.shape[2]])
        #omega_r[(omega_r<.3)&(omega_r>-.3)]=0
        
        Y = prcp[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        prcp_r = r.reshape([prcp.shape[1],prcp.shape[2]])
        #prcp_r[(prcp_r<.3)&(prcp_r>-.3)]=0
        
        Y = chi[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        chi_r = r.reshape([chi.shape[1],chi.shape[2]])
        #chi_r[(chi_r<.3)&(chi_r>-.3)]=0
        
        Y = Uwnd[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        Uwnd_r = r.reshape([Uwnd.shape[1],Uwnd.shape[2]])
        #Uwnd_r[Uwnd_r<.3]=0
        
        Y = Vwnd[ind,...] #select months and years
        Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
        Y = np.squeeze(np.nanmean(Y,1))
        Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]])
        r = Y.T.dot(X) / (np.std(X,0)*np.std(Y,0)*(Y.shape[0]-1))
        Vwnd_r = r.reshape([Vwnd.shape[1],Vwnd.shape[2]])   
        #Vwnd_r[Vwnd_r<.3]=0
        
        Qscale=3
        
        fig = plt.figure();
        ax11=plt.subplot(211);
        ax21=plt.subplot(212);
        m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()          
        im11=m11.pcolormesh(mLons,mLats,maxT_r,shading='flat',cmap='PuOr_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
        im11=m11.pcolormesh(sstLons, sstLats,sst_r,shading='flat',cmap='RdBu_r',latlon=True,vmin=-1.1,vmax=1.1)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");   
        im11=m11.contour(gphLons, gphLats,gph_r,np.arange(-1,1.1,.1), shading='flat',colors='k',latlon=True,alpha=0.75)          
        im21=m21.pcolormesh(gphLons, gphLats,omega_r,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)
        cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
        im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],Uwnd_r[::2,::2],Vwnd_r[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1,scale=5*35, alpha=0.75)
        qk = plt.quiverkey(im21, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(30,20);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/SAM/Crop/seasonal/correlations/'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_PC'+str(PC)+'_yrOffset'+str(yrOffset)+timePeriod+notes+'_temp.png')
        plt.close()

        fig = plt.figure();
        ax11=plt.subplot(211);
        ax21=plt.subplot(212);
        m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()
        im11=m11.pcolormesh(pLons, pLats,prcp_r,shading='flat',cmap=prcp_map,latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");             
        im11=m11.pcolormesh(smLons,smLats,sm_r,shading='flat',cmap='BrBG',latlon=True,vmin=-1.1,vmax=1.1, alpha=1)  
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
        im11=m11.quiver(gphLons[::2,::2],gphLats[::2,::2],Uwnd_r[::2,::2],Vwnd_r[::2,::2],latlon=True, pivot='middle', units='inches',scale=Qscale,headwidth=2.5)
        im21=m21.pcolormesh(gphLons, gphLats,omega_r,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)
        cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
        im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],Uwnd_r[::2,::2],Vwnd_r[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1,scale=5*35, alpha=0.75)
        qk = plt.quiverkey(im21, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(30,20);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/SAM/Crop/seasonal/correlations/'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_PC'+str(PC)+'_yrOffset'+str(yrOffset)+timePeriod+notes+'.png')
        plt.close()


            