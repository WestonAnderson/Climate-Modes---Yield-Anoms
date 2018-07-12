#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 11:21:57 2018

@author: weston

Read in the EOFs and PCs from crop/SST matrices and perform a partial correlation 
between the PCs and various physical quantities.
"""

import numpy as np
import pandas as pd
import time
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits
start = time.clock()
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #

years = [1981,2011]
monCutoff = 1 #not python indexing
notes = '_janSVD_allTrops'
pcType = 'SST' #'Crop' or 'SST' or 'Index'

#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

## colormap function
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/custom_div_cmap.py').read())
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/partial_corr_ind_dep.py').read())
prcp_map = custom_div_cmap(20, mincol='saddlebrown', midcol='white' ,maxcol='CornflowerBlue')
##

cLvls = np.arange(-1,1.1,.1)

seasons = ['DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ']
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr0har'



#~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~#
#read in the netCDF file containing the variables

#read in the netCDF file containing the variables
ncGeo = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/gph.nc',
                        'r',format='NETCDF4') #read in the netCDF file for the variable
geoVar=ncGeo.variables['hgt'][(ncGeo.variables['year'][:]>=years[0]-1)&(ncGeo.variables['year'][:]<years[1]+2),
                              :,np.where(ncGeo.variables['level'][:]==200)[0][0],:,:]

ncSM = netCDF4.Dataset('/Volumes/Data_Archive/Data/SoilMoisture/GLDAS/Noah10_nc/detrended/SoilM.nc','r',format='NETCDF4') 
sm = np.nansum(ncSM.variables['SoilM'][(ncSM.variables['year'][:]>=years[0]-1)&(ncSM.variables['year'][:]<years[1]+2),:,:3,:,:],2)

ncMaxT = netCDF4.Dataset('//Volumes/Data_Archive/Data/Temp/BerkeleyEarth/monthly/Complete_TMAX_LatLong1.nc','r',format='NETCDF4') #read in the netCDF file for the variable
maxT = ncMaxT.variables['temperature'][(ncMaxT.variables['time'][:]>=years[0]-1)&(ncMaxT.variables['time'][:]<years[1]+2),:,:] 

ncUwnd = netCDF4.Dataset('/Volumes/Data_Archive/Data/Wind/detrend/uwnd.nc','r',format='NETCDF4') #read in the netCDF file for the variable
Uwnd=ncUwnd.variables['uwnd'][(ncUwnd.variables['year'][:]>=years[0]-1)&(ncUwnd.variables['year'][:]<years[1]+2),
                              :,np.where(ncUwnd.variables['level'][:]==925)[0][0],:,:]
ncVwnd = netCDF4.Dataset('/Volumes/Data_Archive/Data/Wind/detrend/vwnd.nc','r',format='NETCDF4') #read in the netCDF file for the variable
Vwnd=ncVwnd.variables['vwnd'][(ncVwnd.variables['year'][:]>=years[0]-1)&(ncVwnd.variables['year'][:]<years[1]+2),
                              :,np.where(ncVwnd.variables['level'][:]==925)[0][0],:,:]
#ncSM = netCDF4.Dataset('/Volumes/Data_Archive/Data/Precip/GPCC_nc/detrended/sumP.nc','r',format='NETCDF4') #read in the netCDF file for the variable

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
    prcp[:,m,...] = prcp[:,m,...] - np.nanmean(prcp[:,m,...],axis=0)
    sm[:,m,...] = sm[:,m,...] - np.nanmean(sm[:,m,...],axis=0)
    Uwnd[:,m,...] = Uwnd[:,m,...] - np.nanmean(Uwnd[:,m,...],axis=0)
    Vwnd[:,m,...] = Vwnd[:,m,...] - np.nanmean(Vwnd[:,m,...],axis=0)
maxT[:,...] = maxT[:,...] - np.nanmean(maxT[:,...],axis=0)
maxT = maxT.filled(fill_value=0)
sm = sm.filled(fill_value=0)
#Rearrange to have only one time dimension
geoVar = geoVar.reshape([geoVar.shape[0]*geoVar.shape[1],geoVar.shape[2],geoVar.shape[3]])
prcp = prcp.reshape([prcp.shape[0]*prcp.shape[1],prcp.shape[2],prcp.shape[3]])
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


if pcType is 'Crop':
    iodPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/crop_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy') 
    tavPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/crop_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy') 
    naoPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/NAO/crop_DJFwheat_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+'_janSVD.npy') 
    ensoPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/all crops/crop_all crops_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy')    
elif pcType is 'SST':
    iodPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/sst_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy') 
    tavPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/sst_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy') 
    naoPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/NAO/gph_DJFwheat_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+'_janSVD.npy') 
    ensoPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/all crops/sst_all crops_'+
              ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy')              
elif pcType is 'Index':
    iod = pd.read_csv('/Volumes/Data_Archive/Data/IOD/dmi_monthly.csv',index_col=0)
    iod = iod.loc[years[0]:years[1]][['Jul','Aug','Sep']].mean(1)
    enso = pd.read_csv('/Volumes/Data_Archive/Data/ENSO/NinoIndex34.csv',index_col=0,skiprows=1)
    enso = enso.loc[years[0]:years[1]]['12']#'OND' season
    amm = pd.read_csv('/Volumes/Data_Archive/Data/AMM/amm.csv',index_col=0)
    amm = amm.loc[years[0]:years[1]][['Apr','May','Jun','Jul']].mean(1)
    tsa = pd.read_csv('/Volumes/Data_Archive/Data/TSA/TSA.csv',index_col=0)
    tsa = tsa.loc[years[0]:years[1]][['Apr','May','Jun','Jul']].mean(1)
    nao =pd.read_csv('/Volumes/Data_Archive/Data/NAO/nao_station_seasonal.csv',index_col=0)
    nao.index = nao.index.values-1 #because I'm using the DJF season, and it's counted in the begining of the year
    nao = nao.loc[years[0]:years[1]]['DJF']
else:
    print('error w/ PC type')
    exit

#~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~#           
if (pcType == 'Crop')|(pcType=='SST'):
    modes = ['NAO','IOD','ENSO','AMM']
    pcs = np.vstack([naoPCs[:,0],iodPCs[:,0],-ensoPCs[:,0]+ensoPCs[:,1],tavPCs[:,0]-tavPCs[:,1]]) 
elif pcType is 'Index':
    modes = ['NAO','DMI','ENSO','AMM','TSA']
    pcs = np.vstack([nao,iod,enso,amm,tsa])
    
corners = {'ENSO':[-60,60,-179,179],'NAO':[20,80,-50,50],'AMM':[-60,30,-80,30],'IOD':[-40,30,20,160]}

#Regressions by season
for iS in [0]:#[-5,-3,-1,0,1,3,5,7,11,12,13,15]:
    yrOffset = iS//12

    #Now make the full correlation maps       
    #create the index to select. Note variables have a year of padding on both sides
    # so the 0 year offset will start in year two of the data
    ind = np.tile(np.arange((iS%12-1),(iS%12+2)),pcs.shape[1])+ \
          np.repeat(np.arange(0,years[1]-years[0]+1)+yrOffset+1,3)*12

    X = pcs.T
    
    Y = sstAnom[ind,...] #select months and years
    Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
    Y = np.squeeze(np.nanmean(Y,1)) #average across the season
    Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]]) #unravel to a vector
    
    R = partial_corr(Y,X)
    sst_r = R.T.reshape([R.shape[1],sstAnom.shape[1],sstAnom.shape[2]]) #reshape
    sst_t = sst_r*np.sqrt((pcs.shape[1]-2-(R.shape[1]-1)/(1-sst_r**2)))
    sst_r[sst_r==0]=np.nan
    sst_stpl = sst_r.copy()
    sst_stpl[np.abs(sst_t)>1.696]=np.nan #1.696 is the t statistic for 90% (two tailed) significance in a time series of this length

    Y = geoVar[ind,...] #select months and years
    Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
    Y = np.squeeze(np.nanmean(Y,1)) #average across the season
    Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]]) #unravel to a vector
    R = partial_corr(Y,X)
    gph_r = R.T.reshape([R.shape[1],geoVar.shape[1],geoVar.shape[2]]) #reshape
    gph_t = gph_r*np.sqrt((pcs.shape[1]-2-(R.shape[1]-1)/(1-gph_r**2)))
    gph_all = gph_r.copy()
    gph_r[np.abs(gph_t)<1.696]=np.nan
    
    Y = maxT[ind,...] #select months and years
    Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
    Y = np.squeeze(np.nanmean(Y,1)) #average across the season
    Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]]) #unravel to a vector
    R = partial_corr(Y,X)
    maxT_r = R.T.reshape([R.shape[1],maxT.shape[1],maxT.shape[2]]) #reshape
    maxT_t = maxT_r*np.sqrt((pcs.shape[1]-2-(R.shape[1]-1)/(1-maxT_r**2)))
    maxT_r[maxT_r==0]=np.nan
    maxT_stpl = maxT_r.copy()
    maxT_stpl[np.abs(maxT_t)>1.696]=np.nan
 
    lim = np.int(ind.shape[0]-np.ceil(np.sum(ind>=sm.shape[0])/3)*3)
    Y = sm[ind[:lim],...] #select months and years
    Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
    Y = np.squeeze(np.nanmean(Y,1)) #average across the season
    Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]]) #unravel to a vector
    R = partial_corr(Y,X[:Y.shape[0]])
    sm_r = R.T.reshape([R.shape[1],sm.shape[1],sm.shape[2]]) #reshape
    sm_t = sm_r*np.sqrt((pcs.shape[1]-2-(R.shape[1]-1)/(1-sm_r**2)))
    sm_r[sm_r==0]=np.nan
    sm_stpl = sm_r.copy()
    sm_stpl[np.abs(sm_t)>1.696]=np.nan

    Y = prcp[ind,...] #select months and years
    Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
    Y = np.squeeze(np.nanmean(Y,1)) #average across the season
    Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]]) #unravel to a vector
    R = partial_corr(Y,X)
    prcp_r = R.T.reshape([R.shape[1],prcp.shape[1],prcp.shape[2]]) #reshape
    oceans = mpl_toolkits.basemap.maskoceans(pLons-180,pLats,prcp_r[0,...])
    lndMsk = oceans.mask*1
    lndMsk,lons = mpl_toolkits.basemap.shiftgrid(0,lndMsk,pLons[0,:]-180)
    prcp_r = prcp_r*lndMsk
    prcp_t = prcp_r*np.sqrt((pcs.shape[1]-2-(R.shape[1]-1)/(1-prcp_r**2)))
    prcp_r[prcp_r==0]=np.nan
    prcp_stpl=prcp_r.copy()
    prcp_stpl[np.abs(prcp_t)>1.696]=np.nan
    

    Y = Uwnd[ind,...] #select months and years
    Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
    Y = np.squeeze(np.nanmean(Y,1)) #average across the season
    Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]]) #unravel to a vector
    R = partial_corr(Y,X)
    Uwnd_r = R.T.reshape([R.shape[1],Uwnd.shape[1],Uwnd.shape[2]]) #reshape
    Uwnd_t = Uwnd_r*np.sqrt((pcs.shape[1]-2-(R.shape[1]-1)/(1-Uwnd_r**2)))
    Uwnd_all = Uwnd_r.copy()
    
    Y = Vwnd[ind,...] #select months and years
    Y = Y.reshape([Y.shape[0]//3,3,Y.shape[1],Y.shape[2]]) #reindex
    Y = np.squeeze(np.nanmean(Y,1)) #average across the season
    Y = np.reshape(Y,[Y.shape[0],Y.shape[1]*Y.shape[2]]) #unravel to a vector
    R = partial_corr(Y,X)
    Vwnd_r = R.T.reshape([R.shape[1],Vwnd.shape[1],Vwnd.shape[2]]) #reshape
    Vwnd_t = Vwnd_r*np.sqrt((pcs.shape[1]-2-(R.shape[1]-1)/(1-Vwnd_r**2)))
    Vwnd_all = Vwnd_r.copy()
    Vwnd_r[(np.abs(Vwnd_t)<1.696)&(np.abs(Uwnd_t)<1.696)]=np.nan
    Uwnd_r[(np.abs(Vwnd_t)<1.696)&(np.abs(Uwnd_t)<1.696)]=np.nan

#    
    
    Qscale=3
    
    for iM in range(np.size(modes)):
        
        fig = plt.figure();
        ax11=plt.subplot(111);
        m11 = Basemap(projection='cyl',llcrnrlat=corners[modes[iM]][0],urcrnrlat=corners[modes[iM]][1],
                      llcrnrlon=corners[modes[iM]][2],urcrnrlon=corners[modes[iM]][3],resolution='c',ax=ax11)
        m11.drawcoastlines();m11.drawcountries()
        sst_r_temp = np.append(sst_r[iM,:,90:],sst_r[iM,:,:90],1)
        im11=plt.contourf(sstLons-180, sstLats,sst_r_temp,cLvls,shading='flat',cmap='RdBu_r',latlon=True)  
        lonStpl,latStpl = m11(sstLons[np.isfinite(sst_stpl[iM,...])]-180,sstLats[np.isfinite(sst_stpl[iM,...])]); 
        m11.plot(lonStpl,latStpl,'ko',markersize=.5)
#        cb11 = m11.colorbar(im11,"right"); 
#        cb11.ax.set_yticklabels(cb11.ax.get_yticklabels(), fontsize=24)

#        if modes[iM] is 'NAO':
#            im11=m11.pcolormesh(mLons,mLats,maxT_r[iM,...],shading='flat',cmap='PuOr_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)  
#            lonStpl,latStpl = m11(mLons[np.isfinite(maxT_stpl[iM,...])],mLats[np.isfinite(maxT_stpl[iM,...])]); 
#            m11.plot(lonStpl,latStpl,'ko',markersize=.2)
#        else:
        im11=plt.contourf(smLons,smLats,sm_r[iM,...],cLvls,shading='flat',cmap='BrBG',latlon=True)  
        lonStpl,latStpl = m11(smLons[np.isfinite(sm_stpl[iM,...])],smLats[np.isfinite(sm_stpl[iM,...])]); 
        m11.plot(lonStpl,latStpl,'ko',markersize=.25)
    
        cb11 = m11.colorbar(im11,"right"); 
        cb11.ax.set_yticklabels(cb11.ax.get_yticklabels(), fontsize=24)
        cb11.ax.set_yticklabels(cb11.ax.get_yticklabels(), fontsize=20)
        im11=m11.quiver(gphLons[::2,::2],gphLats[::2,::2],Uwnd_all[iM,::2,::2],Vwnd_all[iM,::2,::2],
                        latlon=True, color='0.5',pivot='middle', units='inches',scale=Qscale,headwidth=1.5)
        im11=m11.quiver(gphLons[::2,::2],gphLats[::2,::2],Uwnd_r[iM,::2,::2],Vwnd_r[iM,::2,::2],
                        latlon=True, pivot='middle', units='inches',scale=Qscale,headwidth=1.5)

        im11=m11.contour(gphLons, gphLats,gph_all[iM,...],cLvls, shading='flat',colors='darkgrey',latlon=True)          
        im11=m11.contour(gphLons, gphLats,gph_r[iM,...],cLvls, shading='flat',colors='k',latlon=True)    

        fig.set_size_inches(30,10);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/partial correlations/'+
                    pcType+'/'+modes[iM]+'_'+ENSOyr+'_'+seasons[iS%12]+'_yrOffset'+str(yrOffset)+notes+'.png')
        plt.close()
            


"""
        
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()          



        im11=m11.pcolormesh(mLons,mLats,maxT_r[iM,...],shading='flat',cmap='PuOr_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)  
        lonStpl,latStpl = m11(mLons[np.isfinite(maxT_stpl[iM,...])],mLats[np.isfinite(maxT_stpl[iM,...])]); 
        m11.plot(lonStpl,latStpl,'ko',markersize=.2)
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");
      
    #    im21=m21.pcolormesh(gphLons, gphLats,omega_r,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)
    #    cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
    #    im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],Uwnd_r[::2,::2],Vwnd_r[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1,scale=5*35, alpha=0.75)
    #    qk = plt.quiverkey(im21, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(30,20)
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/partial correlations/'+pcType+'/'+modes[iM]+'_'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_yrOffset'+str(yrOffset)+notes+'_temp.png')
        plt.close()
    
    
        fig = plt.figure();
        ax11=plt.subplot(211);
        ax21=plt.subplot(212);
        m11 = Basemap(projection='kav7',lon_0=0,ax=ax11); m11.drawcoastlines();m11.drawcountries()
        m21 = Basemap(projection='kav7',lon_0=0,ax=ax21); m21.drawcoastlines();m21.drawcountries()
        im11=m11.pcolormesh(pLons, pLats,prcp_r[iM,...],shading='flat',cmap=prcp_map,latlon=True,vmin=-1.1,vmax=1.1)  
        lonStpl,latStpl = m11(pLons[np.isfinite(prcp_stpl[iM,...])],pLats[np.isfinite(prcp_stpl[iM,...])]); 
        m11.plot(lonStpl,latStpl,'ko',markersize=.3)
        cb11 = m11.colorbar(im11,"bottom", size="5%", pad="2%");             
    #    im21=m21.pcolormesh(gphLons, gphLats,omega_r,shading='flat',cmap='RdYlBu_r',latlon=True,vmin=-1.1,vmax=1.1, alpha=0.75)
    #    cb21 = m21.colorbar(im21,"bottom", size="5%", pad="2%");
    #    im21=m21.quiver(gphLons[::2,::2],gphLats[::2,::2],Uwnd_r[::2,::2],Vwnd_r[::2,::2],latlon=True, units='width',width=0.003,headwidth=.5,headlength=1,scale=5*35, alpha=0.75)
    #    qk = plt.quiverkey(im21, 0.642, 0.5, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure',labelsep=0.05,fontproperties={'size': 20})   
        fig.set_size_inches(30,20);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/partial correlations/'+pcType+'/'+modes[iM]+'_'+str(lvl)+ENSOyr+'_'+seasons[iS%12]+'_yrOffset'+str(yrOffset)+notes+'.png')
        plt.close()
"""