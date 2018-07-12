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

years = [1981,2011]
PCs = [1] #not python indexing
monCutoff = 1 #not python indexing
lvl = 200 # geopotential height lvl
sigLvl = .1682 #.995, .1682     #vertical velocity
timePeriod =''# ' 1960-2010' #'' or ' 1960-2010'
notes = '_janSVD'
nao_crops = 'wheat' #'all crops' #
NAOseason = 'DJF' #DJF, JJA, MAM etc
Tvar = 'TMAX' #'TMIN' or 'TMAX'
NAOind = 'NAO' #NAO, NAOg pr EA_WR or WP
iS=0 #season number to use for climate variables
#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

## colormap function
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/custom_div_cmap.py').read())
prcp_map = custom_div_cmap(20, mincol='saddlebrown', midcol='white' ,maxcol='CornflowerBlue')
##

seasons = ['DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ']
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr0har'
numYrs = years[1]-years[0]
levels = np.arange(-.75,.85,.1)

#~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~#
#read in the netCDF file containing the variables

#read in the netCDF file containing the variables
ncGeo = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/gph.nc',
                        'r',format='NETCDF4') #read in the netCDF file for the variable
geoVar=ncGeo.variables['hgt'][(ncGeo.variables['year'][:]>=years[0]-1)&(ncGeo.variables['year'][:]<years[1]+2),
                              :,np.where(ncGeo.variables['level'][:]==lvl)[0][0],:,:]

ncSM = netCDF4.Dataset('/Volumes/Data_Archive/Data/SoilMoisture/GLDAS/Noah10_nc/detrended/SoilM.nc','r',format='NETCDF4') 
sm = np.nansum(ncSM.variables['SoilM'][(ncSM.variables['year'][:]>=years[0]-1)&(ncSM.variables['year'][:]<years[1]+2),:,:3,:,:],2)

ncMaxT = netCDF4.Dataset('//Volumes/Data_Archive/Data/Temp/BerkeleyEarth/monthly/Complete_'+Tvar+'_LatLong1.nc','r',format='NETCDF4') #read in the netCDF file for the variable
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

cPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/crop_'+NAOseason+nao_crops+
              '_'+ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy') 
sstPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/gph_'+NAOseason+nao_crops+
              '_'+ENSOyr+'_PCs_'+str(years[0])+'-'+str(years[1])+notes+'.npy')             


#~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~##~#~#~#~#~#~#~~#~#~#~#~#~#~#           
for PC in PCs:
    cPC = cPCs[:,PC-1]
    sstPC = sstPCs[:,PC-1]
    #Regressions by season
    yrOffset = iS//12

    #Read in the NAO time series to plot against
    if NAOind is 'NAO':
        NAOdf = pd.read_csv('/Volumes/Data_Archive/Data/NAO/nao_station_seasonal.csv',index_col=0)
        NAO = NAOdf.loc[years[0]:years[1]][seasons[iS%12]].values
    elif NAOind is 'NAOg':
        NAOdf = pd.read_csv('/Volumes/Data_Archive/Data/NAO/nao_g.csv',index_col=0)
        NAO = NAOdf.loc[years[0]+1:years[1]+1][seasons[iS%12]].values
    elif NAOind is 'EA_WR':
        NAOdf = pd.read_csv('/Volumes/Data_Archive/Data/EA_WR/ea.data.csv',index_col=0)
        NAO = NAOdf.loc[years[0]+1:years[1]+1][seasons[iS%12]].values
    elif NAOind is 'WP':
        NAOdf = pd.read_csv('/Volumes/Data_Archive/Data/WP/wp.data.csv',index_col=0)
        NAO = NAOdf.loc[years[0]+1:years[1]+1][seasons[iS%12]].values
    elif NAOind is 'EA':
        NAOdf = pd.read_csv('/Volumes/Data_Archive/Data/EA/ea.csv',index_col=0)
        NAO = NAOdf.loc[years[0]+1:years[1]+1][seasons[iS%12]].values
    elif NAOind is 'SCAND':
        NAOdf = pd.read_csv('/Volumes/Data_Archive/Data/SCAND/scand.csv',index_col=0)
        NAO = NAOdf.loc[years[0]+1:years[1]+1][seasons[iS%12]].values
    else:
        print('ERROR: enter a valid NAO index choice')
        break
    
    fig = plt.figure();ax=plt.subplot(111)
    ax.plot(range(years[0],years[1]+1),NAO/np.std(NAO),'b',label='NAO')
    ax.plot(range(years[0],years[1]+1),sstPC/np.std(sstPC),'k',label='GPH PC'+str(PC))
    ax.plot(range(years[0],years[1]+1),cPC/np.std(cPC),'--k',label='crop yield PC'+str(PC))
    ax.text(2005,2.7,NAOind+'-crop PC: '+str(np.round(np.corrcoef(cPC,NAO)[0][1],4)))
    ax.text(2005,2.9,NAOind+'-sst PC: '+str(np.round(np.corrcoef(sstPC,NAO)[0][1],4)))
    ax.text(2005,3.1,'sst PC-crop PC: '+str(np.round(np.corrcoef(cPC,sstPC)[0][1],4)))
    ax.legend()
    ax.set_ylim(-2.5,2.5);
    fig.set_size_inches(8,5);
    fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/NAO/SST/NAO Index/'+nao_crops+'_gphPC'+str(PC)+ENSOyr+'_'+seasons[iS%12]+'_yrOffset'+str(yrOffset)+timePeriod+notes+NAOind+'.png')
    plt.close()
