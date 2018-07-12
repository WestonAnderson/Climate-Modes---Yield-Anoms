#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 16:18:00 2018

@author: weston
"""
import time
import numpy as np
import pandas as pd
import netCDF4
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
start = time.clock()
month = 10 #note that this is in python indexing, but harvest dates below aren't
latMax = 5.1; latMin = -5.1; lonMax = 240.1; lonMin = 190.1 #nino 3.4 box
latMaxENSO = 20; latMinENSO = -20; lonMaxENSO = 360; lonMinENSO = 0#define the lat lons used in the first SVD
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']


#Read in ENSO PC data
sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
ensoPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/all crops/sst_all crops_'+
          'Jan-DecENSOyr0har_PCs_1981-2011_janSVD_allTrops.npy')  
pc = -ensoPCs[:,0]+ensoPCs[:,1]

ensoEOFs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/all crops/sst_all crops_'+
          'Jan-DecENSOyr0har_EOFs_1981-2011_janSVD_allTrops.npy') 
ENSOsst = -ensoEOFs[:,0]*np.std(ensoPCs[:,0])+ensoEOFs[:,1]*np.std(ensoPCs[:,1])
ENSOsst = np.squeeze(ENSOsst)
ensoLats=sstData.variables['Y'][(sstData.variables['Y'][:]<=latMaxENSO)&(sstData.variables['Y'][:]>=latMinENSO)]
ensoLons=sstData.variables['X'][(sstData.variables['X'][:]<=lonMaxENSO)&(sstData.variables['X'][:]>=lonMinENSO)]
ENSOsst = ENSOsst.reshape([12,ensoLats.shape[0],ensoLons.shape[0]])
ENSOsst =ENSOsst[:,:,np.where((ensoLons<=lonMax)&(ensoLons>=lonMin))[0]]
ENSOsst =ENSOsst[:,np.where((ensoLats<=latMax)&(ensoLats>=latMin))[0],:]
EOFind = np.nanmean(np.nanmean(ENSOsst,1),1)


#Read in ENSO information for later
enso34tab = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/ENSO/NinoIndex34.csv', 
                     header=1)
ensoAll = np.array(enso34tab.loc[str(1950):str(2013)])
ensoAll = ensoAll[:,0:12]
ensoThresh = -np.std(ensoAll[:,month])/2.
coldEvents = np.ones([1,48])*np.nan
coldYears = np.zeros(1,dtype=int)
lstYr = 0
for row in range(1,ensoAll.shape[0]):

        if ensoAll[row,month]<=ensoThresh: #mark the threshold variables
            if lstYr == row-1: #if last year was a cold year skip this one, since it will be included as the t+1 
                if ensoAll[lstYr,month]<=ensoThresh:
                    lstYr=row;continue
                else: continue 
            event = np.ravel(ensoAll[row-1:row+3,:]) 
            coldEvents = np.append(coldEvents,[event],axis=0) 
            coldYears = np.append(coldYears,[1950+row],axis=0) #cold years are the ENSO years, not the t-1 yrs 
            lstYr = row
coldEvents = coldEvents[1:,:] #strip the empty first event
coldYears = coldYears[1:] #strip the empty first event    
    
ENSOmean = np.mean(coldEvents,axis=0).T


vLines = [4.5,4.5+12,4.5+24,4.5+36]

fig = plt.figure(); 
ax1=plt.subplot(111);
ax1.plot(coldEvents.T, ls='--', color='grey')
ax1.set_xticks(np.arange(0,48,3));
ax1.set_xticklabels([months[x] for x in np.arange(0,48,3)%12])
ax1.set_xlim(0,41);
#ax1.vlines(vLines,-3,3,linewidth=1,linestyle='dotted')
for iy in range(coldYears.size):
    ax1.text(iy*2.75+3,2.6,str(str(coldYears[iy])+','),
         horizontalalignment='right',color='grey')
ax1.fill_between(np.arange(0,48,1),ENSOmean,where=ENSOmean<=0,color='cornflowerblue', alpha=.75,)
ax1.fill_between(np.arange(0,48,1),ENSOmean,where=ENSOmean>=-.07,color='firebrick', alpha=.75,)
ax1.plot(range(12,24),EOFind*.5,'--k',linewidth=3)
ax1.set_ylabel('Oceanic Nino Index', size=12)
ax1.set_ylim(-2.5,2.5)
plt.show()