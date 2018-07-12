#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 14:09:52 2018

@author: weston
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.basemap import maskoceans
from unidecode import unidecode
import pandas as pd
from scipy import signal
import netCDF4
import matplotlib
plt.ioff()

notes = '_janSVD_allTrops'
latMax = 20.01; latMin = -20.01; lonMax = 360; lonMin = 0

enso = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+
        '/all crops/sst_all crops_Jan-DecENSOyr0har_EOFs_1981-2011_janSVD_allTrops.npy')

pcs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+
        '/all crops/sst_all crops_Jan-DecENSOyr0har_PCs_1981-2011_janSVD_allTrops.npy')


#Read in the SST anomaly data
sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
#T is months since 1960-01-01, and I want to start at 1950 with one year lag
sstLats=sstData.variables['Y'][(sstData.variables['Y'][:]<latMax)&(sstData.variables['Y'][:]>=latMin)]
sstLons=sstData.variables['X'][(sstData.variables['X'][:]<lonMax)&(sstData.variables['X'][:]>=lonMin)]
sstLons,sstLats=np.meshgrid(sstLons,sstLats)

enso = np.reshape(enso, [12,sstLons.shape[0],sstLats.shape[1],624])
enso = (-enso[:,:,:,0]*np.std(pcs[:,0]))*.5+(enso[:,:,:,1]*np.std(pcs[:,1]))*.5
enso[(enso<0.01)&(enso>-0.01)]=np.nan

sstLvls = np.arange(-1.4,1.6,.2)

fig = plt.figure();
ax1=plt.subplot(411);ax2=plt.subplot(412)
ax3=plt.subplot(413);ax4=plt.subplot(414)
ax1.set_ylabel('1 JFM           ',rotation=0,fontweight='bold',color='firebrick')
ax2.set_ylabel('2 AMJ           ',rotation=0,fontweight='bold',color='firebrick')
ax3.set_ylabel('3 JAS           ',rotation=0,fontweight='bold',color='cornflowerblue')
ax4.set_ylabel('4 OND           ',rotation=0,fontweight='bold',color='cornflowerblue')
m1 = Basemap(projection='cyl',llcrnrlat=-20,urcrnrlat=20,\
            llcrnrlon=0,urcrnrlon=360,lon_0=180,resolution='c',ax=ax1)
m1.drawcoastlines()
m1 = ax1.contourf(sstLons,sstLats,np.nanmean(enso[:3,...],0),cmap='RdBu_r',
                latlon=True,levels=sstLvls,extend='both')
m2 = Basemap(projection='cyl',llcrnrlat=-20,urcrnrlat=20,\
            llcrnrlon=0,urcrnrlon=360,lon_0=180,resolution='c',ax=ax2)
m2.drawcoastlines()
m2 = ax2.contourf(sstLons,sstLats,np.nanmean(enso[3:6,...],0),cmap='RdBu_r',
                latlon=True,levels=sstLvls,extend='both')
m3 = Basemap(projection='cyl',llcrnrlat=-20,urcrnrlat=20,\
            llcrnrlon=0,urcrnrlon=360,lon_0=180,resolution='c',ax=ax3)
m3.drawcoastlines()
m3 = ax3.contourf(sstLons,sstLats,np.nanmean(enso[6:9,...],0),cmap='RdBu_r',
                latlon=True,levels=sstLvls,extend='both')
m4 = Basemap(projection='cyl',llcrnrlat=-20,urcrnrlat=20,\
            llcrnrlon=0,urcrnrlon=360,lon_0=180,resolution='c',ax=ax4)
m4.drawcoastlines()
im4 = ax4.contourf(sstLons,sstLats,np.nanmean(enso[9:,...],0),cmap='RdBu_r',
                latlon=True,levels=sstLvls,extend='both')
fig.set_size_inches(8,4);
plt.show()

"""
fig = plt.figure();
ax1=plt.subplot(111)
m4 = Basemap(projection='cyl',llcrnrlat=-20,urcrnrlat=20,\
            llcrnrlon=200,urcrnrlon=280,lon_0=180,resolution='c',ax=ax1)
m4.drawcoastlines()
im4 = ax1.contourf(sstLons,sstLats,np.nanmean(enso[9:,...],0),cmap='RdBu_r',
                latlon=True,levels=sstLvls,extend='both')
cb4=m4.colorbar(im4,"bottom",pad=0.1); 

plt.show()
"""
