#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 08:22:35 2018

@author: weston
"""

import numpy as np
import pandas as pd
import time
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
start = time.clock()

ensoMons =[5,6,7,8,9,10,11,12,13,14,15,16]
ENyrs = np.array([1982, 1986, 1991, 1994, 1997, 2002, 2004, 2006, 2009])  #EN,  # 2009 not possible with SWC included
LNyrs = np.array([1983, 1988, 1995, 1998, 2007, 2010])

yrMin = 1981
yrMax = 2011
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
latMax = 5; latMin = -5; lonMax = 270; lonMin = 160

#Read in the SST anomaly data
sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
#T is months since 1960-01-01, and I want to start at 1950 with one year lag
sstLats=sstData.variables['Y'][(sstData.variables['Y'][:]<latMax)&(sstData.variables['Y'][:]>latMin)]
sstLons=sstData.variables['X'][(sstData.variables['X'][:]<lonMax)&(sstData.variables['X'][:]>lonMin)]
sstLons,sstLats=np.meshgrid(sstLons,sstLats)

sstAnom=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-int(yrMin)))&
                        (sstData.variables['T'][:]<=12*((2+int(yrMax))-1960)),:,
                        (sstData.variables['Y'][:]<latMax)&(sstData.variables['Y'][:]>latMin),
                        (sstData.variables['X'][:]<lonMax)&(sstData.variables['X'][:]>lonMin)]
#pull only chosen month SST anomalies
sstAnom = np.squeeze(sstAnom);sstAnom=np.ma.filled(sstAnom,0)

ind = np.repeat(np.array(range(1+int(yrMax)-int(yrMin)))*12,np.shape(ensoMons)[0])+\
                np.tile(ensoMons,1+int(yrMax)-int(yrMin))

sstAnom=sstAnom[ind,...]
sstAnom = np.reshape(sstAnom,[1+int(yrMax)-int(yrMin),np.shape(ensoMons)[0],sstAnom.shape[1],sstAnom.shape[2]])

sstAnomEN = sstAnom[ENyrs-yrMin,...]
sstAnomEN = np.nanmean(sstAnomEN,0)


y_ticks = np.arange(0,sstAnomEN.shape[0]*sstAnomEN.shape[1],sstAnomEN.shape[1])+sstAnomEN.shape[0]/2-4
y_labs = [months[m%12] for m in ensoMons]
sstAnomEN = np.reshape(sstAnomEN,[sstAnomEN.shape[0]*sstAnomEN.shape[1],sstAnomEN.shape[2]])
colMax = np.max([np.abs(np.nanmin(sstAnomEN)),np.nanmax(sstAnomEN)])
plt.figure();plt.imshow(sstAnomEN,cmap='RdBu_r',vmin=-colMax,vmax=colMax);plt.colorbar()
plt.yticks(y_ticks,y_labs)
plt.show()


sstAnomLN = sstAnom[LNyrs-yrMin,...]
sstAnomLN = np.nanmean(sstAnomLN,0)


y_ticks = np.arange(0,sstAnomLN.shape[0]*sstAnomLN.shape[1],sstAnomLN.shape[1])+sstAnomLN.shape[0]/2-4
y_labs = [months[m%12] for m in ensoMons]
x_ticks = range()
sstAnomLN = np.reshape(sstAnomLN,[sstAnomLN.shape[0]*sstAnomLN.shape[1],sstAnomLN.shape[2]])
colMax = np.max([np.abs(np.nanmin(sstAnomLN)),np.nanmax(sstAnomLN)])
plt.figure();plt.imshow(sstAnomLN,cmap='RdBu_r',vmin=-colMax,vmax=colMax);plt.colorbar()
plt.yticks(y_ticks,y_labs)
plt.show()
    