#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 17:55:02 2018

@author: weston
"""
from eofs.standard import Eof
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from mpl_toolkits.basemap import Basemap
import matplotlib

yrMin = 1981
yrMax = 2011
ensoMons =[2,3,4,5,6,7,8,9,10,11,12,13]
PCs = 3 #python indexing

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                 Helper Functions                  #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#Read in moving average: moving_average(a, smooth)
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/movingAverage.py').read())
#Read in moving standard deviation: moving_std(a, smooth)
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/movingStd.py').read())
#create surrogate timeseries
#follows from Ebusuzaki (1997) and Schrieber and Shmitz (2000)
#  function call is   ############ surrogate(input TS): ############                  #
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/surrogateData.py').read())
#Read in the shapefiles
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/GADM_dict.py').read())

# Restrict SSTs to tropical Pacific and Indian Ocean
latMax = 5; latMin = -5; lonMax = 270; lonMin = 160
#latMax = 15; latMin = -15; lonMax = 360; lonMin = 0
#latMax = 5; latMin = -5; lonMax = 240; lonMin = 190
#latMax = 15; latMin = -15; lonMax = 360; lonMin = 0

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

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

sstAnom=sstAnom[np.repeat(np.array(range(1+int(yrMax)-int(yrMin)))*12,np.shape(ensoMons)[0])+
                np.tile(ensoMons,1+int(yrMax)-int(yrMin)),...]
sstAnom = np.reshape(sstAnom,[1+int(yrMax)-int(yrMin),np.shape(ensoMons)[0],sstAnom.shape[1],sstAnom.shape[2]])
#sstAnom = np.nanmean(sstAnom,1)

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   Calculate the EOF                       #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#             

#reshape each dataset
sstAnomRav=np.reshape(sstAnom, [sstAnom.shape[0],sstAnom.shape[1]*sstAnom.shape[2]*sstAnom.shape[3]],'C')

# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslatSST = np.cos(np.deg2rad(sstLats))
wgts = np.sqrt(coslatSST)
# add on a series of ones to indicate that the production anomalies in the
# EOF should not be production weighted
wghtsRav = np.repeat(np.ravel(wgts,'C'),sstAnom.shape[1]) 

solver = Eof(sstAnomRav, weights=wghtsRav)

eofs = solver.eofsAsCorrelation(neofs=3)
pcs = solver.pcs(npcs=3, pcscaling=1)

eofs = np.reshape(eofs,sstAnom[:3,:].shape)


for k in range(PCs):
    eof_sst = np.squeeze(eofs[k,...])
    eof_sst = eof_sst[:,::-1,:]
    y_ticks = np.arange(0,eof_sst.shape[0]*eof_sst.shape[1],eof_sst.shape[1])+eof_sst.shape[0]/2
    y_labs = [months[m%12] for m in ensoMons]
    eof_sst = np.reshape(eof_sst,[eof_sst.shape[0]*eof_sst.shape[1],eof_sst.shape[2]])
    colMax = np.max([np.abs(np.nanmin(eof_sst)),np.nanmax(eof_sst)])
    plt.figure();plt.imshow(eof_sst,cmap='RdBu_r',vmin=-colMax,vmax=colMax);plt.colorbar()
    plt.yticks(y_ticks,y_labs)