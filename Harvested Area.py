#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 11:10:12 2018

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
import matplotlib as mpl
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #
thresh = 0.0001 #0.005 = 0.5%
#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#Correct the names to match the mapping names 
cWght_f = netCDF4.Dataset('/Volumes/Data_Archive/Data/M3_NCDF/5ArcMin/HA_Yld/wheat_HarvAreaYield_NetCDF/wheat_AreaYieldProduction.nc')
mMsk = np.load(str('/Volumes/Data_Archive/Data/CropMaps_ComonGrid/AllArea/HarvestedArea/5arcMin_ensemble/ensemble_m.npy'))
wMsk = np.load(str('/Volumes/Data_Archive/Data/CropMaps_ComonGrid/AllArea/HarvestedArea/5arcMin_ensemble/ensemble_w.npy'))
sMsk = np.load(str('/Volumes/Data_Archive/Data/CropMaps_ComonGrid/AllArea/HarvestedArea/5arcMin_ensemble/ensemble_s.npy'))
lats = cWght_f['latitude'][:]
lons = cWght_f['longitude'][:] 
cLats = np.append((lats[0]-lats[1])+lats[0], lats)
cLats = cLats +(lats[1]-lats[0])/2
dlon = (lons[1]-lons[0])*(np.pi/180.) #difference in longitude (in radians) remains constant
R = 6371000. #Radius of the earth in m
cellArea = 2*np.pi*(R**2)*np.abs(np.sin(cLats[1:]*(np.pi/180.))-
                           np.sin(cLats[:-1]*(np.pi/180.)))*dlon #in m
cellArea = np.repeat(cellArea[:,np.newaxis],lons.shape,axis=1)
haThresh = (cellArea*thresh)*1./10000 #converting the threshhold to hectares
cMsk = mMsk+sMsk+wMsk
ha = cMsk.copy()
cMsk_r = cMsk.copy()
cMsk[cMsk>=haThresh]=np.nan
cMsk[cMsk<haThresh]=1;
cMsk_r[cMsk_r<haThresh]=np.nan
cMsk_r[cMsk_r>=haThresh]=1
lons,lats=np.meshgrid(lons,lats)
cMsk = maskoceans(lons, lats, cMsk)
ha = maskoceans(lons, lats, ha)
cMsk_r = maskoceans(lons, lats, cMsk_r)     

thisCMAP = 'BrBG'
clrBr = np.arange(-0.5,0.5+.01,.01)
norm = Normalize(vmin=-0.5, vmax=0.5, clip=False)
mapper=cm.ScalarMappable(norm=norm, cmap=thisCMAP)


fig = plt.figure();ax1=plt.subplot(111)
##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##
# Plot the observations first
##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##
plt.sca(ax1)
m = Basemap(projection='cyl',llcrnrlat=-45,urcrnrlat=60,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c',ax=ax1)
m1 = m.pcolor(lons,lats,ha,cmap='YlGnBu',zorder=1,latlon=True) 
m.pcolor(lons,lats,cMsk,cmap='tab20c_r',zorder=2,latlon=True) 
#Plot everything at the country level first from FAO data before going into subnational data
world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                        name='cnts', drawbounds=True)
USshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/USA_adm_shp/USA_adm1', 
            name='states', drawbounds=True)
BRshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/BRA_adm_shp/BRA_adm1', 
                name='states', drawbounds=True)
ARshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/ARG_adm_shp/ARG_adm1', 
                name='states', drawbounds=True)
CNshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CHN_adm_shp/CHN_adm1', 
                name='states', drawbounds=True)
AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
            name='states', drawbounds=True)
CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                name='states', drawbounds=True)
INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
            name='states', drawbounds=True)
MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                name='states', drawbounds=True)
cbar = m.colorbar(m1, location='right')

fig.set_size_inches(50,15)
plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/Writeup/wheat_maize_soy_HA.png')


