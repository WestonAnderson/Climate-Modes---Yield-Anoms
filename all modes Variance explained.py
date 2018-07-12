#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 13:55:33 2017

@author: weston
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from unidecode import unidecode
import pandas as pd
import netCDF4
import numpy as np
import matplotlib
from mpl_toolkits.basemap import maskoceans
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #

SVDmon = '_janSVD_allTrops'
NAOmons = 'DJF'
mod = 'SVD'
monCutoff = 1 #not python indexing
thresh = 0.0001 #0.0001 = 0.01%
years=[1981,2010]
modes = ['ENSO+IOD+TAV+NAO','ENSO']
varType = 'linAlg' #'linAlg' or 'comp'

#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr0har'

exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/GADM_dict.py').read())


#Define a colorbar
clrMax = .5;clrMin=0

thisCMAP =  'YlOrBr'# matplotlib.colors.Colormap(Margot2_4.mpl_colormap)
clrBr = np.arange(clrMin,clrMax+.01,.01)
norm = Normalize(vmin=clrMin, vmax=clrMax, clip=False)
mapper=cm.ScalarMappable(norm=norm, cmap=thisCMAP)


#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   Read in the data                #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*# 
# First read in each SVD reconstruction for each mode of variability
# Then calculate the variance explained for ENSO, ENSO+IOD, ENSO+IOD+TAV and ENSO+IOD+TAV+NAO

   
## create an empty variance explained DF ## 
ensoVar = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/'+
                           'Climate Modes - Yield Anoms/stats regressions/all cropsVar1961_2012_'+ENSOyr+'_SVD2'+SVDmon+'.csv')   

ensoVar['index'] = [x.split('_')[0]+'_'+unidecode(x.split('_')[1].upper()) for x in ensoVar['index']]
varExp = pd.DataFrame(data={'location':ensoVar['index']})
varExp.set_index('location',inplace=True,drop=True)
varExp['var'] = np.nan #variance
varExp['ENSO_abs'] = np.nan #variance of one mode
varExp['clim_abs'] = np.nan
varExp['IOD_abs'] = np.nan 
varExp['TAV_abs'] = np.nan 
varExp['NAO_abs'] = np.nan 

varExp['ENSO'] = np.nan #variance of one mode
varExp['IOD'] = np.nan 
varExp['TAV'] = np.nan 
varExp['NAO'] = np.nan 
varExp['ENSO+IOD'] = np.nan #variance of combined modes
varExp['ENSO+NAO'] = np.nan 
varExp['ENSO+IOD+TAV'] = np.nan
varExp['ENSO+IOD+TAV+NAO'] = np.nan
varExp['ENSO_marg'] = np.nan #marginal variance
varExp['IOD_marg'] = np.nan 
varExp['TAV_marg'] = np.nan 
varExp['NAO_marg'] = np.nan 
states = ensoVar['index']
del ensoVar

## IOD ## 
iod = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/SVD/'+
                             'all crops1961_2012_'+ENSOyr+'_iodSVD1'+SVDmon+'.csv')  ## TAV ## 
iod.reset_index(inplace=True);
iod.set_index(['state','year'],inplace=True,drop=True)
iod.rename(columns={'SVD':'IOD'},inplace=True)
## TAV ## 
tav = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/SVD/'+
                             'all crops1961_2012_'+ENSOyr+'_tavSVD2'+SVDmon+'.csv')  
tav.reset_index(inplace=True);
tav.set_index(['state','year'],inplace=True,drop=True)
tav.rename(columns={'SVD':'TAV'},inplace=True)
###### ENSO ######
enso = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/SVD/'+
                             'all crops1961_2012_'+ENSOyr+'_SVD2'+SVDmon+'.csv')   
enso.reset_index(inplace=True);
enso.set_index(['state','year'],inplace=True,drop=True)
enso.rename(columns={'SVD':'ENSO'},inplace=True)
###### NAO ######
nao = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/SVD/'+
                             'DJFwheat1961_2012_'+ENSOyr+'_naoSVD1_janSVD.csv')   
nao.reset_index(inplace=True);
nao.set_index(['state','year'],inplace=True,drop=True)
nao.rename(columns={'SVD':'NAO'},inplace=True)

yldAnoms = pd.concat([enso,nao['NAO'],iod['IOD'],tav['TAV']],axis=1)
yldAnoms.reset_index(inplace=True);
yldAnoms = yldAnoms[(yldAnoms.year>=years[0])&(yldAnoms.year<=years[1])]
yldAnoms['state'] = [x.split('_')[0]+'_'+unidecode(x.split('_')[1].upper()) for x in yldAnoms['state']]
yldAnoms.set_index(['state','year'],inplace=True,drop=True)

for i in range(0,np.size(states)):
    varExp['var'].loc[states[i]] = np.var(yldAnoms['yldAnomGau'][states[i]])
    varExp['ENSO_abs'].loc[states[i]]=np.var(yldAnoms['ENSO'][states[i]])
    varExp['clim_abs'].loc[states[i]]=np.var( yldAnoms[['ENSO','IOD','TAV','NAO']].loc[states[i]].sum(axis=1))
    
    varExp['ENSO'].loc[states[i]]=1-\
                        (np.var(yldAnoms['yldAnomGau'][states[i]]-yldAnoms['ENSO'][states[i]])/
                        np.var(yldAnoms['yldAnomGau'][states[i]]))

    varExp['IOD'].loc[states[i]]=1-\
                        (np.var(yldAnoms['yldAnomGau'][states[i]]-yldAnoms['IOD'][states[i]])/
                        np.var(yldAnoms['yldAnomGau'][states[i]]))
    varExp['ENSO+NAO'].loc[states[i]]=1-\
                        (np.var(yldAnoms['yldAnomGau'][states[i]]-( yldAnoms[['ENSO','NAO']].loc[states[i]].sum(axis=1)))/
                        np.var(yldAnoms['yldAnomGau'][states[i]]))
    varExp['ENSO+IOD'].loc[states[i]]=1-\
                        (np.var(yldAnoms['yldAnomGau'][states[i]]-( yldAnoms[['ENSO','IOD']].loc[states[i]].sum(axis=1)))/
                        np.var(yldAnoms['yldAnomGau'][states[i]]))
    varExp['TAV'].loc[states[i]]=1-\
                        (np.var(yldAnoms['yldAnomGau'][states[i]]-yldAnoms['TAV'][states[i]])/
                        np.var(yldAnoms['yldAnomGau'][states[i]]))
    varExp['ENSO+IOD+TAV'].loc[states[i]]=1-\
                        (np.var(yldAnoms['yldAnomGau'][states[i]]-( yldAnoms[['ENSO','IOD','TAV']].loc[states[i]].sum(axis=1)))/
                        np.var(yldAnoms['yldAnomGau'][states[i]]))
    varExp['NAO'].loc[states[i]]=1-\
                        (np.var(yldAnoms['yldAnomGau'][states[i]]-yldAnoms['NAO'][states[i]])/
                        np.var(yldAnoms['yldAnomGau'][states[i]]))
    varExp['ENSO+IOD+TAV+NAO'].loc[states[i]]=1-\
                        (np.var(yldAnoms['yldAnomGau'][states[i]]-( yldAnoms[['ENSO','IOD','TAV','NAO']].loc[states[i]].sum(axis=1)))/
                        np.var(yldAnoms['yldAnomGau'][states[i]]))

    varExp['NAO_marg'].loc[states[i]]=varExp['NAO'].fillna(0)[states[i]]
    varExp['ENSO_marg'].loc[states[i]]=varExp['ENSO+NAO'][states[i]]-varExp['NAO'].fillna(0)[states[i]]
    varExp['IOD_marg'].loc[states[i]]=varExp['ENSO+IOD'][states[i]]-varExp['ENSO'].fillna(0)[states[i]]
    varExp['TAV_marg'].loc[states[i]]=varExp['ENSO+IOD+TAV'][states[i]]-varExp['ENSO+IOD'].fillna(0)[states[i]]

for i in range(0,np.size(states)):
    varExp['var'].loc[states[i]] = (np.linalg.norm(yldAnoms['yldAnomGau'][states[i]])**2)
    varExp['ENSO_abs'].loc[states[i]]=np.var(yldAnoms['ENSO'][states[i]])
    varExp['clim_abs'].loc[states[i]]=np.var( yldAnoms[['ENSO','IOD','TAV','NAO']].loc[states[i]].sum(axis=1))
    
    varExp['ENSO'].loc[states[i]]=((np.linalg.norm(yldAnoms['ENSO'][states[i]])**2)/
                                    (np.linalg.norm(yldAnoms['yldAnomGau'][states[i]])**2))
    varExp['IOD'].loc[states[i]]=((np.linalg.norm(yldAnoms['IOD'][states[i]])**2)/
                                    (np.linalg.norm(yldAnoms['yldAnomGau'][states[i]])**2))
    varExp['ENSO+NAO'].loc[states[i]]=((np.linalg.norm(yldAnoms[['ENSO','NAO']].loc[states[i]].sum(axis=1))**2)/
                                    (np.linalg.norm(yldAnoms['yldAnomGau'][states[i]])**2))
    varExp['ENSO+IOD'].loc[states[i]]=((np.linalg.norm(yldAnoms[['ENSO','IOD']].loc[states[i]].sum(axis=1))**2)/
                                    (np.linalg.norm(yldAnoms['yldAnomGau'][states[i]])**2))
    varExp['TAV'].loc[states[i]]=((np.linalg.norm(yldAnoms['TAV'][states[i]])**2)/
                                    (np.linalg.norm(yldAnoms['yldAnomGau'][states[i]])**2))
    varExp['ENSO+IOD+TAV'].loc[states[i]]=((np.linalg.norm(yldAnoms[['ENSO','IOD','TAV']].loc[states[i]].sum(axis=1))**2)/
                                    (np.linalg.norm(yldAnoms['yldAnomGau'][states[i]])**2))
    varExp['NAO'].loc[states[i]]=((np.linalg.norm(yldAnoms['NAO'][states[i]])**2)/
                                    (np.linalg.norm(yldAnoms['yldAnomGau'][states[i]])**2))
    varExp['ENSO+IOD+TAV+NAO'].loc[states[i]]=((np.linalg.norm(yldAnoms[['ENSO','IOD','TAV','NAO']].loc[states[i]].sum(axis=1))**2)/
                                    (np.linalg.norm(yldAnoms['yldAnomGau'][states[i]])**2))
                            
for crop in ['maize','soy','wheat']:        
        for mode in modes:
            cWght_f = netCDF4.Dataset('/Volumes/Data_Archive/Data/M3_NCDF/5ArcMin/HA_Yld/wheat_HarvAreaYield_NetCDF/wheat_AreaYieldProduction.nc')
            cMsk = np.load(str('/Volumes/Data_Archive/Data/CropMaps_ComonGrid/AllArea/HarvestedArea/5arcMin_ensemble/ensemble_'+crop[0].lower()+'.npy'))
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
            cMsk_r = cMsk.copy()
            cMsk[cMsk>=haThresh]=np.nan
            cMsk[cMsk<haThresh]=1;
            cMsk_r[cMsk_r<haThresh]=np.nan
            cMsk_r[cMsk_r>=haThresh]=1
            lons,lats=np.meshgrid(lons,lats)
            cMsk_all = cMsk.copy()
            cMsk = maskoceans(lons, lats, cMsk)
            cMsk_r = maskoceans(lons, lats, cMsk_r)
            
            #plot    
            fig = plt.figure()
            m = Basemap(projection='kav7',lon_0=0, resolution='l')
            
            #Plot everything at the country level first from FAO data before going into subnational data
            world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                    name='cnts', drawbounds=True)
            m.pcolor(lons,lats,cMsk_r,cmap='Dark2_r',zorder=1.2,latlon=True) 
            names = [];colors = {};i=0
            for shape_dict in m.cnts_info:
                names.append(shape_dict['NAME'])
                if ((shape_dict['NAME']=='Armenia')|(shape_dict['NAME']=='Azerbaijan')|(shape_dict['NAME']=='Belarus')|
                (shape_dict['NAME']=='Estonia')|(shape_dict['NAME']=='Kazakhstan')|(shape_dict['NAME']=='Kyrgyzstan')|
                (shape_dict['NAME']=='Latvia')|(shape_dict['NAME']=='Lithuania')|(shape_dict['NAME']=='Moldova')|
                (shape_dict['NAME']=='Russia')|(shape_dict['NAME']=='Tajikistan')|(shape_dict['NAME']=='Ukraine')|
                (shape_dict['NAME']=='Uzbekistan')):
                        colors[shape_dict['NAME']]=varExp.loc[unidecode(str(crop+'_USSR'))][mode]
                elif ((shape_dict['NAME']=='Belgium')|(shape_dict['NAME']=='Luxembourg')):
                    if crop=='soy':continue
                    else:colors[shape_dict['NAME']]=varExp.loc[unidecode(str(crop+'_BELGIUM-LUXEMBOURG'))][mode]
                elif ((shape_dict['NAME']=='Croatia')|(shape_dict['NAME']=='Bosnia and Herz.')|(shape_dict['NAME']=='Slovenia')
                    |(shape_dict['NAME']=='Macedonia')|(shape_dict['NAME']=='Serbia')):
                        colors[shape_dict['NAME']]=varExp.loc[unidecode(str(crop+'_YUGOSLAVIA'))][mode]
                elif len(np.where(np.array(states)==str(crop+'_'+shape_dict['NAME'].upper()))[0]) != 0:
                    state_color = varExp.loc[unidecode(str(crop+'_'+shape_dict['NAME'].upper()))][mode]
                    colors[shape_dict['NAME']]=state_color
                #else: colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.cnts):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]),edgecolor='k',zorder=1.3)
                ax.add_patch(poly)
            mapper.set_array(clrBr)    
        
  
            # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
            USshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/USA_adm_shp/USA_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if (shape_dict['HASC_1']=='US.AK')|(shape_dict['HASC_1']=='US.DC')|(shape_dict['HASC_1']=='US.DE')|(shape_dict['HASC_1']=='US.HI'):
                    continue#colors[shape_dict['HASC_1']]=0
                elif len(np.where(np.array(states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = varExp.loc[str(crop+'_'+shape_dict['NAME_1'].upper())][mode]
                    colors[shape_dict['HASC_1']]=state_color
                #else: colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                ax.add_patch(poly)
            mapper.set_array(clrBr)
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                ax.add_patch(poly)
               
            BRshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/BRA_adm_shp/BRA_adm1', 
                            name='states', drawbounds=True)
            # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(states)==str(crop+'_'+shape_dict['HASC_1'][3:5]))[0]) != 0:
                    state_color = varExp.loc[str(crop+'_'+shape_dict['HASC_1'][3:5])][mode]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                ax.add_patch(poly) 
            mapper.set_array(clrBr)
            
            ARshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/ARG_adm_shp/ARG_adm1', 
                            name='states', drawbounds=True)
            # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(states)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
                    state_color = varExp.loc[unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))][mode]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                ax.add_patch(poly) 
            mapper.set_array(clrBr)
            
            CNshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CHN_adm_shp/CHN_adm1', 
                            name='states', drawbounds=True)
            # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if ((shape_dict['NAME_1']=='Chongqing')):
                        colors[shape_dict['HASC_1']]=varExp.loc[unidecode(str(crop+'_SICHUAN'))][mode]
                elif len(np.where(np.array(states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = varExp.loc[str(crop+'_'+shape_dict['NAME_1'].upper())][mode]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                ax.add_patch(poly) 
            mapper.set_array(clrBr)
            
            if crop is 'wheat':    
                CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                                name='states', drawbounds=True)
                # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
                names = [];colors = {};i=0
                for shape_dict in m.states_info:
                    names.append(shape_dict['HASC_1'])
                    if len(np.where(np.array(states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                        state_color = varExp.loc[str(crop+'_'+shape_dict['NAME_1'].upper())][mode]
                        colors[shape_dict['HASC_1']]=state_color
                    else: continue#colors[shape_dict['HASC_1']]=0
                ax = plt.gca() # get current axes instance
                for nshape, seg in enumerate(m.states):
                    if names[nshape] not in colors: continue
                    poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                    ax.add_patch(poly)
                
                AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                            name='states', drawbounds=True)
                names = [];colors = {};i=0
                for shape_dict in m.states_info:
                    names.append(shape_dict['HASC_1'])
                    if len(np.where(np.array(states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                        state_color = varExp.loc[str(crop+'_'+shape_dict['NAME_1'].upper())][mode]
                        colors[shape_dict['HASC_1']]=state_color
                    else: continue#colors[shape_dict['HASC_1']]=0
                ax = plt.gca() # get current axes instance
                for nshape, seg in enumerate(m.states):
                    if names[nshape] not in colors: continue
                    poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                    ax.add_patch(poly)
            
                INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                            name='states', drawbounds=True)
                names = [];colors = {};i=0
                for shape_dict in m.states_info:
                    names.append(shape_dict['HASC_1'])
                    if ((shape_dict['NAME_1']=='Telangana')):
                        colors[shape_dict['HASC_1']]=varExp.loc[unidecode(str(crop+'_Andhra Pradesh'.upper()))][mode]
                    elif ((shape_dict['NAME_1']=='Chhattisgarh')):
                        colors[shape_dict['HASC_1']]=varExp.loc[unidecode(str(crop+'_Madhya Pradesh'.upper()))][mode]
                    elif ((shape_dict['NAME_1']=='Jharkhand')):
                        colors[shape_dict['HASC_1']]=varExp.loc[unidecode(str(crop+'_Bihar'.upper()))][mode]
                    elif len(np.where(np.array(states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                        state_color = varExp.loc[unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))][mode]
                        colors[shape_dict['HASC_1']]=state_color
                    else: continue#colors[shape_dict['HASC_1']]=0
                ax = plt.gca() # get current axes instance
                for nshape, seg in enumerate(m.states):
                    if names[nshape] not in colors: continue
                    poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                    ax.add_patch(poly)
                mapper.set_array(clrBr)
                
            if crop is 'maize':
                MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                            name='states', drawbounds=True)
                names = [];colors = {};i=0
                for shape_dict in m.states_info:
                    names.append(shape_dict['HASC_1'])
                    if len(np.where(np.array(states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                        state_color = varExp.loc[str(crop+'_'+shape_dict['NAME_1'].upper())][mode]
                        colors[shape_dict['HASC_1']]=state_color
                    else: continue#colors[shape_dict['HASC_1']]=0
                ax = plt.gca() # get current axes instance
                for nshape, seg in enumerate(m.states):
                    if names[nshape] not in colors: continue
                    poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                    ax.add_patch(poly)
                mapper.set_array(clrBr)
            
                INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                            name='states', drawbounds=True)
                names = [];colors = {};i=0
                for shape_dict in m.states_info:
                    names.append(shape_dict['HASC_1'])
                    if ((shape_dict['NAME_1']=='Telangana')):
                        colors[shape_dict['HASC_1']]=varExp.loc[unidecode(str(crop+'_Andhra Pradesh'.upper()))][mode]
                    elif ((shape_dict['NAME_1']=='Chhattisgarh')):
                        colors[shape_dict['HASC_1']]=varExp.loc[unidecode(str(crop+'_Madhya Pradesh'.upper()))][mode]
                    elif len(np.where(np.array(states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                        state_color = varExp.loc[unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))][mode]
                        colors[shape_dict['HASC_1']]=state_color
                    else: continue#colors[shape_dict['HASC_1']]=0
                ax = plt.gca() # get current axes instance
                for nshape, seg in enumerate(m.states):
                    if names[nshape] not in colors: continue
                    poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                    ax.add_patch(poly)
                mapper.set_array(clrBr)
                
            m.pcolor(lons,lats,cMsk_all,cmap='binary',zorder=1.5,latlon=True) 
            m.pcolor(lons,lats,cMsk,cmap='tab20c_r',zorder=2,latlon=True) 
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
            if crop is 'wheat':
                CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                                name='states', drawbounds=True)
                AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                            name='states', drawbounds=True)
                INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                            name='states', drawbounds=True)
            if crop is 'maize':
                MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                            name='states', drawbounds=True)
                INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                            name='states', drawbounds=True)
            cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
            cb = ColorbarBase(cax,cmap=thisCMAP,norm=norm,orientation='horizontal')
            mapper.set_array(clrBr);
            fig.set_size_inches(40,30)
            plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/variance explained/'+mode+crop.lower()+'_'+ENSOyr+'_'+varType+'.png')
            plt.close()
            