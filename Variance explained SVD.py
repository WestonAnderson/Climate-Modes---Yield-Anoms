#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 13:55:33 2017

@author: weston
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:33:29 2017

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
thresh = 0.0001 #0.005 = 0.5%

#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#


months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr0har'

exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/GADM_dict.py').read())



clrMax = .5;clrMin=0

thisCMAP =  'YlOrBr'# matplotlib.colors.Colormap(Margot2_4.mpl_colormap)
clrBr = np.arange(clrMin,clrMax+.01,.01)
norm = Normalize(vmin=clrMin, vmax=clrMax, clip=False)
mapper=cm.ScalarMappable(norm=norm, cmap=thisCMAP)

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   Read in the data                #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#    
## ENSO ## 
enso = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/'+
                           'Climate Modes - Yield Anoms/stats regressions/all cropsVar1961_2012_'+ENSOyr+'_SVD2'+SVDmon+'.csv')   
enso=enso[np.isfinite(enso.SVD_E_var)]
enso['index']=[x.upper() for x in enso['index']]    
## IOD ## 
iod = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/'+
                           'Climate Modes - Yield Anoms/stats regressions/all cropsVar1961_2012_'+ENSOyr+'_iodSVD2'+SVDmon+'.csv')   
iod=iod[np.isfinite(iod.SVD_E_var)]
iod['index']=[x.upper() for x in iod['index']] 
## TAV ## 
tav = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/'+
                           'Climate Modes - Yield Anoms/stats regressions/all cropsVar1961_2012_'+ENSOyr+'_tavSVD2'+SVDmon+'.csv')   
tav=tav[np.isfinite(tav.SVD_E_var)] #"SVD_E" is correct for Atlantic and IOD also.
tav['index']=[x.upper() for x in tav['index']]  




eof_states = enso['index']
crops = ['WHEAT','MAIZE','SOY']


if clim is 'E':
    SVDcrop = 'all crops'
    PC = 1 #python indexing
    varExp = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/'+
                               'Climate Modes - Yield Anoms/stats regressions/'+SVDcrop+'Var1961_2012_'+ENSOyr+'_SVD'+str(PC+1)+SVDmon+'.csv')   
    varExp=varExp[np.isfinite(varExp.SVD_E_var)]
    varExp['index']=[x.upper() for x in varExp['index']]    
    eof_states = varExp['index']
    crops = ['WHEAT','MAIZE','SOY']

elif clim is 'NAO':
    SVDcrop = 'wheat'
    PC = 0 #python indexing
    varExp = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/'+
                               'Climate Modes - Yield Anoms/stats regressions/'+NAOmons+SVDcrop+'Var1961_2012_'+ENSOyr+'_naoSVD'+str(PC+1)+SVDmon+'.csv')   
    varExp=varExp[np.isfinite(varExp.SVD_NAO_var)]
    varExp['index']=[x.upper() for x in varExp['index']]    
    eof_states = varExp['index']
    crops = ['WHEAT']
    
for crop in crops:
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
                colors[shape_dict['NAME']]=varExp[varExp['index']==unidecode(str(crop+'_USSR'))][mod+'_'+clim+'_var'].values.astype(float)[0]
        elif ((shape_dict['NAME']=='Belgium')|(shape_dict['NAME']=='Luxembourg')):
            if crop=='SOY':continue
            else:colors[shape_dict['NAME']]=varExp[varExp['index']==unidecode(str(crop+'_BELGIUM-LUXEMBOURG'))][mod+'_'+clim+'_var'].values.astype(float)[0]
        elif ((shape_dict['NAME']=='Croatia')|(shape_dict['NAME']=='Bosnia and Herz.')|(shape_dict['NAME']=='Slovenia')
            |(shape_dict['NAME']=='Macedonia')|(shape_dict['NAME']=='Serbia')):
                colors[shape_dict['NAME']]=varExp[varExp['index']==unidecode(str(crop+'_YUGOSLAVIA'))][mod+'_'+clim+'_var'].values.astype(float)[0]
        elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME'].upper()))[0]) != 0:
            state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['NAME'].upper())][mod+'_'+clim+'_var'].values.astype(float)[0]
            colors[shape_dict['NAME']]=state_color
        #else: colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.cnts):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]),edgecolor='k',zorder=1.3)
        ax.add_patch(poly)
    mapper.set_array(clrBr)    
    
    if clim is 'E':   
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        USshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/USA_adm_shp/USA_adm1', 
                    name='states', drawbounds=True)
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if (shape_dict['HASC_1']=='US.AK')|(shape_dict['HASC_1']=='US.DC')|(shape_dict['HASC_1']=='US.DE')|(shape_dict['HASC_1']=='US.HI'):
                continue#colors[shape_dict['HASC_1']]=0
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['NAME_1'].upper())][mod+'_'+clim+'_var'].values.astype(float)[0]
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
            if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['HASC_1'][3:5]))[0]) != 0:
                state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['HASC_1'][3:5])][mod+'_'+clim+'_var'].values.astype(float)[0]
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
            if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['NAME_1'].upper())][mod+'_'+clim+'_var'].values.astype(float)[0]
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
                    colors[shape_dict['HASC_1']]=varExp[varExp['index']==unidecode(str(crop+'_SICHUAN'))][mod+'_'+clim+'_var'].values.astype(float)[0]
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['NAME_1'].upper())][mod+'_'+clim+'_var'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
            ax.add_patch(poly) 
        mapper.set_array(clrBr)
        
        if crop is 'WHEAT':    
            CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                            name='states', drawbounds=True)
            # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['NAME_1'].upper())][mod+'_'+clim+'_var'].values.astype(float)[0]
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
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['NAME_1'].upper())][mod+'_'+clim+'_var'].values.astype(float)[0]
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
                    colors[shape_dict['HASC_1']]=varExp[varExp['index']==unidecode(str(crop+'_Andhra Pradesh'.upper()))][mod+'_'+clim+'_var'].values.astype(float)[0]
                elif ((shape_dict['NAME_1']=='Chhattisgarh')):
                    colors[shape_dict['HASC_1']]=varExp[varExp['index']==unidecode(str(crop+'_Madhya Pradesh'.upper()))][mod+'_'+clim+'_var'].values.astype(float)[0]
                elif ((shape_dict['NAME_1']=='Jharkhand')):
                    colors[shape_dict['HASC_1']]=varExp[varExp['index']==unidecode(str(crop+'_Bihar'.upper()))][mod+'_'+clim+'_var'].values.astype(float)[0]
                elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['NAME_1'].upper())][mod+'_'+clim+'_var'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,color=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                ax.add_patch(poly)
            mapper.set_array(clrBr)
            
        if crop is 'MAIZE':
            MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['NAME_1'].upper())][mod+'_'+clim+'_var'].values.astype(float)[0]
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
                    colors[shape_dict['HASC_1']]=varExp[varExp['index']==unidecode(str(crop+'_Andhra Pradesh'.upper()))][mod+'_'+clim+'_var'].values.astype(float)[0]
                elif ((shape_dict['NAME_1']=='Chhattisgarh')):
                    colors[shape_dict['HASC_1']]=varExp[varExp['index']==unidecode(str(crop+'_Madhya Pradesh'.upper()))][mod+'_'+clim+'_var'].values.astype(float)[0]
                elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = varExp[varExp['index']==str(crop+'_'+shape_dict['NAME_1'].upper())][mod+'_'+clim+'_var'].values.astype(float)[0]
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
    if crop is 'WHEAT':
        CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                        name='states', drawbounds=True)
        AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                    name='states', drawbounds=True)
        INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                    name='states', drawbounds=True)
    if crop is 'MAIZE':
        MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                    name='states', drawbounds=True)
        INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                    name='states', drawbounds=True)
    cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
    cb = ColorbarBase(cax,cmap=thisCMAP,norm=norm,orientation='horizontal')
    mapper.set_array(clrBr);
    fig.set_size_inches(40,30)
    plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/variance explained/'+mod.lower()+clim+crop.lower()+'_'+ENSOyr+'_SVD'+str(PC+1)+'.png')
    plt.close()