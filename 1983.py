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
import matplotlib
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #
year = 1983
notes = '_janSVD_allTrops'
PCs = 1 #python indexing
thresh = 0.0001 #0.005 = 0.5%
crop = 'maize'
#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

latMaxENSO = 20; latMinENSO = -20; lonMaxENSO = 360; lonMinENSO = 0#define the lat lons used in the first SVD
e_latMax = 5; e_latMin = -5; e_lonMax = 270; e_lonMin = 160


spST=['Uruguay','Paraguay','Aguascalientes','Campeche','Colima','Chiapas','Distrito Federal',
'Guanajuato','Guerrero','Buenos Aires','Catamarca','Chaco','Cordoba','Corrientes','Entre Rios','Formosa',
'Jujuy','La Pampa','Misiones','Salta','San Luis','Santa Fe','Santiago del Estero','Tucuman',
'Chubut','La Rioja','Mendoza','Neuquen','Rio Negro','Chubut','Santa Cruz','San Juan','Tierra del Fuego',
'RR','RO','AC','AM','AP','PA','MA','PI','CE','RN','PB','PE','AL','SE','BA','MT',
'GO','DF','MG','ES','RJ','TO','MS','SP','PR','SC','RS',
'Hunan','Fujian','Guangdong','Sichuan','Zhejiang','Jiangsu','Anhui','Chongqing',
'Guangxi','Guizhou','Yunnan','Jiangxi','Myanmar','Laos','Thailand','Cambodia','Vietnam',
'Philipines','Singapore','Indonesia','Papua New Guinea','Malaysia','Australia',
'Durango','Prince Edward Island','Nova Scotia', 'New Brunswick', 'Quebec', 'Ontario', 'Manitoba',
'Coahuila','Coahuila de Zaragoza','Saskatchewan', 'Alberta', 'British Columbia',
'South Africa','Zimbabwe','Botswana','Zambia','Tanzania','Malawi','Mozambique','Zambia',
"Côte d'Ivoire",'Senegal','Gambia','Guinea','Guinea-Bissau','Mali','Sierra Leone',
'Ghana','Liberia','Togo','Benin','Burkina Faso','Nigeria','French Guiana','Angola',
'American Samoa', 'Andorra','Barbados','Belize', 'Bermuda', 
'Costa Rica',  'Cuba','Dominica', 'Dominican Rep.',  'El Salvador',
 'Guam', 'Guatemala', 'Guyana','Haiti', 'Honduras','Jamaica','Nicaragua',
 'Panama','Puerto Rico','Singapore', 'Suriname', 'Tonga', 'Trinidad and Tobago',
 'Venezuela','Colombia', 'Peru','Ecuador','Bolivia']

#Read in the SST anomaly data
sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
#T is months since 1960-01-01, and I want to start at 1950 with one year lag
sstLats=sstData.variables['Y'][:]
sstLons=sstData.variables['X'][:]
sstAnom=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-int(year)))&
                        (sstData.variables['T'][:]<=12*((1+int(year))-1960)),:,:,:]
sstAnom=np.squeeze(sstAnom)
obLats = np.where((sstLats[:]>=e_latMin)&(sstLats[:]<=e_latMax))[0]
obLons = np.where((sstLons[:]>=e_lonMin)&(sstLons[:]<=e_lonMax))[0]
obsHM = np.nanmean(sstAnom[:,obLats[0]:obLats[-1],obLons[0]:obLons[-1]],1)

sstLon = sstLons-180
sstLons,sstLats=np.meshgrid(sstLons,sstLats)

#Read in the crop yield anomaly data
crDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/Dataframes/'+crop+'.csv') 
crDF.state=[x.upper() for x in crDF.state] #make state names the same
crDF = crDF[['state','yldAnomGau']].loc[crDF.year==year]

#Load in the reconstructed crop yield and SST anomalies
ENSOcrp = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+
        '/all crops/crop_all crops_Jan-DecENSOyr0har_recon_2EOFs_1981-2011'+notes+'.pkl')
crpCols = [x for x in ENSOcrp.columns if x[0:len(crop)]==crop]
ENSOcrp = ENSOcrp[crpCols].loc[year]
ENSOcrp=ENSOcrp.reset_index()
ENSOcrp['state']=[x.split('_')[1] for x in ENSOcrp['index']]
ENSOcrp=ENSOcrp.rename(columns={year:'ensoYld'})

crDF=ENSOcrp.merge(crDF)#merge the two crop yield dataframes

ENSOsst = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+
        '/all crops/sst_all crops_Jan-DecENSOyr0har_recon_2EOFs_1981-2011'+notes+'.npy')
ensoLats=sstData.variables['Y'][(sstData.variables['Y'][:]<=latMaxENSO)&(sstData.variables['Y'][:]>=latMinENSO)]
ensoLons=sstData.variables['X'][(sstData.variables['X'][:]<=lonMaxENSO)&(sstData.variables['X'][:]>=lonMinENSO)]
eLats = np.where((ensoLats[:]>=e_latMin)&(ensoLats[:]<=e_latMax))[0]
eLons = np.where((ensoLons[:]>=e_lonMin)&(ensoLons[:]<=e_lonMax))[0]

ensoSSTlons = sstLons[0,eLons[0]:eLons[-1]]
ensoSSTlons=ensoSSTlons%180
ensoSSTlons[ensoSSTlons<150]=180-ensoSSTlons[ensoSSTlons<150]
E_W = np.append(np.append(np.repeat('E',10),[''],0),np.repeat('W',44),0)
ensoLabs = [np.str(int(ensoSSTlons[x]))+' '+E_W[x] for x in range(ensoSSTlons.size)]

ENSOsst = ENSOsst.reshape([ENSOsst.shape[0],12,ensoLats.shape[0],ensoLons.shape[0]])
ensoHM = np.nanmean(ENSOsst[year-1981,:,eLats[0]:eLats[-1],eLons[0]:eLons[-1]],1)

ensoLons = ensoLons
ensoLons,ensoLats=np.meshgrid(ensoLons,ensoLats)

states = [unidecode(x.upper()) for x in crDF.state.values] 
firstStates = [unidecode(x.upper()) for x in spST]
firstStates = list(set(crDF.state.values).intersection(firstStates))

ylds = pd.DataFrame(data={'state':crDF.state,'obsYld':crDF.yldAnomGau,'ensoYld':crDF.ensoYld})
yldsFirst = ylds[ylds.state.isin(firstStates)]
#Correct the names to match the mapping names 
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

thisCMAP = 'BrBG'
clrBr = np.arange(-0.5,0.5+.01,.01)
norm = Normalize(vmin=-0.5, vmax=0.5, clip=False)
mapper=cm.ScalarMappable(norm=norm, cmap=thisCMAP)

y_ticks = np.arange(0,12)
y_labs = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']


##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##
#                               Plot
##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##

fig = plt.figure();
ax1=plt.subplot(222);ax2=plt.subplot(224)
ax3=plt.subplot(2,10,5);ax4=plt.subplot(2,10,15)

#First plot the Hovmoller
ax31 = ax3.imshow(obsHM,cmap='RdBu_r',vmin=-2.1,vmax=2.1,aspect='auto',alpha=0.6);
ax4.imshow(ensoHM,cmap='RdBu_r',vmin=-2.1,vmax=2.1,aspect='auto',alpha=0.6);
ax4.set_xticks(range(ensoHM.shape[1])[0::10])
ax4.set_xticklabels(ensoLabs[0::10],fontsize=25,rotation = -90)
ax3.set_yticks(y_ticks[::2]);
ax3.set_yticklabels(y_labs[::2],size=25)
ax4.set_yticks(y_ticks[::2]);
ax4.set_yticklabels(y_labs[::2],size=25)
ax3.set_xticks([])
fig.set_size_inches(50,15)
plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/specific years/'+crop+
            str(year)+'_PCs'+str(PCs+1)+notes+'.png')

##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##
# Plot the observations IN WINTER / SPRING
##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##
plt.sca(ax1)
m = Basemap(projection='cyl',llcrnrlat=-40,urcrnrlat=55,\
            llcrnrlon=-179,urcrnrlon=179,resolution='c',ax=ax1)
#m.pcolor(sstLons.astype(int),sstLats.astype(int),eof_sst,cmap='RdBu_r',latlon=True,vmin=-0.75,vmax=0.75)
      
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
    (shape_dict['NAME']=='Uzbekistan')|(shape_dict['NAME']=='Georgia')):
            colors[shape_dict['NAME']]=ylds[ylds.state==unidecode(str('USSR'))]['obsYld'].values.astype(float)[0]
    elif ((shape_dict['NAME']=='Belgium')|(shape_dict['NAME']=='Luxembourg')):
        if crop=='soy':continue
        else:colors[shape_dict['NAME']]=ylds[ylds.state==unidecode(str('BELGIUM-LUXEMBOURG'))]['obsYld'].values.astype(float)[0]
    elif ((shape_dict['NAME']=='Slovakia')|(shape_dict['NAME']=='Czechia')):
        if crop=='soy':continue
        else:colors[shape_dict['NAME']]=ylds[ylds.state==unidecode(str('CZECHOSLOVAKIA'))]['obsYld'].values.astype(float)[0]
    elif ((shape_dict['NAME']=='Croatia')|(shape_dict['NAME']=='Bosnia and Herz.')|(shape_dict['NAME']=='Slovenia')
        |(shape_dict['NAME']=='Macedonia')|(shape_dict['NAME']=='Serbia')):
            colors[shape_dict['NAME']]=ylds[ylds.state==unidecode(str('YUGOSLAVIA'))]['obsYld'].values.astype(float)[0]
    if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME']).upper()))[0]) != 0:
        state_color = ylds[np.array(states)==unidecode(str(shape_dict['NAME']).upper())]['obsYld'].values.astype(float)[0]
        colors[shape_dict['NAME']]=state_color
    #else: colors[shape_dict['NAME']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.cnts):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
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
    elif len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
        state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['obsYld'].values.astype(float)[0]
        colors[shape_dict['HASC_1']]=state_color
    #else: colors[shape_dict['HASC_1']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.states):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
    ax.add_patch(poly)
mapper.set_array(clrBr)

BRshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/BRA_adm_shp/BRA_adm1', 
                name='states', drawbounds=True)
# collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
names = [];colors = {};i=0
for shape_dict in m.states_info:
    names.append(shape_dict['HASC_1'])
    #if (shape_dict['HASC_1']=='BR.AM'):
    #    continue#colors[shape_dict['HASC_1']]=0
    if len(np.where(np.array(states)==str(shape_dict['HASC_1'][3:5]))[0]) != 0:
        state_color = ylds[ylds.state==str(shape_dict['HASC_1'][3:5])]['obsYld'].values.astype(float)[0]
        colors[shape_dict['HASC_1']]=state_color
    else: continue#colors[shape_dict['HASC_1']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.states):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
    ax.add_patch(poly) 
mapper.set_array(clrBr)

ARshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/ARG_adm_shp/ARG_adm1', 
                name='states', drawbounds=True)
# collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
names = [];colors = {};i=0
for shape_dict in m.states_info:
    names.append(shape_dict['HASC_1'])
    #if (shape_dict['HASC_1']=='AR.DF')|(shape_dict['HASC_1']=='AR.SC')|(shape_dict['HASC_1']=='AR.TF'):
    #    continue#colors[shape_dict['HASC_1']]=0
    if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
        state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['obsYld'].values.astype(float)[0]
        colors[shape_dict['HASC_1']]=state_color
    else: continue#colors[shape_dict['HASC_1']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.states):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
    ax.add_patch(poly) 
mapper.set_array(clrBr)

CNshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CHN_adm_shp/CHN_adm1', 
                name='states', drawbounds=True)
# collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
names = [];colors = {};i=0
for shape_dict in m.states_info:
    names.append(shape_dict['HASC_1'])
    if ((shape_dict['NAME_1']=='Chongqing')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('SICHUAN'))]['obsYld'].values.astype(float)[0]
    elif len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
        state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['obsYld'].values.astype(float)[0]
        colors[shape_dict['HASC_1']]=state_color
    else: continue#colors[shape_dict['HASC_1']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.states):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
    ax.add_patch(poly) 
mapper.set_array(clrBr)

if crop is 'wheat':
    AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                name='states', drawbounds=True)
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['obsYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly) 
    mapper.set_array(clrBr) 

    CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                    name='states', drawbounds=True)
    # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['obsYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly) 

    INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                name='states', drawbounds=True)
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if ((shape_dict['NAME_1']=='Telangana')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Andhra Pradesh'.upper()))]['obsYld'].values.astype(float)[0]
        elif ((shape_dict['NAME_1']=='Chhattisgarh')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Madhya Pradesh'.upper()))]['obsYld'].values.astype(float)[0]
        elif ((shape_dict['NAME_1']=='Jharkhand')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Bihar'.upper()))]['obsYld'].values.astype(float)[0]
        elif len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['obsYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly)
    mapper.set_array(clrBr)   

if crop is 'maize':
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly)
    MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                    name='states', drawbounds=True)
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['obsYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly)
    mapper.set_array(clrBr)
                
    INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                name='states', drawbounds=True)
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if ((shape_dict['NAME_1']=='Telangana')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Andhra Pradesh'.upper()))]['obsYld'].values.astype(float)[0]
        elif ((shape_dict['NAME_1']=='Chhattisgarh')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Madhya Pradesh'.upper()))]['obsYld'].values.astype(float)[0]
        elif len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['obsYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly)
    mapper.set_array(clrBr) 

#m.pcolor(lons,lats,cMsk_all,cmap='binary',zorder=1.5,latlon=True) 
m.pcolor(lons,lats,cMsk,cmap='tab20c_r',zorder=2,latlon=True) 
#m1=m.pcolormesh(sstLons, sstLats,sstAnomJAS,zorder=1.1,shading='flat',cmap='RdBu_r',latlon=True,vmin=-2.1,vmax=2.1,alpha=0.6)  

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
    AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                name='states', drawbounds=True)
    INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                name='states', drawbounds=True)
if crop is 'maize':
    INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                name='states', drawbounds=True)        


##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##
    # Now plot the reconstruction
##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##~*~##
plt.sca(ax2)

m = Basemap(projection='cyl',llcrnrlat=-40,urcrnrlat=55,\
            llcrnrlon=-179,urcrnrlon=179,resolution='c',ax=ax2)
#m.pcolor(sstLons.astype(int),sstLats.astype(int),eof_sst,cmap='RdBu_r',latlon=True,vmin=-0.75,vmax=0.75)
      
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
    (shape_dict['NAME']=='Uzbekistan')|(shape_dict['NAME']=='Georgia')):
            colors[shape_dict['NAME']]=ylds[ylds.state==unidecode(str('USSR'))]['ensoYld'].values.astype(float)[0]
    elif ((shape_dict['NAME']=='Belgium')|(shape_dict['NAME']=='Luxembourg')):
        if crop=='soy':continue
        else:colors[shape_dict['NAME']]=ylds[ylds.state==unidecode(str('BELGIUM-LUXEMBOURG'))]['ensoYld'].values.astype(float)[0]
    elif ((shape_dict['NAME']=='Slovakia')|(shape_dict['NAME']=='Czechia')):
        if crop=='soy':continue
        else:colors[shape_dict['NAME']]=ylds[ylds.state==unidecode(str('CZECHOSLOVAKIA'))]['ensoYld'].values.astype(float)[0]
    elif ((shape_dict['NAME']=='Croatia')|(shape_dict['NAME']=='Bosnia and Herz.')|(shape_dict['NAME']=='Slovenia')
        |(shape_dict['NAME']=='Macedonia')|(shape_dict['NAME']=='Serbia')):
            colors[shape_dict['NAME']]=ylds[ylds.state==unidecode(str('YUGOSLAVIA'))]['ensoYld'].values.astype(float)[0]
    if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME']).upper()))[0]) != 0:
        state_color = ylds[np.array(states)==unidecode(str(shape_dict['NAME']).upper())]['ensoYld'].values.astype(float)[0]
        colors[shape_dict['NAME']]=state_color
    #else: colors[shape_dict['NAME']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.cnts):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
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
    elif len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
        state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['ensoYld'].values.astype(float)[0]
        colors[shape_dict['HASC_1']]=state_color
    #else: colors[shape_dict['HASC_1']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.states):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
    ax.add_patch(poly)
mapper.set_array(clrBr)

BRshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/BRA_adm_shp/BRA_adm1', 
                name='states', drawbounds=True)
# collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
names = [];colors = {};i=0
for shape_dict in m.states_info:
    names.append(shape_dict['HASC_1'])
    #if (shape_dict['HASC_1']=='BR.AM'):
    #    continue#colors[shape_dict['HASC_1']]=0
    if len(np.where(np.array(states)==str(shape_dict['HASC_1'][3:5]))[0]) != 0:
        state_color = ylds[ylds.state==str(shape_dict['HASC_1'][3:5])]['ensoYld'].values.astype(float)[0]
        colors[shape_dict['HASC_1']]=state_color
    else: continue#colors[shape_dict['HASC_1']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.states):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
    ax.add_patch(poly) 
mapper.set_array(clrBr)

ARshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/ARG_adm_shp/ARG_adm1', 
                name='states', drawbounds=True)
# collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
names = [];colors = {};i=0
for shape_dict in m.states_info:
    names.append(shape_dict['HASC_1'])
    #if (shape_dict['HASC_1']=='AR.DF')|(shape_dict['HASC_1']=='AR.SC')|(shape_dict['HASC_1']=='AR.TF'):
    #    continue#colors[shape_dict['HASC_1']]=0
    if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
        state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['ensoYld'].values.astype(float)[0]
        colors[shape_dict['HASC_1']]=state_color
    else: continue#colors[shape_dict['HASC_1']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.states):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
    ax.add_patch(poly) 
mapper.set_array(clrBr)

CNshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CHN_adm_shp/CHN_adm1', 
                name='states', drawbounds=True)
# collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
names = [];colors = {};i=0
for shape_dict in m.states_info:
    names.append(shape_dict['HASC_1'])
    if ((shape_dict['NAME_1']=='Chongqing')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('SICHUAN'))]['obsYld'].values.astype(float)[0]
    elif len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
        state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['obsYld'].values.astype(float)[0]
        colors[shape_dict['HASC_1']]=state_color
    else: continue#colors[shape_dict['HASC_1']]=0
ax = plt.gca() # get current axes instance
for nshape, seg in enumerate(m.states):
    if names[nshape] not in colors: continue
    poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
    ax.add_patch(poly) 
mapper.set_array(clrBr)

if crop is 'wheat':
    AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                name='states', drawbounds=True)
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['ensoYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly) 
    mapper.set_array(clrBr) 

    CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                    name='states', drawbounds=True)
    # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['ensoYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly) 

    INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                name='states', drawbounds=True)
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if ((shape_dict['NAME_1']=='Telangana')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Andhra Pradesh'.upper()))]['ensoYld'].values.astype(float)[0]
        elif ((shape_dict['NAME_1']=='Chhattisgarh')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Madhya Pradesh'.upper()))]['ensoYld'].values.astype(float)[0]
        elif ((shape_dict['NAME_1']=='Jharkhand')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Bihar'.upper()))]['ensoYld'].values.astype(float)[0]
        elif len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['ensoYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly)
    mapper.set_array(clrBr)   

if crop is 'maize':
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly)
    MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                    name='states', drawbounds=True)
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['ensoYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly)
    mapper.set_array(clrBr)
                
    INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                name='states', drawbounds=True)
    names = [];colors = {};i=0
    for shape_dict in m.states_info:
        names.append(shape_dict['HASC_1'])
        if ((shape_dict['NAME_1']=='Telangana')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Andhra Pradesh'.upper()))]['ensoYld'].values.astype(float)[0]
        elif ((shape_dict['NAME_1']=='Chhattisgarh')):
            colors[shape_dict['HASC_1']]=ylds[ylds.state==unidecode(str('Madhya Pradesh'.upper()))]['ensoYld'].values.astype(float)[0]
        elif len(np.where(np.array(states)==unidecode(str(shape_dict['NAME_1'].upper())))[0]) != 0:
            state_color = ylds[ylds.state==unidecode(str(shape_dict['NAME_1'].upper()))]['ensoYld'].values.astype(float)[0]
            colors[shape_dict['HASC_1']]=state_color
        else: continue#colors[shape_dict['HASC_1']]=0
    ax = plt.gca() # get current axes instance
    for nshape, seg in enumerate(m.states):
        if names[nshape] not in colors: continue
        poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
        ax.add_patch(poly)
    mapper.set_array(clrBr) 

#m.pcolor(lons,lats,cMsk_all,cmap='binary',zorder=1.5,latlon=True) 
m.pcolor(lons,lats,cMsk,cmap='tab20c_r',zorder=2,latlon=True) 
#m1=m.pcolormesh(ensoLons, ensoLats,ENSOsstJAS,zorder=1.1,shading='flat',cmap='RdBu_r',latlon=True,vmin=-2.1,vmax=2.1,alpha=0.6)  


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
    AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                name='states', drawbounds=True)
    INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                name='states', drawbounds=True)
if crop is 'maize':
    INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                name='states', drawbounds=True)        
fig.set_size_inches(50,15)
fig.tight_layout()
plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/specific years/'+crop+
            str(year)+'_PCs'+str(PCs+1)+notes+'.png')
plt.close()

