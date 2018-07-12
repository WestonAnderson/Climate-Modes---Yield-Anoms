#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 13:55:33 2017

@author: weston

11/14/17 - This script complements the regression diagnosis of variance epxlained R script
            It reads in that dataframe, calculates an SVD and adds on that column

 EDITS
1/18/18 - Made the SVD for a matrix that contains all month's SST instead of doing it on the annual average
1/23/18 - Updated the SVD to save the expansion coefficients for the crop time series and SST time series   
         
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
import netCDF4
import matplotlib
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #

yrMin = 1981
yrMax = 2011
NAOseason = 'DJF' #DJF, JJA, MAM etc
PCs = 1 #python indexing
monCutoff = 1
timePeriod =''# ' 1960-2010' #'' or ' 1960-2010'
notes = '_janSVD'
ensonotes = '_janSVD_allTrops'
ensoPCs = 2 #not python indexing
NAO_crop = 'wheat' #'all crops' #'wheat'#'maize'#'soy'#'maize+soy'
thresh = 0.0001  #0.005 = 0.5%
NAOsvdOrder = '1st' #'1st' or '2nd'
gphLvl=500

#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#     NAO definition bounding box
latMax = 81; latMin = 19; lonMax = 41; lonMin =269

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
seasons = {'MAMJJA':[[3,4,5,6,7,8]],'JA':[[6,7]],
           'DJF':[[-1,0,1]],'JFM':[[0,1,2]],'FMA':[[1,2,3]],
           'MAM':[[2,3,4]],'AMJ':[[3,4,5]],'MJJ':[[4,5,6]],'JJA':[[5,6,7]],
           'JAS':[[6,7,8]],'ASO':[[7,8,9]],'SON':[[8,9,10]]}
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr0har'
NAOseas =  np.array(seasons[NAOseason])

exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/project/Climate Modes - Yield Anoms/NAO state list.py').read())


thisCMAP =  'BrBG'# matplotlib.colors.Colormap(Margot2_4.mpl_colormap)
clrBr = np.arange(-0.6,0.6+.01,.01)
norm = Normalize(vmin=-0.6, vmax=0.6, clip=False)
mapper=cm.ScalarMappable(norm=norm, cmap=thisCMAP)

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   Read in the data                        #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#   
#crop data
wDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/Dataframes/wheat.csv')   
#drop the country statistics that will be included at the state level instead
wDF = wDF.loc[wDF.state!='Value']
wDF.state=[x.upper() for x in wDF.state] #make state names the same
wDF['state'] = ',wheat_'.join(wDF['state'].values).split(',')
wDF['state'].iloc[0] = str('wheat_'+wDF['state'].iloc[0])
#count up total production before dropping any states for having a short production time series and before adjusting the years for the SVD
wDF=wDF[np.isfinite(wDF.yldAnomGau)] #drop incomplete entries
wDF.state[wDF.state=='wheat_NEW SOUTH WALES(B)']='wheat_NEW SOUTH WALES'

mDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/Dataframes/maize.csv') 
mDF.state=[x.upper() for x in mDF.state] #make state names the same
mDF['state'] = ',maize_'.join(mDF['state'].values).split(',')
mDF['state'].iloc[0] = str('maize_'+mDF['state'].iloc[0])
#count up total production before dropping any states for having a short production time series and before adjusting the years for the SVD
mDF=mDF[np.isfinite(mDF.yldAnomGau)] #drop incomplete entries

sDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/Dataframes/soy.csv')    
sDF.state=[x.upper() for x in sDF.state] #make state names the same
sDF['state'] = ',soy_'.join(sDF['state'].values).split(',')
sDF['state'].iloc[0] = str('soy_'+sDF['state'].iloc[0])
#count up total production before dropping any states for having a short production time series and before adjusting the years for the SVD
sDF=sDF[np.isfinite(sDF.yldAnomGau)] #drop incomplete entries

cropDF =wDF.append(mDF.append(sDF,ignore_index=True),ignore_index=True)
cropDF.set_index(['state','year'],drop=True,inplace=True) #re-index the data
obs = cropDF['yldAnomGau']#*cropDF['expectedYldGau']*cropDF['Harvested_Area']
cropDF.reset_index(inplace=True) #re-index the data to alter year
#cropDF.year = cropDF.year - cropDF.yrAdded #adjust years for SVD matrix

CRpctAnom = cropDF.pivot(index='year',columns='state',values='yldAnomGau') #table for SVD
CRpctAnom=CRpctAnom.loc[yrMin:yrMax,:] #crop to the correct years
CRpctAnom=CRpctAnom.dropna(axis=1) #drop incomplete columns

#Duplicate the DF for the NAO analysis with only the NAO states
NAOstates=[x.upper() for x in NAOstates] #make state names the same
wNAOstates = ',wheat_'.join(NAOstates).split(',')
wNAOstates[0] = str('wheat_'+wNAOstates[0])
mNAOstates = ',maize_'.join(NAOstates).split(',')
mNAOstates[0] = str('maize_'+mNAOstates[0])
sNAOstates = ',soy_'.join(NAOstates).split(',')
sNAOstates[0] = str('soy_'+sNAOstates[0])
crNAOstates = wNAOstates+mNAOstates+sNAOstates
msNAOstates = mNAOstates+sNAOstates

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                     NAO yield data                        #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*# 
### Calculate the NAO yield time series (select NAO countries, remove ENSO influence)
if NAOsvdOrder=='1st':
    CRpctAnomLessENSO = CRpctAnom #Dont remove the impact of ENSO
if NAOsvdOrder=='2nd':
    #Read in ENSO reconstructed with 2 EOFs
    ensoYlds = pd.read_pickle('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/'+
                          'SVD/all crops/crop_all crops_'+ENSOyr+'_recon_'+str(ensoPCs)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+ensonotes+'.pkl')
    st_list = list(set(CRpctAnom.columns).intersection(ensoYlds))
    CRpctAnomLessENSO = CRpctAnom - ensoYlds[st_list] #Remove the impact of ENSO
    
#subset the crop dataframe to only be NAO states
if NAO_crop is 'all crops':
    eof_states = list(set(CRpctAnomLessENSO.columns).intersection(crNAOstates))
    crops = ['wheat','maize','soy']
elif NAO_crop is 'wheat':
    eof_states=list(set(CRpctAnomLessENSO.columns).intersection(wNAOstates))#find the states 
    crops = ['wheat']
elif NAO_crop is 'maize':
    eof_states = list(set(CRpctAnomLessENSO.columns).intersection(mNAOstates))
    crops = ['maize']
elif NAO_crop is 'soy':
    eof_states = list(set(CRpctAnomLessENSO.columns).intersection(sNAOstates))
    crops = ['soy']
elif NAO_crop is 'maize+soy':
    eof_states = list(set(CRpctAnomLessENSO.columns).intersection(msNAOstates))
    crops = ['maize','soy']
else:
    print('use valid crop selection')
    exit
CRpctAnomNao = CRpctAnomLessENSO[eof_states]

cropDF.set_index(['state'],drop=True,inplace=True) 
cropDF=cropDF.loc[eof_states]
cropDF = cropDF.loc[cropDF.year>=yrMin]
cropDF = cropDF.loc[cropDF.year<=yrMax]
cropDF.reset_index(inplace=True) #drop the index back to a column to alter it
cropDF.set_index(['state','year'],drop=True,inplace=True) #re-index

#Create and save the crop variance dataframe
cropVarDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/'+
                                'stats regressions/all cropsVar1961_2012_Jan-DecENSOyr0har_janSVD_allTrops_template.csv')
cropVarDF['SVD_NAO_var']=np.nan

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#           Read in the 500 hPa GPH data
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#read in the netCDF file containing the variables
ncGeo = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/gph.nc',
                        'r',format='NETCDF4') #read in the netCDF file for the variable
geoVar=ncGeo.variables['hgt'][(ncGeo.variables['year'][:]>=yrMin-1)&(ncGeo.variables['year'][:]<yrMax+2),
                              :,np.where(ncGeo.variables['level'][:]==gphLvl)[0][0],
                              (ncGeo.variables['latitude'][:]<latMax)&(ncGeo.variables['latitude'][:]>latMin),
                              (ncGeo.variables['longitude'][:]<lonMax)|(ncGeo.variables['longitude'][:]>lonMin)]
#ncGeo = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/slp.nc',
#                        'r',format='NETCDF4') #read in the netCDF file for the variable
#geoVar=ncGeo.variables['slp'][(ncGeo.variables['year'][:]>=yrMin-1)&(ncGeo.variables['year'][:]<yrMax+2),
#                              :,(ncGeo.variables['latitude'][:]<latMax)&(ncGeo.variables['latitude'][:]>latMin),
#                              (ncGeo.variables['longitude'][:]<lonMax)|(ncGeo.variables['longitude'][:]>lonMin)]


geoLons = ncGeo.variables['longitude'][(ncGeo.variables['longitude'][:]<lonMax)|(ncGeo.variables['longitude'][:]>lonMin)]
lonSplt = np.min(np.where(geoLons>lonMin))
geoLons = np.append(geoLons[lonSplt:],geoLons[:lonSplt],0)#reorder for nice plotting
geoLats = ncGeo.variables['latitude'][(ncGeo.variables['latitude'][:]<latMax)&(ncGeo.variables['latitude'][:]>latMin)]
geoLons,geoLats=np.meshgrid(geoLons,geoLats)

#Remove the climatology to get anomalies by month
for m in range(0,12):
    geoVar[:,m,...] = geoVar[:,m,...] - np.nanmean(geoVar[:,m,...],axis=0)
#Rearrange to have only one time dimension
geoVar = geoVar.reshape([geoVar.shape[0]*geoVar.shape[1],geoVar.shape[2],geoVar.shape[3]])
geoVar = np.reshape(geoVar,[geoVar.shape[0],geoVar.shape[1]*geoVar.shape[2]])

indYrs = np.repeat(np.array(range(1,2+int(yrMax)-int(yrMin)))*12,np.shape(NAOseas)[1])
indMons = np.tile(NAOseas,1+int(yrMax)-int(yrMin))
for iSeas in range(NAOseas.shape[0]):
    ind = indYrs+indMons[iSeas]
    if iSeas==0:
        geoVarNAO = geoVar[ind,...]
        geoVarNAO = geoVarNAO.reshape([1+int(yrMax)-int(yrMin),NAOseas.shape[1],geoVarNAO.shape[1]])
        geoVarNAO = np.mean(geoVarNAO,1)
    else:
        geoVarNAO_t = geoVar[ind,...]
        geoVarNAO_t = geoVarNAO_t.reshape([1+int(yrMax)-int(yrMin),NAOseas.shape[1],geoVarNAO_t.shape[1]])
        geoVarNAO_t = np.mean(geoVarNAO_t,1)
        geoVarNAO = np.append(geoVarNAO,geoVarNAO_t,1)

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   Calculate the SVD                       #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#   
# Calculate SVD using the relative yield anomalies in each country
# Standardize the SST and yield datasets using a single value
# for the standardization      
sigProd = np.nanstd(CRpctAnomNao)
sigNAO = np.std(geoVarNAO)
wghtProd = 1./sigProd;
wghtNAO=1./sigNAO
if wghtProd[~np.isfinite(wghtProd)].size > 0:
    wghtProd[~np.isfinite(wghtProd)]=0
if wghtNAO[~np.isfinite(wghtNAO)].size > 0:
    wghtNAO[~np.isfinite(wghtNAO)]=0        


#set the data (c) to be a concatenated matrix of crop production and sst anomalies
c = np.append(geoVarNAO,CRpctAnomNao,1)
c[np.isnan(c)]=0


#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*#
#       Calculate the covariance matrix          # 
#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*# 
#covWghts = np.repeat(wghtsRav[np.newaxis,:],sstAnomRav.shape[0],axis=0)
NAO_c_cov = np.dot(np.transpose(geoVarNAO*wghtNAO),(CRpctAnomNao)*wghtProd)
offset = geoVarNAO.shape[1]       
#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*#
#                   Calculate the EOF            # 
#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*#
# matrices U, Sig, V calculated from data matrix Z and S. The covariance matrix is ZS^T
eofs_nao,eigens,eofs_prod = np.linalg.svd(NAO_c_cov)

Ynao = geoVarNAO*wghtNAO
Yc = CRpctAnomNao.values*wghtProd
Xnao = eofs_nao
Xc = eofs_prod.T
naoPC = Ynao.dot(Xnao)
cPC = Yc.dot(Xc)

np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/crop_'+NAOseason+NAO_crop+'_'+ENSOyr+'_eigens_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',eigens)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/crop_'+NAOseason+NAO_crop+'_'+ENSOyr+'_PCs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',cPC)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/gph_'+NAOseason+NAO_crop+'_'+ENSOyr+'_PCs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',naoPC)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/crop_'+NAOseason+NAO_crop+'_'+ENSOyr+'_EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',eofs_prod)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/gph_'+NAOseason+NAO_crop+'_'+ENSOyr+'_EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',eofs_nao[:,:eofs_prod.shape[0]])


#calculate the reconstructions and the variance explained with up to 10 PCs:
#First calculate the total variance in each dataset
naoTotVar = np.sum(geoVarNAO**2)
cTotVar = np.sum(CRpctAnomNao.values**2)

cVar = []
naoVar = []

#make trailing zeros for eigenvectors
for iPC in range(PCs+1):
    CRpctAnomNao2 = CRpctAnomNao.copy()
    cropDF2 = cropDF.copy()
    cropVarDF2 = cropVarDF.copy()
    #Reconstruct each field using the eofs
    naoRecon = np.zeros(geoVarNAO.shape)*np.nan
    cRecon = np.zeros(CRpctAnomNao.values.shape)*np.nan
    for k in range(naoRecon.shape[0]):
        scaleFactorNAO = np.diag(np.dot(eofs_nao[:,:iPC+1].T,np.repeat(geoVarNAO[k,:,np.newaxis],repeats=iPC+1,axis=1)))
        scaleFactorC = np.diag(np.dot(eofs_prod[:iPC+1,:],np.repeat(CRpctAnomNao.values[k,:,np.newaxis],repeats=iPC+1,axis=1)))
        naoRecon[k,:]= np.dot(eofs_nao[:,:iPC+1],scaleFactorNAO)
        cRecon[k,:]= np.dot(scaleFactorC,eofs_prod[:iPC+1,:])
    cVar.append(np.sum(cRecon**2))
    naoVar.append(np.sum(naoRecon**2))

    cReconDF = pd.DataFrame(cRecon,columns=eof_states,index=list(range(yrMin,yrMax+1)))
    cReconDF.to_pickle('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/crop_'+NAOseason+NAO_crop+'_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.pkl')
    print('wrote recon')
    np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/crop_'+NAOseason+NAO_crop+'_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',cRecon)
    np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/gph_'+NAOseason+NAO_crop+'_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',naoRecon)
    
    naoEOF = naoRecon
    CRpctAnomNao2[:] = cRecon
    #reshape the PC and use the std of the time expansion coefficient to convert it to units of the original matrix
    eof_nao = np.reshape(eofs_nao[:,iPC],[geoLats.shape[0]*NAOseas.shape[0],geoLats.shape[1]])*np.std(naoPC[:,iPC])*0.5 #scaled by 0.5 to match ENSO scaling
    eof_prod = eofs_prod[iPC,:]*np.std(cPC[:,iPC])*0.5 #scaled by 0.5 to match ENSO scaling
    #reorder for nice plotting, the same is done with the longitudes
    eof_nao=np.append(eof_nao[:,lonSplt:],eof_nao[:,:lonSplt],1)

    #Add the SVD values back into the fulld ataframe
    CRpctAnomNao2=CRpctAnomNao2.stack().reset_index()
    CRpctAnomNao2.columns = ['year','state','SVD']
    CRpctAnomNao2.set_index(['state','year'],drop=True,inplace=True) 
    cropDF2['SVD'] = CRpctAnomNao2['SVD']
    cropDF2.reset_index(inplace=True) #re-index the data to alter year
#    cropDF2.year = cropDF2.year + cropDF2.yrAdded #recover the original year
    cropDF2.set_index(['state','year'],drop=True,inplace=True) #re-index
    pd.DataFrame.to_csv(cropDF2,'/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/'+
                        'SVD/'+NAOseason+NAO_crop+'1961_2012_'+ENSOyr+'_naoSVD'+str(iPC+1)+notes+'.csv')    
    #Make the SVD an estimate of production instead of yield anomalies
    SVD = cropDF2['SVD']#*cropDF2['expectedYldGau']*cropDF2['Harvested_Area']

    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
    #                   Plot the Regions                #
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#         
    for crop in crops:#['wheat','maize','soy']:
        
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
        
        eof_states2 = eof_states.copy()
        #Correct the names to match the mapping names 
        
        fig = plt.figure()
        m = Basemap(projection='kav7',lon_0=0)
        #m.pcolor(sstLons.astype(int),sstLats.astype(int),eof_sst,cmap='RdBu_r',latlon=True,vmin=-0.05,vmax=0.05)
   
        
        SVDcol = pd.DataFrame(data={'state':eof_states,'SVD':eof_prod})
#       print(SVDcol)
#       break
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
                    #print(SVDcol[SVDcol.state==unidecode(str(crop+'_USSR'))]['SVD'].values.astype(float)[0])
                    colors[shape_dict['NAME']]=SVDcol[SVDcol.state==unidecode(str(crop+'_USSR'))]['SVD'].values.astype(float)[0]
            elif ((shape_dict['NAME']=='Slovakia')|(shape_dict['NAME']=='Czechia')):
                if crop=='soy':continue
                else:colors[shape_dict['NAME']]=SVDcol[SVDcol.state==unidecode(str(crop+'_CZECHOSLOVAKIA'))]['SVD'].values.astype(float)[0]
            elif ((shape_dict['NAME']=='Belgium')|(shape_dict['NAME']=='Luxembourg')):
                if crop=='soy':continue
                else:colors[shape_dict['NAME']]=SVDcol[SVDcol.state==unidecode(str(crop+'_BELGIUM-LUXEMBOURG'))]['SVD'].values.astype(float)[0]
            elif ((shape_dict['NAME']=='Croatia')|(shape_dict['NAME']=='Bosnia and Herz.')|(shape_dict['NAME']=='Slovenia')
                |(shape_dict['NAME']=='Macedonia')|(shape_dict['NAME']=='Serbia')):
                    colors[shape_dict['NAME']]=SVDcol[SVDcol.state==unidecode(str(crop+'_YUGOSLAVIA'))]['SVD'].values.astype(float)[0]
            elif len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME'].upper())))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['NAME']]=state_color
            #else: colors[shape_dict['HASC_1']]=0
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
            elif len(np.where(np.array(eof_states2)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            #else: colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
            ax.add_patch(poly)
        mapper.set_array(clrBr)
        
        if crop is 'wheat':
            CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                            name='states', drawbounds=True)
            # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states2)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                ax.add_patch(poly)  

        m.pcolor(lons,lats,cMsk_all,cmap='binary',zorder=1.5,latlon=True) 
        m.pcolor(lons,lats,cMsk,cmap='tab20c_r',zorder=2,latlon=True) 
        world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                name='cnts', drawbounds=True)
        USshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/USA_adm_shp/USA_adm1', 
                    name='states', drawbounds=True)
        if crop is 'wheat':
            CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                            name='states', drawbounds=True)      
        
        m.contour(geoLons, geoLats,eof_nao,np.arange(-1.5,1.6,.2), shading='flat',colors='k',latlon=True,alpha=0.75)          

        cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
        cb = ColorbarBase(cax,cmap=thisCMAP,norm=norm,orientation='horizontal')
        mapper.set_array(clrBr);
        fig.set_size_inches(30,20)
        plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/NAO SVD/'+NAOseason+NAO_crop+'SVD_'+crop+'_PC'+str(iPC+1)+ENSOyr+notes+'_NAO'+NAOsvdOrder+'.png')
        plt.close()

      
        
#estimate the variance for each PC
naoVar = 100*(np.array(naoVar[1:])-np.array(naoVar[:-1]))/naoTotVar
cVar = 100*(np.array(cVar[1:])-np.array(cVar[:-1]))/cTotVar

fig = plt.figure(); ax1 = plt.subplot(311); ax2 = plt.subplot(312); ax3 = plt.subplot(313)
ax1.set_title('Percent Variance Explained');
ax1.plot(range(1,PCs+1),eigens[0:PCs]*100/np.sum(eigens),'ko')
ax2.plot(range(1,naoVar.shape[0]+1),naoVar,'ko')
ax3.plot(range(1,cVar.shape[0]+1),cVar,'ko')
ax1.set_ylabel('cross-covariance');ax2.set_ylabel('SLP variance')
ax3.set_ylabel('Crop yield variance');plt.xlabel('Mode')
fig.set_size_inches(4,8);
fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/NAO SVD/'+NAOseason+NAO_crop+'_EigenVectors_'+ENSOyr+'_toPC'+str(iPC)+notes+'_NAO'+NAOsvdOrder+'.png')
plt.close() 
  
