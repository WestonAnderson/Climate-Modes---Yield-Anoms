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
from scipy import signal
import netCDF4
import matplotlib
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #
yrMin = 1981
yrMax = 2011
notes = '_janSVD_allTrops'
monCutoff = 1
PCs = 2 #python indexing
timePeriod =''# ' 1960-2010' #'' or ' 1960-2010'
thresh = 0.0001 #0.005 = 0.5%
window = 5 #for the variability. smoothing window will be = (window*2 +1)
lvl = 1000
NAOseason = 'DJF' #DJF, JJA, MAM etc
NAO_crop = 'wheat' #'all crops' #'wheat'#'maize'#'soy'#'maize+soy'
#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

if notes == '_marSVD': ensoMons =[2,3,4,5,6,7,8,9,10,11,12,13]
elif notes == '_janSVD': ensoMons =[0,1,2,3,4,5,6,7,8,9,10,11]
elif notes == '_janSVD_allTrops': ensoMons =[3,4,5,6]
elif notes == '_janSVD_allTrops_25N25S': ensoMons =[3,4,5,6]
elif notes == '_junSVD': ensoMons =[5,6,7,8,9,10,11,12,13,14,15,16]
elif notes == '_decSVD': ensoMons =[-1,0,1,2,3,4,5,6,7,8,9,10]
elif notes == '_sepSVD': ensoMons =[-4,-3,-2,-1,0,1,2,3,4,5,6,7]
elif notes == '_16monSVD': ensoMons =[-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11]
elif notes == '_18monSVD': ensoMons =[-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]
elif notes == '_24monSVD': ensoMons =[-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr0har'

# Restrict SSTs to tropical Pacific and Indian Ocean
latMax = 20; latMin = -20; lonMax = 300; lonMin = 20
latMaxENSO = 20.01; latMinENSO = -20.01; lonMaxENSO = 360; lonMinENSO = 0#define the lat lons used in the first SVD

#Read in the moving standard deviation function
    #The 'smooth' argument is half the window size, so for an 11year running std
    # you would enter smooth=5, so that the window = 5*2+1
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/movingStd.py').read())
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/project/Climate Modes - Yield Anoms/TAV state list.py').read())

thisCMAP =  'BrBG'# matplotlib.colors.Colormap(Margot2_4.mpl_colormap)
clrBr = np.arange(-0.6,0.6+.01,.01)
norm = Normalize(vmin=-0.6, vmax=0.6, clip=False)
mapper=cm.ScalarMappable(norm=norm, cmap=thisCMAP)

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   Read in the data                        #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#         
#crop data first
wDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/Dataframes/wheat.csv')   
#drop the country statistics that will be included at the state level instead
wDF = wDF[wDF.year>=yrMin]
wDF = wDF[wDF.year<=yrMax]
wDF.state=[x.upper() for x in wDF.state] #make state names the same
wDF=wDF[np.isfinite(wDF.yldAnomGau)] #drop incomplete entries
wDF['state'] = ',wheat_'.join(wDF['state'].values).split(',')
wDF['state'].iloc[0] = str('wheat_'+wDF['state'].iloc[0])
#count up total production before dropping any states for having a short production time series and before adjusting the years for the SVD
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
CRpctAnom = cropDF.pivot(index='year',columns='state',values='yldAnomGau') #table for SVD
CRpctAnom=CRpctAnom.loc[yrMin:yrMax,:] #crop to the correct years
CRpctAnom=CRpctAnom.dropna(axis=1) #drop incomplete columns

#save combined crop dataframe
cropDF.set_index(['state','year'],drop=True,inplace=True) #re-index the data
cropDF['SVD'] = np.nan #create to fill later

#Read in the SST anomaly data
sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
#T is months since 1960-01-01, and I want to start at 1950 with one year lag
sstLats=sstData.variables['Y'][(sstData.variables['Y'][:]<=latMax)&(sstData.variables['Y'][:]>=latMin)]
sstLons=sstData.variables['X'][((sstData.variables['X'][:]>=lonMax)&(sstData.variables['X'][:]<=360))|
                              ((sstData.variables['X'][:]<=lonMin)&(sstData.variables['X'][:]>=0))]
sstAnom=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-int(yrMin)))&
                        (sstData.variables['T'][:]<=12*((2+int(yrMax))-1960)),:,
                        (sstData.variables['Y'][:]<=latMax)&(sstData.variables['Y'][:]>=latMin),
                        ((sstData.variables['X'][:]>=lonMax))|((sstData.variables['X'][:]<=lonMin))]

#read in the netCDF file containing the variables
ncGeo = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/gph.nc',
                        'r',format='NETCDF4') #read in the netCDF file for the variable
geoVar=ncGeo.variables['hgt'][(ncGeo.variables['year'][:]>=yrMin)&(ncGeo.variables['year'][:]<yrMax+2),
                              :,np.where(ncGeo.variables['level'][:]==lvl)[0][0],
                              (ncGeo.variables['latitude'][:]<=latMax)&(ncGeo.variables['latitude'][:]>=latMin),
                              ((ncGeo.variables['longitude'][:]>=lonMax))|((ncGeo.variables['longitude'][:]<=lonMin))]
ncP = netCDF4.Dataset('/Volumes/Data_Archive/Data/Precip/GPCP/detrended/GPCP_precip.nc','r',format='NETCDF4')
prcp = ncP.variables['precip'][(ncP.variables['year'][:]>=yrMin)&(ncP.variables['year'][:]<yrMax+2),:,
                     (ncP.variables['latitude'][:]<=latMax)&(ncP.variables['latitude'][:]>=latMin),
                     ((ncP.variables['longitude'][:]>=lonMax))|((ncP.variables['longitude'][:]<=lonMin))]
#ncSM = netCDF4.Dataset('/Volumes/Data_Archive/Data/SoilMoisture/GLDAS/Noah10_nc/detrended/SoilM.nc','r',format='NETCDF4') 
#sm = np.nansum(ncSM.variables['SoilM'][(ncSM.variables['year'][:]>=yrMin)&(ncSM.variables['year'][:]<yrMax+2),:,:3,:,:],2)
ncMaxT = netCDF4.Dataset('//Volumes/Data_Archive/Data/Temp/BerkeleyEarth/monthly/Complete_TMAX_LatLong1.nc','r',format='NETCDF4') #read in the netCDF file for the variable
maxT = ncMaxT.variables['temperature'][(ncMaxT.variables['time'][:]>=yrMin)&(ncMaxT.variables['time'][:]<yrMax+2),
                        (ncMaxT.variables['latitude'][:]<=latMax)&(ncMaxT.variables['latitude'][:]>=latMin),
                        ((ncMaxT.variables['longitude'][:]>=lonMax))|((ncMaxT.variables['longitude'][:]<=lonMin))]



for m in range(0,12):
    geoVar[:,m,...] = geoVar[:,m,...] - np.nanmean(geoVar[:,m,...],axis=0)
    prcp[:,m,...] = prcp[:,m,...] - np.nanmean(prcp[:,m,...],axis=0)
#    sm[:,m,...] = sm[:,m,...] - np.nanmean(sm[:,m,...],axis=0)
maxT[:,...] = maxT[:,...] - np.nanmean(maxT[:,...],axis=0)
 
geoVar = np.reshape(geoVar,[geoVar.shape[0]*geoVar.shape[1],geoVar.shape[2],geoVar.shape[3]])
prcp = prcp.reshape([prcp.shape[0]*prcp.shape[1],prcp.shape[2],prcp.shape[3]])
#sm= sm.reshape([sm.shape[0]*sm.shape[1],sm.shape[2],sm.shape[3]])

geoLats=ncGeo.variables['latitude'][(ncGeo.variables['latitude'][:]<=latMax)&(ncGeo.variables['latitude'][:]>=latMin)]
geoLons=ncGeo.variables['longitude'][((ncGeo.variables['longitude'][:]>=lonMax))|((ncGeo.variables['longitude'][:]<=lonMin))]

#smLats=ncSM.variables['latitude'][(ncSM.variables['latitude'][:]<=latMax)&(ncSM.variables['latitude'][:]>=latMin)]
#smLons=ncSM.variables['longitude'][((ncSM.variables['longitude'][:]>=lonMax))|((ncSM.variables['longitude'][:]<=lonMin))]

mLats=ncMaxT.variables['latitude'][(ncMaxT.variables['latitude'][:]<=latMax)&(ncMaxT.variables['latitude'][:]>=latMin)]
mLons=ncMaxT.variables['longitude'][((ncMaxT.variables['longitude'][:]>=lonMax))|((ncMaxT.variables['longitude'][:]<=lonMin))]

pLats=ncP.variables['latitude'][(ncP.variables['latitude'][:]<=latMax)&(ncP.variables['latitude'][:]>=latMin)]
pLons=ncP.variables['longitude'][((ncP.variables['longitude'][:]>=lonMax))|((ncP.variables['longitude'][:]<=lonMin))]

mLon = mLons-180
mLons, mLats = np.meshgrid(mLons,mLats)
#smLon = smLons-180
#smLons, smLats = np.meshgrid(smLons,smLats)
pLon = pLons-180
pLons, pLats = np.meshgrid(pLons,pLats)
sstLon = sstLons-180

#sstLons[sstLons<180]=sstLons[sstLons<180]+360
sstLons,sstLats=np.meshgrid(sstLons,sstLats)
geoLon = geoLons-180
geoLons,geoLats=np.meshgrid(geoLons,geoLats)

sstAnom[sstAnom==-999]=0
#pull only chosen month SST anomalies
geoVar = np.squeeze(geoVar);
sstAnom = np.squeeze(sstAnom);
sstAnom=np.ma.filled(sstAnom,0)
maxT =np.ma.filled(maxT,0)
#sm =np.ma.filled(sm,0)

sstAnom=sstAnom[np.repeat(np.array(range(1+int(yrMax)-int(yrMin)))*12,np.shape(ensoMons)[0])+
                np.tile(ensoMons,1+int(yrMax)-int(yrMin)),...]
sstAnom = np.reshape(sstAnom,[1+int(yrMax)-int(yrMin),np.shape(ensoMons)[0],sstAnom.shape[1],sstAnom.shape[2]])

geoVar=geoVar[np.repeat(np.array(range(1+int(yrMax)-int(yrMin)))*12,np.shape(ensoMons)[0])+
                np.tile(ensoMons,1+int(yrMax)-int(yrMin)),...]
geoVar = np.reshape(geoVar,[1+int(yrMax)-int(yrMin),np.shape(ensoMons)[0],geoVar.shape[1],geoVar.shape[2]])

prcp=prcp[np.repeat(np.array(range(1+int(yrMax)-int(yrMin)))*12,np.shape(ensoMons)[0])+
                np.tile(ensoMons,1+int(yrMax)-int(yrMin)),...]
prcp = np.reshape(prcp,[1+int(yrMax)-int(yrMin),np.shape(ensoMons)[0],prcp.shape[1],prcp.shape[2]])

maxT=maxT[np.repeat(np.array(range(1+int(yrMax)-int(yrMin)))*12,np.shape(ensoMons)[0])+
                np.tile(ensoMons,1+int(yrMax)-int(yrMin)),...]
maxT = np.reshape(maxT,[1+int(yrMax)-int(yrMin),np.shape(ensoMons)[0],maxT.shape[1],maxT.shape[2]])

#sm=sm[np.repeat(np.array(range(1+int(yrMax)-int(yrMin)))*12,np.shape(ensoMons)[0])+
#                np.tile(ensoMons,1+int(yrMax)-int(yrMin)),...]
#sm = np.reshape(sm,[1+int(yrMax)-int(yrMin),np.shape(ensoMons)[0],sm.shape[1],sm.shape[2]])


#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#   Remove the influence of ENSO and NAO from crop yields               #
#    and remove the influence of ENSO from SST anomalies                #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#read in the SST matrix
if (notes == '_janSVD_allTrops')|(notes == '_janSVD_allTrops_25N25S'):
  naoLoad = '_janSVD'
  ENSOsst = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/all crops/sst_all crops_'+ENSOyr+'_recon_2EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy')
  ensoLats=sstData.variables['Y'][(sstData.variables['Y'][:]<=latMaxENSO)&(sstData.variables['Y'][:]>=latMinENSO)]
  ensoLons=sstData.variables['X'][(sstData.variables['X'][:]<=lonMaxENSO)&(sstData.variables['X'][:]>=lonMinENSO)]
  ENSOsst = ENSOsst.reshape([ENSOsst.shape[0],12,ensoLats.shape[0],ensoLons.shape[0]])
  tLons = sstData.variables['X'][:]
  tLats = sstData.variables['Y'][(sstData.variables['Y'][:]<=latMax)&(sstData.variables['Y'][:]>=latMin)]
  ENSOsst =ENSOsst[:,:,:,np.where((ensoLons>=lonMax)|(ensoLons<=lonMin))[0]]
  ENSOsst =ENSOsst[:,:,np.where((ensoLats<=latMax)&(ensoLats>=latMin))[0],:]
  sstAnom = sstAnom - ENSOsst[:,ensoMons,...]
else:
    print('\n\n no ENSO SSTs removed \n\n')

tLons, tLats = np.meshgrid(tLons,tLats)
zeroTemp = tLons*np.nan

#remove the influence from the crop yield matrix
ENSOcr = pd.read_pickle('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/all crops/crop_all crops_'+ENSOyr+
     '_recon_2EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.pkl')
#ensoCols = list(set(CRpctAnom.columns).intersection(naoW.columns))
naoW = pd.read_pickle('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/NAO/crop_'+NAOseason+
                   NAO_crop+'_'+ENSOyr+'_recon_1EOFs_'+str(yrMin)+'-'+str(yrMax)+naoLoad+'.pkl')
naoCols = list(set(CRpctAnom.columns).intersection(naoW.columns))

CRpctAnom=CRpctAnom-ENSOcr #remove ENSO
CRpctAnom[naoCols]=CRpctAnom[naoCols]-naoW[naoCols] #remove NAO
CRpctAnom=CRpctAnom.dropna(axis=1) #drop incomplete columns

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#

wTAVstates=[x.upper() for x in TAVstates] #make state names the same
mTAVstates=[x.upper() for x in TAVstates] 
sTAVstates=[x.upper() for x in TAVstates] 
wTAVstates = ',wheat_'.join(wTAVstates).split(',')
wTAVstates[0] = str('wheat_'+wTAVstates[0])
mTAVstates = ',maize_'.join(mTAVstates).split(',')
mTAVstates[0] = str('maize_'+mTAVstates[0])
sTAVstates = ',soy_'.join(sTAVstates).split(',')
sTAVstates[0] = str('soy_'+sTAVstates[0])
crTAVstates = wTAVstates+mTAVstates+sTAVstates

eof_states = list(set(CRpctAnom.columns).intersection(crTAVstates))

CRpctAnom = CRpctAnom[eof_states]

cropDF.reset_index(inplace=True) #re-index the data to alter year
cropDF.set_index(['state'],drop=True,inplace=True)#re-index to subset
cropDF=cropDF.loc[eof_states]#subset
cropDF = cropDF.loc[cropDF.year>=yrMin]
cropDF = cropDF.loc[cropDF.year<=yrMax]
cropDF.reset_index(inplace=True) #drop the index back to a column to alter it
cropDF.set_index(['state','year'],drop=True,inplace=True) #re-index

#Create and save the crop variance dataframe
cropVarDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/'+
                                'stats regressions/all cropsVar1961_2012_'+ENSOyr+notes+'.csv')
cropVarDF['SVD_TAV_var']=np.nan


#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   Calculate the SVD                       #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#   
# Calculate SVD using the relative yield anomalies in each country
# Standardize the SST and yield datasets using a single value
# for the standardization      
wghtProd = 1./np.nanstd(CRpctAnom)
wghtSST = 1./np.std(sstAnom)
wghtGPH = 1./np.std(geoVar)
#wghtSM = 1./np.std(sm)
wghtP = 1./np.std(prcp)
wghtMaxT = 1./np.std(maxT)

if wghtProd[~np.isfinite(wghtProd)].size > 0:
    wghtProd[~np.isfinite(wghtProd)]=0
if wghtSST[~np.isfinite(wghtSST)].size > 0:
    wghtSST[~np.isfinite(wghtSST)]=0        
#if wghtSM[~np.isfinite(wghtSM)].size > 0:
#    wghtSM[~np.isfinite(wghtSM)]=0  
if wghtMaxT[~np.isfinite(wghtMaxT)].size > 0:
    wghtMaxT[~np.isfinite(wghtMaxT)]=0  
    
#reshape each dataset
sstAnomRav=np.reshape(sstAnom, [sstAnom.shape[0],sstAnom.shape[1]*sstAnom.shape[2]*sstAnom.shape[3]],'C')
geoVarRav = np.reshape(geoVar, [geoVar.shape[0],geoVar.shape[1]*geoVar.shape[2]*geoVar.shape[3]],'C')
maxTRav = np.reshape(maxT, [maxT.shape[0],maxT.shape[1]*maxT.shape[2]*maxT.shape[3]],'C')
prcpRav = np.reshape(prcp, [prcp.shape[0],prcp.shape[1]*prcp.shape[2]*prcp.shape[3]],'C')
#smRav = np.reshape(sm, [sm.shape[0],sm.shape[1]*sm.shape[2]*sm.shape[3]],'C')

# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslatSST = np.cos(np.deg2rad(sstLats))
wgtsSST = wghtSST*np.sqrt(coslatSST)
coslatGPH = np.cos(np.deg2rad(geoLats))
wgtsGPH = wghtGPH*np.sqrt(coslatGPH)
coslatP = np.cos(np.deg2rad(pLats))
wgtsP = wghtP*np.sqrt(coslatP)
coslatMaxT = np.cos(np.deg2rad(mLats))
wgtsMaxT = wghtMaxT*np.sqrt(coslatMaxT)
#coslatSM = np.cos(np.deg2rad(smLats))
#wgtsSM = wghtSM*np.sqrt(coslatSM)
# add on a series of ones to indicate that the production anomalies in the
# EOF should not be production weighted
wghtsRavSST = np.repeat(np.ravel(wgtsSST,'C'),sstAnom.shape[1]) ; 
wghtsRavGPH = np.repeat(np.ravel(wgtsGPH,'C'),geoVar.shape[1]) ; 
wghtsRavP = np.repeat(np.ravel(wgtsP,'C'),prcp.shape[1]) ; 
#wghtsRavSM = np.repeat(np.ravel(wgtsSM,'C'),sm.shape[1]) ; 
wghtsRavMaxT = np.repeat(np.ravel(wgtsMaxT,'C'),maxT.shape[1]) ; 

climRav = sstAnomRav#np.append(prcpRav,np.append(maxTRav,np.append(sstAnomRav,geoVarRav,1),1),1)
wghtsRav = wghtsRavSST#np.append(wghtsRavP,np.append(wghtsRavMaxT,np.append(wghtsRavSST,wghtsRavGPH,0),0),0)

#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*#
#       Calculate the covariance matrix          # 
#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*# 
#covWghts = np.repeat(wghtsRav[np.newaxis,:],sstAnomRav.shape[0],axis=0)
SST_c_cov = np.dot(np.transpose(climRav*wghtsRav),(CRpctAnom)*wghtProd)
offset = sstAnomRav.shape[1]       
#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*#
#                   Calculate the EOF            # 
#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*#
# matrices U, Sig, V calculated from data matrix Z and S. The covariance matrix is ZS^T
eofs_sst,eigens,eofs_prod = np.linalg.svd(SST_c_cov)

Ysst = climRav*wghtsRav
Yc = CRpctAnom.values*wghtProd
Xsst = eofs_sst
Xc = eofs_prod.T
sstPC = Ysst.dot(Xsst)
cPC = Yc.dot(Xc)

plt.figure();plt.plot(100*(eigens[:10]**2/np.sum(eigens**2)),'o');plt.show()


np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/crop_'+ENSOyr+'_eigens_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',eigens)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/crop_'+ENSOyr+'_PCs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',cPC)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/sst_'+ENSOyr+'_PCs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',sstPC)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/crop_'+ENSOyr+'_EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',eofs_prod)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/sst_'+ENSOyr+'_EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',eofs_sst[:,:eofs_prod.shape[0]])


#calculate the reconstructions and the variance explained with up to 10 PCs:
#First calculate the total variance in each dataset
sstTotVar = np.sum(sstAnomRav**2)
cTotVar = np.sum(CRpctAnom.values**2)

cVar = []
sstVar = []
#make trailing zeros for eigenvectors
for iPC in range(1,PCs+1):
    CRpctAnom2 = CRpctAnom.copy()
    cropDF2 = cropDF.copy()
    cropVarDF2 = cropVarDF.copy()
    #Reconstruct each field using the eofs
    sstRecon = np.zeros(sstAnomRav.shape)*np.nan
    cRecon = np.zeros(CRpctAnom.values.shape)*np.nan
    for k in range(sstRecon.shape[0]):
        scaleFactorSST = np.diag(np.dot(eofs_sst[:,:iPC+1].T,np.repeat(sstAnomRav[k,:,np.newaxis],repeats=iPC+1,axis=1)))
        scaleFactorC = np.diag(np.dot(eofs_prod[:iPC+1,:],np.repeat(CRpctAnom.values[k,:,np.newaxis],repeats=iPC+1,axis=1)))
        sstRecon[k,:]= np.dot(eofs_sst[:,:iPC+1],scaleFactorSST)
        cRecon[k,:]= np.dot(scaleFactorC,eofs_prod[:iPC+1,:])
    cVar.append(np.sum(cRecon**2))
    sstVar.append(np.sum(sstRecon**2))

    cReconDF = pd.DataFrame(cRecon,columns=eof_states,index=list(range(yrMin,yrMax+1)))
    cReconDF.to_pickle('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/crop_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.pkl')
    np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/crop_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',cRecon)
    np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/sst_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',sstRecon)
    
    CRpctAnom2[:] = cRecon
    #reshape the PC and use the std of the time expansion coefficient to convert it to units of the original matrix
    eof_sst = np.reshape(eofs_sst[:,iPC],sstAnom[0,:].shape)*np.std(sstPC[:,iPC])
    eof_prod = eofs_prod[iPC,:]*np.std(cPC[:,iPC])

    #Add the SVD values back into the fulld ataframe
    CRpctAnom2=CRpctAnom2.stack().reset_index()
    CRpctAnom2.columns = ['year','state','SVD']
    CRpctAnom2.set_index(['state','year'],drop=True,inplace=True) 
    cropDF2['SVD'] = CRpctAnom2['SVD']
    cropDF2.reset_index(inplace=True) #drop the index back to a column to alter it
    cropDF2.set_index(['state','year'],drop=True,inplace=True) #re-index
    pd.DataFrame.to_csv(cropDF2,'/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/'+
                        'SVD/all crops1961_2012_'+ENSOyr+'_tavSVD'+str(iPC+1)+notes+'.csv')
    #Make the SVD an estimate of production instead of yield anomalies
    SVD = cropDF2['SVD']#*cropDF2['expectedYldGau']*cropDF2['Harvested_Area']
        

    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
    #                   Plot the Regions                #
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#                 
    for crop in ['maize','wheat','soy']:
        continue
        eof_states2 = eof_states.copy()
        eof_statesW = [unidecode(x.upper()) for x in eof_states]
        #Correct the names to match the mapping names 
        
        fig = plt.figure()
        m = Basemap(projection='kav7',lon_0=0)
        #m.pcolor(sstLons.astype(int),sstLats.astype(int),eof_sst,cmap='RdBu_r',latlon=True,vmin=-0.75,vmax=0.75)
        
        m1=m.pcolormesh(sstLons, sstLats,np.mean(eof_sst,0),shading='flat',cmap='RdBu_r',latlon=True,vmin=-1.1,vmax=1.1)  
        cb = m.colorbar(m1,"bottom", size="5%", pad="8%");
          
        SVDcol = pd.DataFrame(data={'state':eof_states2,'SVD':eof_prod})

        #Plot everything at the country level first from FAO data before going into subnational data
        world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                name='cnts', drawbounds=True)
        names = [];colors = {};i=0
        for shape_dict in m.cnts_info:
            names.append(shape_dict['NAME'])
            if len(np.where(np.array(eof_statesW)==unidecode(str(crop+'_'+shape_dict['NAME'].upper())))[0]) != 0:
                state_color = SVDcol[eof_statesW==unidecode(str(crop+'_'+shape_dict['NAME']).upper())]['SVD'].values.astype(float)[0]
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
            elif len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
                state_color = SVDcol[SVDcol.state==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))]['SVD'].values.astype(float)[0]
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
            if len(np.where(np.array(eof_states2)==str(crop+'_'+shape_dict['HASC_1'][3:5]))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['HASC_1'][3:5])]['SVD'].values.astype(float)[0]
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
            if len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
                state_color = SVDcol[SVDcol.state==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))]['SVD'].values.astype(float)[0]
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
                if len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k',zorder=1.3)
                ax.add_patch(poly)
            mapper.set_array(clrBr)     
        
        cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
        cb = ColorbarBase(cax,cmap=thisCMAP,norm=norm,orientation='horizontal')
        mapper.set_array(clrBr);
        fig.set_size_inches(30,20)
        plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/TAV SVD/'+crop+'_PC'+str(iPC+1)+ENSOyr+notes+'.png')
        plt.close()

#    CRpctAnom2 = CRpctAnom.copy()
#    cropDF2 = cropDF.copy()
#    cropVarDF2 = cropVarDF.copy()
#    #Reconstruct each field using the eofs
#    sstRecon = np.zeros(sstAnomRav.shape)*np.nan
#    cRecon = np.zeros(CRpctAnom.values.shape)*np.nan
#    for k in range(sstRecon.shape[0]):
#        scaleFactorSST = np.diag(np.dot(eofs_sst[:,:iPC+1].T,np.repeat(sstAnomRav[k,:,np.newaxis],repeats=iPC+1,axis=1)))
#        scaleFactorC = np.diag(np.dot(eofs_prod[:iPC+1,:],np.repeat(CRpctAnom.values[k,:,np.newaxis],repeats=iPC+1,axis=1)))
#        sstRecon[k,:]= eofs_sst[:,:iPC+1].dot(scaleFactorSST)
#        cRecon[k,:]= scaleFactorC.dot(eofs_prod[:iPC+1,:])
#    cVar.append(np.sum(cRecon**2))
#    sstVar.append(np.sum(sstRecon**2))
#  
#    np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/crop_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',cRecon)
#    np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/TAV/sst_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',sstRecon)
#            
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
    # Plot the combination of PCs    
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
    if(iPC==1):
        if (notes=='_janSVD_allTrops_25N25S'):
            e1fac=1*0.5 #scaled by 0.5 to match ENSO scaling
            e2fac=-1*0.5
        elif (notes=='_janSVD_allTrops'):
            e1fac=1*0.5 #scaled by 0.5 to match ENSO scaling
            e2fac=-1*0.5            
        else:
            e1fac=1*0.5
            e2fac=1*0.5
        sstEOF = sstRecon
        subNatEOF = cRecon
        eof_sst = e1fac*np.reshape(eofs_sst[:,iPC-1],sstAnom[0,:].shape)*np.std(sstPC[:,iPC-1])+\
                  e2fac*np.reshape(eofs_sst[:,iPC],sstAnom[0,:].shape)*np.std(sstPC[:,iPC])

        eofTemp = zeroTemp.copy()
        eofTemp[:,tLons[0,:]<=lonMin]=np.mean(eof_sst[:,:,sstLons[0,:]<=lonMin],0)
        eofTemp[:,tLons[0,:]>=lonMax]=np.mean(eof_sst[:,:,sstLons[0,:]>=lonMax],0)
        
        eof_prod = e1fac*eofs_prod[iPC-1,:]*np.std(cPC[:,iPC-1])+\
                   e2fac*eofs_prod[iPC,:]*np.std(cPC[:,iPC])

        #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
        #                   Plot the Regions                #
        #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#         
        for crop in ['maize','wheat','soy']:
            cWght_f = netCDF4.Dataset('/Volumes/Data_Archive/Data/M3_NCDF/5ArcMin/HA_Yld/wheat_HarvAreaYield_NetCDF/wheat_AreaYieldProduction.nc')
            cMsk = np.load(str('/Volumes/Data_Archive/Data/CropMaps_ComonGrid/AllArea/HarvestedArea/5arcMin_ensemble/ensemble_'+crop[0].lower()+'.npy'))
            lats = cWght_f['latitude'][:]
            lons = cWght_f['longitude'][:]
            
            cLats = np.append((lats[0]-lats[1])+lats[0], lats)
            cLats = cLats +(lats[1]-lats[0])/2
            dlon = (lons[1]-lons[0])*(np.pi/180.) #difference in longitude (in radians) remains constant
            R = 6371000. #Radius of the earth in m
            cellArea = 2*np.pi*(R**2)*np.abs(np.sin(cLats[1:]*(np.pi/180.))- 
                                       np.sin(cLats[:-1]*(np.pi/180.)))*dlon #in m^2
            cellArea = np.repeat(cellArea[:,np.newaxis],lons.shape,axis=1)
            haThresh = (cellArea*thresh)*1./10000 #converting m^2 to ha for the threshold to compare with the ha map
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
            eof_statesW = [unidecode(x.upper()) for x in eof_states]
            #Correct the names to match the mapping names   
            
            fig = plt.figure()
            m = Basemap(projection='kav7',lon_0=0)
            #m.pcolor(sstLons.astype(int),sstLats.astype(int),eof_sst,cmap='RdBu_r',latlon=True,vmin=-0.75,vmax=0.75)
                          
            SVDcol = pd.DataFrame(data={'state':eof_states2,'SVD':eof_prod})
#            with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
#                print(SVDcol)
#            break
            #Plot everything at the country level first from FAO data before going into subnational data
            world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                    name='cnts', drawbounds=True)
            m.pcolor(lons,lats,cMsk_r,cmap='Dark2_r',zorder=1.2,latlon=True) 
            names = [];colors = {};i=0
            for shape_dict in m.cnts_info:
                names.append(shape_dict['NAME'])
                if len(np.where(np.array(eof_statesW)==unidecode(str(crop+'_'+shape_dict['NAME']).upper()))[0]) != 0:
                    state_color = SVDcol[np.array(eof_statesW)==unidecode(str(crop+'_'+shape_dict['NAME']).upper())]['SVD'].values.astype(float)[0]
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
                elif len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))]['SVD'].values.astype(float)[0]
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
                if len(np.where(np.array(eof_states2)==str(crop+'_'+shape_dict['HASC_1'][3:5]))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['HASC_1'][3:5])]['SVD'].values.astype(float)[0]
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
                if len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))]['SVD'].values.astype(float)[0]
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
                    if len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
                        state_color = SVDcol[SVDcol.state==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))]['SVD'].values.astype(float)[0]
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
            m1=m.pcolormesh(tLons,tLats,eofTemp,zorder=1.1,shading='flat',cmap='RdBu_r',latlon=True,vmin=-1.25,vmax=1.25)  
            cb = m.colorbar(m1,"bottom", size="5%", pad="8%");

            world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                    name='cnts', drawbounds=True)
            USshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/USA_adm_shp/USA_adm1', 
                        name='states', drawbounds=True)
            BRshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/BRA_adm_shp/BRA_adm1', 
                            name='states', drawbounds=True)
            ARshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/ARG_adm_shp/ARG_adm1', 
                            name='states', drawbounds=True)
            if crop is 'maize':
                MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                            name='states', drawbounds=True)
            
            cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
            cb = ColorbarBase(cax,cmap=thisCMAP,norm=norm,orientation='horizontal')
            mapper.set_array(clrBr);
            fig.set_size_inches(30,20)
            plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/TAV SVD/'+crop+'_uptoPC'+str(iPC+1)+ENSOyr+notes+'.png')
            plt.close()


        
#estimate the variance for each PC
sstVar = 100*(np.array(sstVar[1:])-np.array(sstVar[:-1]))/sstTotVar
cVar = 100*(np.array(cVar[1:])-np.array(cVar[:-1]))/cTotVar

fig = plt.figure(); ax1 = plt.subplot(311); ax2 = plt.subplot(312); ax3 = plt.subplot(313)
ax1.set_title('Percent Variance Explained');
ax1.plot(range(1,11),eigens[0:10]*100/np.sum(eigens),'ko')
ax2.plot(range(1,sstVar.shape[0]+1),sstVar,'ko')
ax3.plot(range(1,cVar.shape[0]+1),cVar,'ko')
ax1.set_ylabel('cross-covariance');ax2.set_ylabel('SST variance')
ax3.set_ylabel('Crop yield variance');plt.xlabel('Mode')
fig.set_size_inches(4,8);
fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/TAV SVD/EigenVectors_'+ENSOyr+'_toPC'+str(iPC+1)+notes+'.png')
plt.close() 

