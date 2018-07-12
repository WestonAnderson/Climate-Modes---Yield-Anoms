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
PCs = 0 #python indexing
timePeriod =''# ' 1960-2010' #'' or ' 1960-2010'
thresh = 0.0001 #0.005 = 0.5%
window = 5 #for the variability. smoothing window will be = (window*2 +1)

NAOseason = 'DJF' #DJF, JJA, MAM etc
NAO_crop = 'wheat' #'all crops' #'wheat'#'maize'#'soy'#'maize+soy'
#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

if notes == '_marSVD': ensoMons =[2,3,4,5,6,7,8,9,10,11,12,13]
elif notes == '_janSVD': ensoMons =[0,1,2,3,4,5,6,7,8,9,10,11]
elif notes == '_janSVD_allTrops': ensoMons =[6,7,8]
elif notes == '_janSVD_allTrops_25N25S': ensoMons =[6,7,8]
elif notes == '_junSVD': ensoMons =[5,6,7,8,9,10,11,12,13,14,15,16]
elif notes == '_decSVD': ensoMons =[-1,0,1,2,3,4,5,6,7,8,9,10]
elif notes == '_sepSVD': ensoMons =[-4,-3,-2,-1,0,1,2,3,4,5,6,7]
elif notes == '_16monSVD': ensoMons =[-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11]
elif notes == '_18monSVD': ensoMons =[-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]
elif notes == '_24monSVD': ensoMons =[-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr0har'

# Restrict SSTs to tropical Pacific and Indian Ocean
latMax = 20; latMin = -20; lonMax = 140; lonMin = 35
latMaxENSO = 20.01; latMinENSO = -20.01; lonMaxENSO = 360; lonMinENSO = 0#define the lat lons used in the first SVD

#Read in the moving standard deviation function
    #The 'smooth' argument is half the window size, so for an 11year running std
    # you would enter smooth=5, so that the window = 5*2+1
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/movingStd.py').read())
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/project/Climate Modes - Yield Anoms/IOD state list.py').read())

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
sstLons=sstData.variables['X'][(sstData.variables['X'][:]<=lonMax)&(sstData.variables['X'][:]>=lonMin)]
sstAnom=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-int(yrMin)))&
                        (sstData.variables['T'][:]<=12*((2+int(yrMax))-1960)),:,
                        (sstData.variables['Y'][:]<=latMax)&(sstData.variables['Y'][:]>=latMin),
                        (sstData.variables['X'][:]<=lonMax)&(sstData.variables['X'][:]>=lonMin)]
#sstAnom = signal.detrend(sstAnom,0,'linear')

sstLon = sstLons-180
sstLons,sstLats=np.meshgrid(sstLons,sstLats)

sstAnom[sstAnom==-999]=0
#pull only chosen month SST anomalies
sstAnom = np.squeeze(sstAnom);sstAnom=np.ma.filled(sstAnom,0)

sstAnom=sstAnom[np.repeat(np.array(range(1+int(yrMax)-int(yrMin)))*12,np.shape(ensoMons)[0])+
                np.tile(ensoMons,1+int(yrMax)-int(yrMin)),...]
sstAnom = np.reshape(sstAnom,[1+int(yrMax)-int(yrMin),np.shape(ensoMons)[0],sstAnom.shape[1],sstAnom.shape[2]])
#sstAnom = np.nanmean(sstAnom,1)


#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#   Remove the influence of ENSO and NAO from crop yields               #
#    and remove the influence of ENSO from SST anomalies                #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#read in the SST matrix
if (notes == '_janSVD_allTrops')|(notes == '_janSVD_allTrops_25N25S'):
  naoLoad = '_janSVD'
  ENSOsst = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+
                    '/all crops/sst_all crops_'+ENSOyr+'_recon_2EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy')
  ensoLats=sstData.variables['Y'][(sstData.variables['Y'][:]<=latMaxENSO)&(sstData.variables['Y'][:]>=latMinENSO)]
  ensoLons=sstData.variables['X'][(sstData.variables['X'][:]<=lonMaxENSO)&(sstData.variables['X'][:]>=lonMinENSO)]
  ENSOsst = ENSOsst.reshape([ENSOsst.shape[0],12,ensoLats.shape[0],ensoLons.shape[0]])
  ENSOsst =ENSOsst[:,:,:,np.where((ensoLons<=lonMax)&(ensoLons>=lonMin))[0]]
  ENSOsst =ENSOsst[:,:,np.where((ensoLats<=latMax)&(ensoLats>=latMin))[0],:]
  sstAnom = sstAnom - ENSOsst[:,ensoMons,...]
else:
    print('\n\n no ENSO SSTs removed \n\n')



#remove the influence from the crop yield matrix
ENSOcr = pd.read_pickle('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/all crops/crop_all crops_'+ENSOyr+
     '_recon_2EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.pkl')
#ensoCols = list(set(CRpctAnom.columns).intersection(naoW.columns))
naoW = pd.read_pickle('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/NAO/crop_'+NAOseason+
                   NAO_crop+'_'+ENSOyr+'_recon_1EOFs_'+str(yrMin)+'-'+str(yrMax)+naoLoad+'.pkl')
naoCols = list(set(CRpctAnom.columns).intersection(naoW.columns))

CRpctAnom=CRpctAnom-ENSOcr #remove ENSO
CRpctAnom[naoCols]=CRpctAnom[naoCols]-naoW[naoCols] #remove NAO

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#

wIODstates=[x.upper() for x in wIODstates] #make state names the same
mIODstates=[x.upper() for x in mIODstates] #make state names the same
sIODstates=[x.upper() for x in sIODstates] #make state names the same
wIODstates = ',wheat_'.join(wIODstates).split(',')
wIODstates[0] = str('wheat_'+wIODstates[0])
mIODstates = ',maize_'.join(mIODstates).split(',')
mIODstates[0] = str('maize_'+mIODstates[0])
sIODstates = ',soy_'.join(sIODstates).split(',')
sIODstates[0] = str('soy_'+sIODstates[0])
crIODstates = wIODstates+mIODstates+sIODstates

eof_states = list(set(CRpctAnom.columns).intersection(crIODstates))

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
cropVarDF['SVD_IOD_var']=np.nan


#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   Calculate the SVD                       #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#   
# Calculate SVD using the relative yield anomalies in each country
# Standardize the SST and yield datasets using a single value
# for the standardization      
sigProd = np.nanstd(CRpctAnom)
sigSST = np.std(sstAnom)
wghtProd = 1./sigProd;
wghtSST=1./sigSST
if wghtProd[~np.isfinite(wghtProd)].size > 0:
    wghtProd[~np.isfinite(wghtProd)]=0
if wghtSST[~np.isfinite(wghtSST)].size > 0:
    wghtSST[~np.isfinite(wghtSST)]=0        

#reshape each dataset
sstAnomRav=np.reshape(sstAnom, [sstAnom.shape[0],sstAnom.shape[1]*sstAnom.shape[2]*sstAnom.shape[3]],'C')

#set the data (c) to be a concatenated matrix of crop production and sst anomalies
c = np.append(sstAnomRav,CRpctAnom,1)
c[np.isnan(c)]=0

# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslatSST = np.cos(np.deg2rad(sstLats))
wgts = wghtSST*np.sqrt(coslatSST)
# add on a series of ones to indicate that the production anomalies in the
# EOF should not be production weighted
wghtsRavSST = np.repeat(np.ravel(wgts,'C'),sstAnom.shape[1]) ; 

#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*#
#       Calculate the covariance matrix          # 
#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*# 
#covWghts = np.repeat(wghtsRav[np.newaxis,:],sstAnomRav.shape[0],axis=0)
SST_c_cov = np.dot(np.transpose(sstAnomRav*wghtsRavSST),(CRpctAnom)*wghtProd)
offset = sstAnomRav.shape[1]       
#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*#
#                   Calculate the EOF            # 
#^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*#
# matrices U, Sig, V calculated from data matrix Z and S. The covariance matrix is ZS^T
eofs_sst,eigens,eofs_prod = np.linalg.svd(SST_c_cov)

Ysst = sstAnomRav*wghtsRavSST
Yc = CRpctAnom.values*wghtProd
Xsst = eofs_sst
Xc = eofs_prod.T
sstPC = Ysst.dot(Xsst)
cPC = Yc.dot(Xc)


np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/crop_'+ENSOyr+'_eigens_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',eigens)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/crop_'+ENSOyr+'_PCs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',cPC)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/sst_'+ENSOyr+'_PCs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',sstPC)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/crop_'+ENSOyr+'_EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',eofs_prod)
np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/sst_'+ENSOyr+'_EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',eofs_sst[:,:eofs_prod.shape[0]])


#calculate the reconstructions and the variance explained with up to 10 PCs:
#First calculate the total variance in each dataset
sstTotVar = np.sum(sstAnomRav**2)
cTotVar = np.sum(CRpctAnom.values**2)

cVar = []
sstVar = []
#make trailing zeros for eigenvectors
for iPC in range(PCs+1):
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
    cReconDF.to_pickle('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/crop_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.pkl')
    np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/crop_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',cRecon)
    np.save('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/IOD/sst_'+ENSOyr+'_recon_'+str(iPC+1)+'EOFs_'+str(yrMin)+'-'+str(yrMax)+notes+'.npy',sstRecon)
    
    CRpctAnom2[:] = cRecon
    #reshape the PC and use the std of the time expansion coefficient to convert it to units of the original matrix
    eof_sst = np.reshape(eofs_sst[:,iPC],sstAnom[0,:].shape)*np.std(sstPC[:,iPC])*0.5 #scaled by 0.5 to match ENSO scaling
    eof_prod = eofs_prod[iPC,:]*np.std(cPC[:,iPC])*0.5 #scaled by 0.5 to match ENSO scaling

    eof_sst = -1*eof_sst
    eof_prod = -1*eof_prod

    #Add the SVD values back into the fulld ataframe
    CRpctAnom2=CRpctAnom2.stack().reset_index()
    CRpctAnom2.columns = ['year','state','SVD']
    CRpctAnom2.set_index(['state','year'],drop=True,inplace=True) 
    cropDF2['SVD'] = CRpctAnom2['SVD']
    cropDF2.reset_index(inplace=True) #drop the index back to a column to alter it
    cropDF2.set_index(['state','year'],drop=True,inplace=True) #re-index
    pd.DataFrame.to_csv(cropDF2,'/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/'+
                        'SVD/all crops1961_2012_'+ENSOyr+'_iodSVD'+str(iPC+1)+notes+'.csv')
    #Make the SVD an estimate of production instead of yield anomalies
    SVD = cropDF2['SVD']#*cropDF2['expectedYldGau']*cropDF2['Harvested_Area']

    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
    #                   Plot the Regions                #
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#         
    for crop in ['wheat','maize','soy']:
        eof_states2 = eof_states.copy()
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
        
        fig = plt.figure()
        m = Basemap(projection='kav7',lon_0=0)
        #m.pcolor(sstLons.astype(int),sstLats.astype(int),eof_sst,cmap='RdBu_r',latlon=True,vmin=-0.75,vmax=0.75)
              
        SVDcol = pd.DataFrame(data={'state':eof_states2,'SVD':eof_prod})
        
        
        #Plot everything at the country level first from FAO data before going into subnational data
        world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                name='cnts', drawbounds=True)
        m.pcolor(lons,lats,cMsk_r,cmap='Dark2_r',zorder=1.2,latlon=True)    
        
        names = [];colors = {};i=0
        for shape_dict in m.cnts_info:
            names.append(shape_dict['NAME'])
            if len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME'].upper())))[0]) != 0:
                state_color = SVDcol[SVDcol.state==unidecode(str(crop+'_'+shape_dict['NAME'].upper()))]['SVD'].values.astype(float)[0]
                colors[shape_dict['NAME']]=state_color
            #else: colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.cnts):
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
            if (shape_dict['HASC_1']=='CN.XZ')|(shape_dict['HASC_1']=='CN.TJ'):
                continue#colors[shape_dict['HASC_1']]=0
            elif len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
                state_color = SVDcol[SVDcol.state==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper()))]['SVD'].values.astype(float)[0]
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
        
            INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if ((shape_dict['NAME_1']=='Telangana')):
                    colors[shape_dict['HASC_1']]=SVDcol[SVDcol.state==unidecode(str(crop+'_Andhra Pradesh'.upper()))]['SVD'].values.astype(float)[0]
                elif ((shape_dict['NAME_1']=='Chhattisgarh')):
                    colors[shape_dict['HASC_1']]=SVDcol[SVDcol.state==unidecode(str(crop+'_Madhya Pradesh'.upper()))]['SVD'].values.astype(float)[0]
                elif ((shape_dict['NAME_1']=='Jharkhand')):
                    colors[shape_dict['HASC_1']]=SVDcol[SVDcol.state==unidecode(str(crop+'_Bihar'.upper()))]['SVD'].values.astype(float)[0]
                elif len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
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
        
            INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if ((shape_dict['NAME_1']=='Telangana')):
                    colors[shape_dict['HASC_1']]=SVDcol[SVDcol.state==unidecode(str(crop+'_Andhra Pradesh'.upper()))]['SVD'].values.astype(float)[0]
                elif ((shape_dict['NAME_1']=='Chhattisgarh')):
                    colors[shape_dict['HASC_1']]=SVDcol[SVDcol.state==unidecode(str(crop+'_Madhya Pradesh'.upper()))]['SVD'].values.astype(float)[0]
                elif len(np.where(np.array(eof_states2)==unidecode(str(crop+'_'+shape_dict['NAME_1'].upper())))[0]) != 0:
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
        eof_plt = np.mean(eof_sst,0)
        eof_plt = np.ma.masked_where(eof_plt==0,eof_plt)
        m1=m.pcolormesh(sstLons, sstLats,eof_plt,zorder=1.1,shading='flat',cmap='RdBu_r',latlon=True,vmin=-1.1,vmax=1.1)  
        cb = m.colorbar(m1,"bottom", size="5%", pad="1%");plt.show()
        world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                name='cnts', drawbounds=True)  
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
        cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
        cb = ColorbarBase(cax,cmap=thisCMAP,norm=norm,orientation='horizontal')
        mapper.set_array(clrBr);
        fig.set_size_inches(30,20)
        plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/IOD SVD/'+crop+'_PC'+str(iPC+1)+ENSOyr+notes+'.png')
        plt.close()

    
#plt.figure();plt.plot(100*(eigens[:10]**2/np.sum(eigens**2)),'o');plt.show()
#          
#        
##estimate the variance for each PC
#sstVar = 100*(np.array(sstVar[1:])-np.array(sstVar[:-1]))/sstTotVar
#cVar = 100*(np.array(cVar[1:])-np.array(cVar[:-1]))/cTotVar
#
#fig = plt.figure(); ax1 = plt.subplot(311); ax2 = plt.subplot(312); ax3 = plt.subplot(313)
#ax1.set_title('Percent Variance Explained');
#ax1.plot(range(1,PCs+1),eigens[0:PCs]*100/np.sum(eigens),'ko')
#ax2.plot(range(1,sstVar.shape[0]+1),sstVar,'ko')
#ax3.plot(range(1,cVar.shape[0]+1),cVar,'ko')
#ax1.set_ylabel('cross-covariance');ax2.set_ylabel('SST variance')
#ax3.set_ylabel('Crop yield variance');plt.xlabel('Mode')
#fig.set_size_inches(4,8);
#fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/IOD SVD/EigenVectors_'+ENSOyr+'_toPC'+str(iPC)+notes+'.png')
#plt.close() 
#
