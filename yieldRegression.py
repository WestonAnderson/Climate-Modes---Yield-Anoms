# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 13:20:47 2015
@author: westonanderson

Construct regression models using state-wise data and aggregated anomalies.
The mean of each anomaly (including sum precip) is used across the state to match yields

precipitation from the winter season is used directly rather than linearly
relating it first to soil moisture in the spring

EDIT 05/12/16: CA and MX added, but only yield, production and harvested area data calculated
                The climate data require a subnational grid to aggregate onto, and so are skipped
                because the MIRCA data doesn't include them.
                
EDIT 06/02/16: fixed the cropland mask portion of the script so that variables could be binned
                in only the cropland regions.

EDIT 06/15/16: added the potential for soil water content (SWC) to be turned on. Intended for BR only

EDIT 07/05/16: added factors to convert everything to standard units of
    YIELD: kg/ha
    HARVESTED AREA: ha
    PRODUCTION: kg
    
EDIT 07/08/16: removed the ENSO specific variables, made the whole script transferable to a remote machine

EDIT 07/15/16: added an absolute production anomaly, included gaussian smoothing to calculate yield anomalies

EDIT: 10/2/17: Converting this script to be used in the global modes of climate variability. This script
                is still only for the subnational statistics data
                
EDIT: 10/31/17 - SPOOKY - Updated to pull harvest days from Sacks (2010)
                        - Correct Argentina and Brazil year reported (which is planting)
                            to be for harvesting.
EDIT: 11/17/17 - The spooky update had a typo. Fixed the harvest/plant mismatch issue
                    and also fixed the EN/LN year variable to match      

EDIT: 1/18/18 - Included climate variables averaged over the 3 months prior to harvest
"""

import numpy as np
import pandas as pd
import time
from scipy import signal
from scipy import ndimage
from scipy import stats
import netCDF4
import copy
start = time.clock()

def moving_average(a, smooth) :
    n=smooth*2+1
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

crops = ['soy']#'wheat','maize','soy','winterwheat','springwheat']#,
stORreg = 'state' #state or region or, if using US full dataset, stat_long
countries = ['IND']#
dt = '' # '_dt' or ''
dep = 1 #SM depth for beginning of the flowering season soil moisture 0: 0-10cm, 1:10-40cm, 2:40-100, 3:100-200m
notes = '' #If notes are 'long', then it's for long US record
loc = 'local' #remote or local
thresh = 0.001 #0.005 = 0.5%

stYr = [1950];endYr = [2010] #start and end years for the climate analysis

if loc is 'remote':
    path = '/home/weston/'
elif loc is 'local':
    path = '/Volumes/Data_Archive/'

FAO = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/Yield/HistoricalData/FAOSTAT/FAOSTAT_2014.csv',header=0)
MircaCodes = np.load('/Volumes/Data_Archive/Data/CropCalendar/CCCmask/1.0Deg/UnitCode/SubnatUnitCode.npy')
CntryCode = np.load('/Volumes/Data_Archive/Data/CropCalendar/CCCmask/UnitCode/CntryUnitCode.npy', encoding='latin1')
adminCodes = pd.read_table('/Volumes/Data_Archive/Data/CropCalendar/MIRCA2000/unit_name_FAO.txt',sep='\t')
#re-parse the unit code names
for n in range(0,np.size(adminCodes['UnitName'])):
    if np.size(adminCodes['UnitName'][n].split('_'))>1:
        newName = adminCodes['UnitName'][n].split('_')[1]
        adminCodes['UnitName'][n]=newName
msk = netCDF4.Dataset('/Volumes/Data_Archive/Data/landSeaMask/1.0Deg/lsmask.nc')
LSmsk = msk.variables['mask'][0,...];
LSmsk = np.append(LSmsk[:,180:],LSmsk[:,:180],axis=1)



for crop in crops:
    print( crop)
    if crop=='springwheat':
        cropNam = 'w_sp'
        ccCrop = 'Wheat.Main'
        climCrop = 'winterwheat'
    elif crop=='winterwheat':
        cropNam = 'w_win'
        ccCrop = 'Wheat.Main'
        climCrop = 'winterwheat'
    elif crop=='maize':
        ccCrop = 'Maize'
        cropNam = crop[0]
        climCrop = 'maize'
    elif crop=='soy':
        ccCrop = 'Soybeans'
        cropNam = crop[0]
        climCrop = 'soy'
    elif crop=='wheat':
        ccCrop = 'Wheat.Main'
        cropNam = crop[0]
        climCrop = 'winterwheat'

       
    #read in land-sea masks
    msk_gauss = netCDF4.Dataset(str(path+'Data/landSeaMask/gaussGrid/lsmask.nc'))
    msk_1 = netCDF4.Dataset(str(path+'Data/landSeaMask/1.0Deg/lsmask.nc'))
    LSmskGauss = msk_gauss.variables['land'][0,...]
    LSmsk1deg = msk_1.variables['mask'][0,...];
                    
    #Read in the climate data
    smNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/globalCorrelatedRisks/IntegratedVars'+dt+'/preHarvest/SoilM'+climCrop+'.nc')
    kddNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/globalCorrelatedRisks/IntegratedVars'+dt+'/preHarvest/KDD'+climCrop+'.nc')
    sumpNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/globalCorrelatedRisks/IntegratedVars'+dt+'/preHarvest/sumP'+climCrop+'.nc')

    #<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#
    # calculate the crop mask
    lats = sumpNC.variables['latitude'][:]
    lons = sumpNC.variables['longitude'][:]
    cMsk = np.load(str('/Volumes/Data_Archive/Data/CropMaps_ComonGrid/AllArea/HarvestedArea/gpcc_grid_ensemble/ensemble_'+climCrop[0]+'.npy'))
    cLats = np.append((lats[0]-lats[1])+lats[0], lats)
    cLats = cLats +(lats[1]-lats[0])/2
    dlon = (lons[1]-lons[0])*(np.pi/180.) #difference in longitude (in radians) remains constant
    R = 6371000. #Radius of the earth in m
    cellArea = 2*np.pi*(R**2)*np.abs(np.sin(cLats[1:]*(np.pi/180.))-
                               np.sin(cLats[:-1]*(np.pi/180.)))*dlon #in m
    cellArea = np.repeat(cellArea[:,np.newaxis],lons.shape,axis=1)
    haThresh = (cellArea*thresh)*1./10000 #converting the threshhold to hectares
    cMsk[cMsk<haThresh]=0;cMsk[cMsk>=haThresh]=1

    ENyrs = [1953, 1957, 1963, 1965, 1968, 1972, 1976, 1979, 1982, 1986, 1991, 1994, 1997, 2002, 2004, 2006, 2009] #EN
    LNyrs = [1954, 1964, 1970, 1973, 1983, 1988, 1995, 1998, 2007, 2010] #LN0_2

    yearNC = np.array(range(1900,2016))
    colNams = [ 'country','country name','region','flwrSeas',\
                'yield','Production','Harvested_Area',\
                'yldAnom','yldAnomFD','yldAnomGau','yldAnomSmooth5','yldAnomSmooth9',\
                'yldAnomAbs','yldAnomSmooth5Abs','yldAnomSmooth9Abs','yldAnomGauAbs',\
                'expectedYld','expectedYldSmooth5','expectedYldSmooth9','expectedYldGau',\
                'prodAnomSmooth5','prodAnomSmooth9','prodAnomGau','directProdAnom','directProdExp',\
                'expectedProd','expectedProdSmooth5','expectedProdSmooth9','expectedProdGau',\
                'GDD','KDD','KDD_sum','sumPan','sumPan_sum','sumP','sumP2','sumPanLag','PET','P_E',\
                'VPD','soilM','soilMseed','EN','LN','LN2']
        
    mismatch = []
    
    for country in countries:
        print( country)
        yldFactor = 1.;prdFactor = 1.;HAfactor = 1.
        if 'sumpAnom_lag' in locals(): 
            del sumpAnom_lag
        #define the country parameters
        if country is 'US':
            cntName = 'United States'
            dataFile = 'USDA_NASS'
            breakPt = 30 #number of years in to put the breakpoint for detrending
            
            if crop=='wheat':
                prodFactor = 1000./36.7440 #Bu to tonnes to kg
                yldFactor= (1000.*2.47)/36.7440 #Bu/Acre to Kg/ha
                HAfactor = 1./2.47 #acres to ha
                growMon = [np.array([4,5,6])] #3,7for temperature
            elif crop=='winterwheat':
                prodFactor = 1000./36.7440
                yldFactor = (1000.*2.47)/36.7440
                HAfactor = 1./2.47
                growMon = [np.array([4,5,6])] #3,7for temperature
            elif crop=='springwheat':
                prodFactor = 1000./36.7440
                yldFactor = (1000.*2.47)/36.7440
                HAfactor = 1./2.47
                growMon = [np.array([5,6,7])] #3,7for temperature
            elif crop=='maize':
                prodFactor = 1000./39.3680
                yldFactor = (1000.*2.47)/39.3680
                HAfactor = 1./2.47
                growMon = [np.array([6,7,8])] #9,10for temperature
            elif crop=='soy':
                prodFactor = 1000./36.7440
                yldFactor = (1000.*2.47)/36.7440
                HAfactor = 1./2.47
                growMon = [np.array([6,7,8])] #9,10for temperature
            yldSt = [1900]; 
            yldEnd =[2015]
            
        elif country is 'AR':
            cntName = 'Argentina'
            dataFile = 'AR_SIIA'
            prodFactor = 1000. #convert tonnes to kg
            yldFactor = 1000.#convert tonnes/ha to kg/ha
            breakPt = 30 #number of years in to put the breakpoint for detrending
            if crop=='soy':
                growMon = [np.array([12,13,14])] #15,16 for temperature
                yldSt = [1969]
                yldEnd =[2013]
            elif crop=='maize':
                growMon = [np.array([11,12,13])] #14,15
                yldSt = [1969]
                yldEnd =[2013]
            elif (crop=='wheat')|(crop=='springwheat'):
                growMon = [np.array([9,10,11])] #8,12
                yldSt = [1969]
                yldEnd =[2013]
    
        elif country is 'BR':
            cntName = 'Brazil'
            dataFile = 'BR_CONAB'
            prodFactor = 1000.*1000 #convert thousands of tonnes to kg
            yldFactor = 1.#strangely, BR claims to be in (Em kg/ha) but is clearly (kg/ha)
            HAfactor = 1000.#convert thousands of ha to ha
            if crop=='soy':
                growMon = [np.array([13,14,15])] #
                yldSt = [1978]
                yldEnd = [2013]
                breakPt = 7 #number of years in to put the breakpoint for detrending
            elif (crop=='maize1')|(crop=='maize'):
                growMon = [np.array([9,10,11,12,13,14,15])] #
                yldSt = [1978]
                yldEnd = [2013]
                breakPt = 12
            elif crop=='maize2':
                growMon = [np.array([13,14,15,16])] #
                yldSt =  [1978]
                yldEnd = [2013]
                breakPt = 0
            elif (crop=='wheat')|(crop=='springwheat'):
                growMon = [np.array([9,10,11])]#8,12
                yldSt = [1978]
                yldEnd = [2013]
                breakPt = 7
                
        elif country is 'MX':
            cntName = 'Mexico'
            prodFactor = 1000. #convert tonnes to kg:
            yldFactor = 1000. #convert tonnes/ha to kg/ha             
            dataFile = 'MX_SAGARPA'
            yldSt = [1980]
            yldEnd =[2012]
            growMon = [np.array([5,6,7])]
            breakPt = 0 #number of years in to put the breakpoint for detrending   
           
        elif country is 'CA':
            cntName = 'Canada'
            dataFile = 'CA_CANSIM'
            prodFactor = 1000. #convert tonnes to kg
            yldFactor = 1000. #convert tonnes/ha to kg/ha   
            yldSt = [1950]
            yldEnd =[2013] 
            growMon = [np.array([5,6,7])]
            breakPt = 30 #number of years in to put the breakpoint for detrending   
        
        elif country is 'AU':
            cntName = 'Australia'
            dataFile = 'AUS_ABS'
            prodFactor = 1000. #convert tonnes to kg:
            yldFactor = 1000. #convert tonnes/ha to kg/ha      
            yldSt = [1950]
            yldEnd =[2012]
            growMon = [np.array([9,10,11])]
            breakPt = 30 #number of years in to put the breakpoint for detrending   

        elif country is 'CHN':
            cntName = 'China'
            dataFile = 'CHN'
            yldSt = [1949]
            yldEnd =[2013]
            prodFactor = 10000*1000. #10k tonnes to kg
            HAfactor = 666.7 #10k "management units" (mu) to ha
            yldFactor = (10000/666.66) # kg/mu to kg/ha
            if (crop=='wheat')|(crop=='maize')|(crop=='winterwheat'): 
                growMon = [np.array([4,5,6])]
            elif crop=='soy':
                growMon = [np.array([7,8,9])]
            breakPt = 30 #number of years in to put the breakpoint for detrending        
    
        elif country is 'IND':
            cntName = 'India'
            prodFactor = 1000.*1000. #convert 000's of tonnes to kg:
            yldFactor = 1. #kg/ha already
            HAfactor = 1000. #000's of ha to ha
            dataFile = 'IND'
            yldSt = [1977]
            yldEnd =[2012]
            growMon = [np.array([2,3,4])]
            breakPt = 0 #number of years in to put the breakpoint for detrending               

        elif country is 'EU':
            cntName = 'European Union'
            dataFile = 'FAOSTAT/EU'
            if(dataFile=='EU'):
                prodFactor = 1000.*1000. #convert 000's of tonnes to kg:
                yldFactor = 1000. #000's of tonnes /000's of HA to  kg/ha
                HAfactor = 1000. #000's of ha to ha
            if(dataFile=='FAOSTAT/EU'):                
                prodFactor = 1000. #convert tonnes to kg:
                yldFactor = .1 #hg/ha to kg/ha
                HAfactor = 1. #already in ha
            prcp_lag = ['Y'] #whether to include lagged precip in the regression
            yldSt = [1960]
            yldEnd =[2013] 
            growMon = [np.array([5,6,7])]
            breakPt = 20 #number of years in to put the breakpoint for detrending  

 
        #Read in admin codes and parse. This avoids confusing states like 'Distrito Federal' between countries
        adminCodes = pd.read_table(path+'Data/CropCalendar/MIRCA2000/unit_name_state.txt',sep='\t')
        #re-parse the unit code names
        for n in range(0,np.size(adminCodes['UnitName'])):
            if np.size(adminCodes['UnitName'][n].split('_'))>1:
                newName = adminCodes['UnitName'][n].split('_')[1]
                adminCodes['UnitName'][n]=newName

        #Calculate the climate anomalies for each subnational unit  
        ##################### sumP #####################
        sump = sumpNC.variables['sumP'][:]
        sump = sump - np.nanmean(sump,0)
        sumpYrs = sumpNC.variables['year'][:]
        sump = np.append(sump[:,:,180:],sump[:,:,:180],axis=2)
        sump=sump*cMsk
         
        ##################### kdd #####################
        kdd = kddNC.variables['KDD'][:]
        kdd = kdd - np.nanmean(kdd,0)
        kdd = kdd[:,::-1,:]
        kddYrs = kddNC.variables['year'][:] 
        kdd=kdd*cMsk
        
        ##################### SoilM #####################
        sm = smNC.variables['SoilM'][:]
        nans = np.zeros([sm.shape[0],30,360])*np.nan
        sm = np.append(sm[:,::-1,:],nans,axis=1)
        sm = sm - np.nanmean(sm,0)
        smYrs = smNC.variables['year'][:]
        sm=sm*cMsk
        
        #read in the yield data
        CNTRYyld = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/Yield/HistoricalData/'+dataFile+'/'+country+'_'+crop+'_'+stORreg+'.csv',
                                      header=0); #CNTRYyld = CNTRYyld.sort()
        #read in the harvested area data
        CNTRYha = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/HarvestedArea/HistoricalData/'+dataFile+'/'+country+'_'+crop+'_'+stORreg+'.csv',
                                      header=0); #CNTRYha = CNTRYha.sort()
        #read in the total production data
        CNTRYp = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/production/HistoricalData/'+dataFile+'/'+country+'_'+crop+'_'+stORreg+'.csv',
                                      header=0) #CNTRYp = CNTRYp.sort();
        #NaNs coded as 0 in the data
        CNTRYyld[CNTRYyld==0]=np.nan
        CNTRYha[CNTRYha==0]=np.nan
        
        states = CNTRYyld.columns
        #create the pandas dataframe to fill, use a multi index of [year, state]
        # To get all state/year combinations 
        index = pd.MultiIndex.from_product([yearNC,states], names=['year', 'state'])
        df = pd.DataFrame(index=index, columns=colNams)    
        
        CCC_file = netCDF4.Dataset('/Volumes/Data_Archive/Data/CropCalendar/Sacks/1.0Deg/'+ccCrop+'.crop.calendar.fill.nc')
        CC_tab = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/CropCalendar/Sacks/All_data_with_climate.csv')
        CCC = CCC_file['harvest'][:]
        
        #Region
        regions = np.load('/Volumes/Data_Archive/Data/SREX_shapefiles/numpy objects/refRegionsMod5.npy')
    
        for stNam in states:
            stYield = CNTRYyld[str(yldSt[0]):str(yldEnd[0])][stNam]
            HA = np.array(CNTRYha[str(yldSt[0]):str(yldEnd[0])][stNam])
            Prod = np.array(CNTRYp[str(yldSt[0]):str(yldEnd[0])][stNam])
            Prod=Prod[Prod!=0];
            Prod=Prod[~np.isnan(Prod)]; HA = HA[~np.isnan(HA)]; stYield = stYield.dropna();
            stYield = stYield*yldFactor
            Prod = Prod*prodFactor
            HA = HA*HAfactor 
            if (len(stYield) <= breakPt)|(len(Prod) <= breakPt)|(len(stYield)<9):
                continue
            if np.nanmax(100*(stYield.values*HA - Prod)/Prod)>2.:
                print( stNam)
                print( 'something is wrong, check units of yield, HA and production, avg percent error:')
                print( np.nanmean(100*(stYield.values*HA - Prod)/Prod))
                    #break
                    #if the difference in yield*HA and production is more than 2%, then something is wrong
            stYrs = stYield.index.year
            if ((country=='AR')|(country=='BR')):#Correct planting year to be harvest year if country is Argentina
                if (crop is 'maize')|(crop is 'soy'):
                    print('WARNING: reported date (planting) moved forward one year to harvest')
                    stYrs = stYield.index.year+1
            stYldAbs = np.array(signal.detrend(stYield,bp=breakPt))
            stYld = np.array(signal.detrend(stYield,bp=breakPt)/
                    (stYield-signal.detrend(stYield,bp=breakPt))) #%yield anom = (yld obs - yld exp)/yld exp
            expYld = (stYield-signal.detrend(stYield,bp=breakPt))
            stYldFD = np.array(stYield[1:])-np.array(stYield[:-1]) 
            stYldSm9 = np.array(stYield[4:-4]-moving_average(np.array(stYield),4))/moving_average(np.array(stYield),4)
            stYldSm5 = np.array(stYield[2:-2]-moving_average(np.array(stYield),2))/moving_average(np.array(stYield),2)   
            stYldSm9Abs= np.array(stYield[4:-4]-moving_average(np.array(stYield),4))
            stYldSm5Abs = np.array(stYield[2:-2]-moving_average(np.array(stYield),2))
            stYldGau = (stYield.values-ndimage.filters.gaussian_filter1d(stYield.values,3))/ndimage.filters.gaussian_filter1d(stYield.values,3)
            select = list(zip(stYrs,np.repeat(stNam,stYrs.size)))
            
            df.loc[select,'country']=np.repeat(country,stYrs.size)
            df.loc[select,'country name']=np.repeat(cntName,stYrs.size)
            df.loc[select,'yield']=np.array(stYield)
            df.loc[select,'Production']= Prod
            df.loc[select,'directProdAnom']= (Prod - ndimage.filters.gaussian_filter1d(Prod,3))
            df.loc[select,'directProdExp']= ndimage.filters.gaussian_filter1d(Prod,3)
            df.loc[select,'ProdAnomAbs']= stYldAbs*HA
            df.loc[select,'prodAnomGau']= (stYield.values-ndimage.filters.gaussian_filter1d(stYield.values,3))*HA
            df.loc[select[2:-2],'prodAnomSmooth5']= stYldSm5Abs*HA[2:-2]
            df.loc[select[4:-4],'prodAnomSmooth9']= stYldSm9Abs*HA[4:-4]
            df.loc[select,'yldAnom']=stYld
            df.loc[select,'yldAnomGau']=stYldGau
            df.loc[select[2:-2],'yldAnomSmooth5']=stYldSm5
            df.loc[select[4:-4],'yldAnomSmooth9']=stYldSm9
            df.loc[select,'expectedProd']=np.array(expYld)*HA
            df.loc[select,'expectedProdGau']=ndimage.filters.gaussian_filter1d(stYield.values,3)*HA
            df.loc[select[2:-2],'expectedProdSmooth5']=moving_average(np.array(stYield),2)*HA[2:-2]
            df.loc[select[4:-4],'expectedProdSmooth9']=moving_average(np.array(stYield),4)*HA[4:-4]
            df.loc[select,'expectedYld']=np.array(expYld)
            df.loc[select,'expectedYldGau']=ndimage.filters.gaussian_filter1d(stYield.values,3)
            df.loc[select[2:-2],'expectedYldSmooth5']=moving_average(np.array(stYield),2)
            df.loc[select[4:-4],'expectedYldSmooth9']=moving_average(np.array(stYield),4)
            df.loc[select,'yldAnomAbs']=stYldAbs
            df.loc[select,'yldAnomGauAbs']=(stYield.values-ndimage.filters.gaussian_filter1d(stYield.values,3))
            df.loc[select[2:-2],'yldAnomSmooth5Abs']=stYldSm5Abs
            df.loc[select[4:-4],'yldAnomSmooth9Abs']=stYldSm9Abs
            if np.shape(select)[0]>1:
                df.loc[select[1:],'yldAnomFD']=stYldFD
            df.loc[select,'Harvested_Area']=HA  
            
            #First go through to pull data on crop harvest days  
            #limit the analysis to land areas, so that when you sum of subnational areas you don't get values over the ocean                   
            MircaCodes = MircaCodes*LSmsk
            MircaCodesTemp = copy.deepcopy(MircaCodes)
            if (country is 'BR')|(country is 'US'):#Do we need title case?
                code = adminCodes['UnitCode'][adminCodes['UnitName']==stNam.title()]
            else:
                code = adminCodes['UnitCode'][adminCodes['UnitName']==stNam]
            if code.size==0:
                mismatch.append(stNam)
                continue
            for num in code:
                MircaCodesTemp[MircaCodes==num]=1
                MircaCodesTemp[MircaCodesTemp!=1]=0
                #Finally, calculate the julian harvest day
            if np.sum(MircaCodesTemp>0):
                har = stats.mode(CCC[MircaCodesTemp==1],nan_policy='omit') [0]
                reg = stats.mode(regions[MircaCodesTemp==1])[0][0]
                smState = np.squeeze(np.nanmean(sm[:,MircaCodesTemp==1],axis=1))
                sumpState = np.squeeze(np.nanmean(sump[:,MircaCodesTemp==1],axis=1))
                kddState = np.squeeze(np.nanmean(kdd[:,MircaCodesTemp==1],axis=1))
            else: 
                har = np.nan;
                reg=np.nan
                smState = np.repeat(np.nan,smYrs.shape[0])
                sumpState = np.repeat(np.nan,sumpYrs.shape[0])
                kddState = np.repeat(np.nan,kddYrs.shape[0])
            
            df.loc[list(zip(smYrs,np.repeat(stNam,smYrs.size))),'SoilM']=smState
            df.loc[list(zip(sumpYrs,np.repeat(stNam,sumpYrs.size))),'sumP']=sumpState
            df.loc[list(zip(kddYrs,np.repeat(stNam,kddYrs.size))),'KDD']=kddState                          
            df.loc[select,'flwrSeas']=np.repeat(growMon[0][1],stYrs.size)
            df.loc[select,'harDay']=np.repeat(har,stYrs.size)
            df.loc[select,'region']=np.repeat(reg,stYrs.size)
          
            
            if(((country=='AR')&(crop=='maize'))|((country=='BR')&(crop=='maize'))|((country=='AR')&(crop=='soy'))|((country=='BR')&(crop=='soy'))):
                ENselect = list(zip(np.array(ENyrs)+1,np.repeat(stNam,np.size(ENyrs))))
                ENselect_1 = list(zip(np.array(ENyrs)+1-1,np.repeat(stNam,np.size(ENyrs))))
                ENselect1 = list(zip(np.array(ENyrs)+1+1,np.repeat(stNam,np.size(ENyrs))))
                LNselect = list(zip(np.array(LNyrs)+1,np.repeat(stNam,np.size(LNyrs))))
                LNselect_1 = list(zip(np.array(LNyrs)+1-1,np.repeat(stNam,np.size(LNyrs))))
                LNselect1 = list(zip(np.array(LNyrs)+1+1,np.repeat(stNam,np.size(LNyrs))))
                LNselect2 = list(zip(np.array(LNyrs)+1+2,np.repeat(stNam,np.size(LNyrs))))
            else:
                ENselect = list(zip(np.array(ENyrs),np.repeat(stNam,np.size(ENyrs))))
                ENselect_1 = list(zip(np.array(ENyrs)-1,np.repeat(stNam,np.size(ENyrs))))
                ENselect1 = list(zip(np.array(ENyrs)+1,np.repeat(stNam,np.size(ENyrs))))
                LNselect = list(zip(np.array(LNyrs),np.repeat(stNam,np.size(LNyrs))))
                LNselect_1 = list(zip(np.array(LNyrs)-1,np.repeat(stNam,np.size(LNyrs))))
                LNselect1 = list(zip(np.array(LNyrs)+1,np.repeat(stNam,np.size(LNyrs))))
                LNselect2 = list(zip(np.array(LNyrs)+2,np.repeat(stNam,np.size(LNyrs))))                
            df.loc[ENselect,'EN']=0;df.loc[ENselect_1,'EN']=-1;df.loc[ENselect1,'EN']=1
            df.loc[LNselect,'LN']=0;df.loc[LNselect_1,'LN']=-1;
            df.loc[LNselect,'LN2']=0;df.loc[LNselect_1,'LN2']=-1;#repeat for when LN2 == LN-1 (so it is overwritten)
            df.loc[LNselect1,'LN']=1
            df.loc[LNselect1,'LN2']=1
            df.loc[LNselect2,'LN2']=2
        df.to_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/byCountry/'+country+'_'+crop+dt+notes+'.csv',',')
            
elapsed = time.clock()-start
print( elapsed)
print('mismatched state codes:')
print(mismatch)
     