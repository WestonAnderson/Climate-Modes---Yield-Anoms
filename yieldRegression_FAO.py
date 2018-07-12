# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 13:20:47 2015
@author: westonanderson

This script calculates yield anomalies from FAO country-level data.
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

crops = ['Wheat'] #['Soybeans','Maize','Wheat']
climCropNames = {'Soybeans':'soy','Maize':'maize','Wheat':'winterwheat'}
sacksCropNames = {'Soybeans':'Soybeans','Maize':'Maize','Wheat':'Wheat.Main'}
dt = '' # '_dt' or ''
thresh = 0.001 #0.005 = 0.5%

ENyrs = [1963, 1965, 1968, 1972, 1976, 1979, 1982, 1986, 1991, 1994, 1997, 2002, 2004, 2006, 2009] #EN 1953, 1957, 
LNyrs = [1964, 1970, 1973, 1983, 1988, 1995, 1998, 2007, 2010] #LN0_2 1954, 

colNams = [ 'country name','region','flwrSeas','yield','Production','Harvested_Area',\
            'yldAnom','yldAnomFD','yldAnomGau','yldAnomSmooth5','yldAnomSmooth9',\
            'yldAnomLinAbs','yldAnomSmooth5Abs','yldAnomSmooth9Abs','yldAnomGauAbs',\
            'expectedYld','expectedYldSmooth5','expectedYldSmooth9','expectedYldGau',\
            'prodAnomSmooth5','prodAnomSmooth9','prodAnomGau','directProdAnom','directProdExp',\
            'expectedProd','expectedProdSmooth5','expectedProdSmooth9','expectedProdGau',\
            'GDD','KDD','KDD_sum','sumPan','sumPan_sum','sumP','sumP2','sumPanLag','PET','P_E',\
            'VPD','soilM','soilMseed','EN','LN','LN2']
    
mismatch = []

years = np.array(range(1900,2016))
stYr = [1950];endYr = [2010] #start and end years for the climate analysis
yldFactor = .1; #FAO yield data is in Hg/HA. Convert to Kg/HA
prdFactor = 1000.; #FAO data is in tonnes, convert to kg
HAfactor = 1. #FAO data already in ha

FAO = pd.DataFrame.from_csv('/Volumes/Data_Archive/Data/Yield/HistoricalData/FAOSTAT/FAOSTAT_2014.csv',header=0)
MircaCodes = np.load('/Volumes/Data_Archive/Data/CropCalendar/CCCmask/1.0Deg/UnitCode/SubnatUnitCode.npy')
CntryCode = np.load('/Volumes/Data_Archive/Data/CropCalendar/CCCmask/UnitCode/CntryUnitCode.npy', encoding='latin1')
adminCodes = pd.read_table('/Volumes/Data_Archive/Data/CropCalendar/MIRCA2000/unit_name_FAO.txt',sep='\t')

#Rename FAO countries to match other sources
FAO['Area'][FAO['Area']=='Bolivia (Plurinational State of)']='Bolivia'
FAO['Area'][FAO['Area']=='Venezuela (Bolivarian Republic of)']='Venezuela'
FAO['Area'][FAO['Area']=='Ethiopia PDR']='Ethiopia'
FAO['Area'][FAO['Area']=='Democratic Republic of the Congo']="Dem. Rep. Congo"
FAO['Area'][FAO['Area']=='Dominican Republic']="Dominican Rep."
FAO['Area'][FAO['Area']=='Central African Republic']='Central African Rep.'
FAO['Area'][FAO['Area']=='Sudan (former)']='Sudan'
FAO['Area'][FAO['Area']=='United Republic of Tanzania']='Tanzania'
FAO['Area'][FAO['Area']=='Republic of Korea']='South Korea'
FAO['Area'][FAO['Area']=="Democratic People's Republic of Korea"]='North Korea'
FAO['Area'][FAO['Area']=='Viet Nam']='Vietnam'
FAO['Area'][FAO['Area']=='Iran (Islamic Republic of)']='Iran'
FAO['Area'][FAO['Area']=='Syrian Arab Republic']='Syria'
FAO['Area'][FAO['Area']=='Serbia']='Serbia and Montenegro'
FAO['Area'][FAO['Area']=='Yugoslav SFR']='Yugoslavia'
FAO['Area'][FAO['Area']=='Georgia']='Georgia (cntry)'

#Remove countries that will be included at the state level

#re-parse the unit code names
for n in range(0,np.size(adminCodes['UnitName'])):
    if np.size(adminCodes['UnitName'][n].split('_'))>1:
        newName = adminCodes['UnitName'][n].split('_')[0]
        adminCodes['UnitName'][n]=newName
msk = netCDF4.Dataset('/Volumes/Data_Archive/Data/landSeaMask/1.0Deg/lsmask.nc')
LSmsk = msk.variables['mask'][0,...];
LSmsk = np.append(LSmsk[:,180:],LSmsk[:,:180],axis=1)

for crop in crops:
    #Read in the climate data
    smNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/globalCorrelatedRisks/IntegratedVars'+dt+'/preHarvest/SoilM'+climCropNames[crop]+'.nc')
    kddNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/globalCorrelatedRisks/IntegratedVars'+dt+'/preHarvest/KDD'+climCropNames[crop]+'.nc')
    sumpNC = netCDF4.Dataset('/Volumes/Data_Archive/Results/globalCorrelatedRisks/IntegratedVars'+dt+'/preHarvest/sumP'+climCropNames[crop]+'.nc')

    #<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>##<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#
    # calculate the crop mask
    lats = sumpNC.variables['latitude'][:]
    lons = sumpNC.variables['longitude'][:]
    cMsk = np.load(str('/Volumes/Data_Archive/Data/CropMaps_ComonGrid/AllArea/HarvestedArea/gpcc_grid_ensemble/ensemble_'+climCropNames[crop][0]+'.npy'))
    cLats = np.append((lats[0]-lats[1])+lats[0], lats)
    cLats = cLats +(lats[1]-lats[0])/2
    dlon = (lons[1]-lons[0])*(np.pi/180.) #difference in longitude (in radians) remains constant
    R = 6371000. #Radius of the earth in m
    cellArea = 2*np.pi*(R**2)*np.abs(np.sin(cLats[1:]*(np.pi/180.))-
                               np.sin(cLats[:-1]*(np.pi/180.)))*dlon #in m
    cellArea = np.repeat(cellArea[:,np.newaxis],lons.shape,axis=1)
    haThresh = (cellArea*thresh)*1./10000 #converting the threshhold to hectares
    cMsk[cMsk<haThresh]=0;cMsk[cMsk>=haThresh]=1
 
    #Calculate the climate anomalies for each subnational unit  
    ##################### sumP #####################
    sump = sumpNC.variables['sumP'][:]
    sump = sump - np.nanmean(sump,0)
    sumpYrs = sumpNC.variables['year'][:]
    sump = np.append(sump[:,:,180:],sump[:,:,:180],axis=2)
    sump=np.array(sump*cMsk)
     
    ##################### kdd #####################
    kdd = kddNC.variables['KDD'][:]
    kdd = kdd - np.nanmean(kdd,0)
    kdd = kdd[:,::-1,:]
    kddYrs = kddNC.variables['year'][:] 
    kdd=np.array(kdd*cMsk)
    
    ##################### SoilM #####################
    sm = smNC.variables['SoilM'][:]
    nans = np.zeros([sm.shape[0],30,360])*np.nan
    sm = np.append(sm[:,::-1,:],nans,axis=1)
    sm = sm - np.nanmean(sm,0)
    smYrs = smNC.variables['year'][:]
    sm=np.array(sm*cMsk)
    
    CNTRYp = FAO[(FAO.Element=='Production')&(FAO.Item==crop)][['Area','Year','Value']]
    CNTRYyld = FAO[(FAO.Element=='Yield')&(FAO.Item==crop)][['Area','Year','Value']]
    CNTRYha = FAO[(FAO.Element=='Area harvested')&(FAO.Item==crop)][['Area','Year','Value']]
    CNTRYyld.set_index('Year',inplace=True);CNTRYha.set_index('Year',inplace=True);CNTRYp.set_index('Year',inplace=True)

    #Combine countries
    ########################################## USSR  ##########################################
    USSRp = CNTRYp[(CNTRYp.Area=='Azerbaijan')|(CNTRYp.Area=='Armenia')|(CNTRYp.Area=='Belarus')|(CNTRYp.Area=='Estonia')|
            (CNTRYp.Area=='Kazakhstan')|(CNTRYp.Area=='Kyrgyzstan')|(CNTRYp.Area=='Latvia')|(CNTRYp.Area=='Georgia (cntry)')|
            (CNTRYp.Area=='Lithuania')|(CNTRYp.Area=='Republic of Moldova')|(CNTRYp.Area=='Russian Federation')|
            (CNTRYp.Area=='Tajikistan')|(CNTRYp.Area=='Ukraine')|(CNTRYp.Area=='Uzbekistan')].sum(level='Year').Value
            
    USSRha = CNTRYha[(CNTRYha.Area=='Azerbaijan')|(CNTRYha.Area=='Armenia')|(CNTRYha.Area=='Belarus')|(CNTRYha.Area=='Estonia')|
            (CNTRYha.Area=='Kazakhstan')|(CNTRYha.Area=='Kyrgyzstan')|(CNTRYha.Area=='Latvia')|(CNTRYha.Area=='Georgia (cntry)')|
            (CNTRYha.Area=='Lithuania')|(CNTRYha.Area=='Republic of Moldova')|(CNTRYha.Area=='Russian Federation')|
            (CNTRYha.Area=='Tajikistan')|(CNTRYha.Area=='Ukraine')|(CNTRYha.Area=='Uzbekistan')].sum(level='Year').Value
    
    USSRpLst = pd.DataFrame(data={'Year':np.arange(1992,2015),'Area':np.repeat('USSR',np.arange(1992,2015).size),'Value':USSRp})
    USSRhaLst = pd.DataFrame(data={'Year':np.arange(1992,2015),'Area':np.repeat('USSR',np.arange(1992,2015).size),'Value':USSRha})
    USSRyldLst = pd.DataFrame(data={'Year':np.arange(1992,2015),'Area':np.repeat('USSR',np.arange(1992,2015).size),'Value':USSRp*10000/USSRha}) #factor of 10000 is to convert from tonnes to Hg to keep units consistent

    USSRpLst.set_index('Year',inplace=True)
    USSRhaLst.set_index('Year',inplace=True)
    USSRyldLst.set_index('Year',inplace=True)

    CNTRYp=CNTRYp.append(USSRpLst)
    CNTRYha=CNTRYha.append(USSRhaLst)
    CNTRYyld=CNTRYyld.append(USSRyldLst)
 
    
    ########################################## Yugoslavia  ##########################################
    YugoP = CNTRYp[(CNTRYp.Area=='Croatia')|(CNTRYp.Area=='Slovenia')|(CNTRYp.Area=='Bosnia and Herzegovina')|
            (CNTRYp.Area=='The former Yugoslav Republic of Macedonia')|(CNTRYp.Area=='Serbia and Montenegro')].sum(level='Year').Value

    YugoHA = CNTRYha[(CNTRYha.Area=='Croatia')|(CNTRYha.Area=='Slovenia')|(CNTRYha.Area=='Bosnia and Herzegovina')|
            (CNTRYha.Area=='The former Yugoslav Republic of Macedonia')|(CNTRYha.Area=='Serbia and Montenegro')].sum(level='Year').Value
  
    YugoPlst = pd.DataFrame(data={'Year':np.arange(1992,2015),'Area':np.repeat('Yugoslavia',np.arange(1992,2015).size),'Value':YugoP})
    YugoHAlst = pd.DataFrame(data={'Year':np.arange(1992,2015),'Area':np.repeat('Yugoslavia',np.arange(1992,2015).size),'Value':YugoHA})
    YugoYLDlst = pd.DataFrame(data={'Year':np.arange(1992,2015),'Area':np.repeat('Yugoslavia',np.arange(1992,2015).size),'Value':YugoP*10000/YugoHA})#factor of 10000 is to convert from tonnes to Hg to keep units consistent

    YugoPlst.set_index('Year',inplace=True)
    YugoHAlst.set_index('Year',inplace=True)
    YugoYLDlst.set_index('Year',inplace=True)    

    CNTRYp=CNTRYp.append(YugoPlst)
    CNTRYha=CNTRYha.append(YugoHAlst)
    CNTRYyld=CNTRYyld.append(YugoYLDlst)

    ########################################## Czechoslovakia  ##########################################
    CzechP = CNTRYp[(CNTRYp.Area=='Slovakia')|(CNTRYp.Area=='Czech Republic')|(CNTRYp.Area=='Czechia')].sum(level='Year').Value
    CzechHA = CNTRYha[(CNTRYha.Area=='Slovakia')|(CNTRYha.Area=='Czech Republic')|(CNTRYha.Area=='Czechia')].sum(level='Year').Value
  
    CzechPlst = pd.DataFrame(data={'Year':np.arange(1993,2015),'Area':np.repeat('Czechoslovakia',np.arange(1993,2015).size),'Value':CzechP})
    CzechHAlst = pd.DataFrame(data={'Year':np.arange(1993,2015),'Area':np.repeat('Czechoslovakia',np.arange(1993,2015).size),'Value':CzechHA})
    CzechYLDlst = pd.DataFrame(data={'Year':np.arange(1993,2015),'Area':np.repeat('Czechoslovakia',np.arange(1993,2015).size),'Value':CzechP*10000/CzechHA})#factor of 10000 is to convert from tonnes to Hg to keep units consistent

    CzechPlst.set_index('Year',inplace=True)
    CzechHAlst.set_index('Year',inplace=True)
    CzechYLDlst.set_index('Year',inplace=True)    

    CNTRYp=CNTRYp.append(CzechPlst)
    CNTRYha=CNTRYha.append(CzechHAlst)
    CNTRYyld=CNTRYyld.append(CzechYLDlst)
    
    #remove the constituent countries from combined states and from all countries that will be included at the state level from separate data sources
    CNTRYp=CNTRYp[~(CNTRYp.Area=='China')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='China, mainland')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='China, Taiwan Province of')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='China, Hong Kong SAR')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='United States of America')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Argentina')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Brazil')]
    CNTRYha=CNTRYha[~(CNTRYha.Area=='China, mainland')]
    CNTRYha=CNTRYha[~(CNTRYha.Area=='China, Taiwan Province of')]
    CNTRYha=CNTRYha[~(CNTRYha.Area=='China, Hong Kong SAR')]
    CNTRYha=CNTRYha[~(CNTRYha.Area=='United States of America')]
    CNTRYha=CNTRYha[~(CNTRYha.Area=='Argentina')]
    CNTRYha=CNTRYha[~(CNTRYha.Area=='Brazil')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='China, mainland')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='China, Taiwan Province of')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='China, Hong Kong SAR')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='United States of America')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Argentina')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Brazil')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Azerbaijan')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Armenia')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Belarus')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Estonia')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Kazakhstan')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Kyrgyzstan')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Latvia')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Georgia (cntry)')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Lithuania')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Republic of Moldova')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Russian Federation')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Tajikistan')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Ukraine')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Uzbekistan')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Croatia')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Slovenia')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Bosnia and Herzegovina')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='The former Yugoslav Republic of Macedonia')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Serbia and Montenegro')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Slovakia')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Czechia')]
    CNTRYp=CNTRYp[~(CNTRYp.Area=='Czech Republic')]
    if (crop=='Maize'):
        CNTRYp=CNTRYp[~(CNTRYp.Area=='Mexico')]
        CNTRYyld=CNTRYyld[~(CNTRYyld.Area=='Mexico')]
        CNTRYha=CNTRYha[~(CNTRYha.Area=='Mexico')]
    if (crop=='Wheat'):  
        CNTRYp=CNTRYp[~(CNTRYp.Area=='Australia')]
        CNTRYp=CNTRYp[~(CNTRYp.Area=='Canada')]
        CNTRYha=CNTRYha[~(CNTRYha.Area=='Australia')]
        CNTRYha=CNTRYha[~(CNTRYha.Area=='Canada')]
        CNTRYyld=CNTRYyld[~(CNTRYyld.Area=='Australia')]
        CNTRYyld=CNTRYyld[~(CNTRYyld.Area=='Canada')]
    if (crop=='Maize')|(crop=='Wheat'):
        ########################################## Belgium-Luxembourg  ##########################################
        blP = CNTRYp[CNTRYp.Area=='Belgium'].Value+CNTRYp[CNTRYp.Area=='Luxembourg'].Value
        blHA = CNTRYha[CNTRYha.Area=='Belgium'].Value+CNTRYha[CNTRYha.Area=='Luxembourg'].Value
    
        blPlst = pd.DataFrame(data={'Year':np.arange(2000,2015),'Area':np.repeat('Belgium-Luxembourg',np.arange(2000,2015).size),'Value':blP})
        blHAlst = pd.DataFrame(data={'Year':np.arange(2000,2015),'Area':np.repeat('Belgium-Luxembourg',np.arange(2000,2015).size),'Value':blHA})
        blYLDlst = pd.DataFrame(data={'Year':np.arange(2000,2015),'Area':np.repeat('Belgium-Luxembourg',np.arange(2000,2015).size),'Value':blP*10000/blHA})#factor of 10000 is to convert from tonnes to Hg to keep units consistent
    
        blPlst.set_index('Year',inplace=True)
        blHAlst.set_index('Year',inplace=True)
        blYLDlst.set_index('Year',inplace=True)    
    
        CNTRYp=CNTRYp.append(blPlst)
        CNTRYha=CNTRYha.append(blHAlst)
        CNTRYyld=CNTRYyld.append(blYLDlst)   

        CNTRYp=CNTRYp[~(CNTRYp.Area=='India')]
        CNTRYyld=CNTRYyld[~(CNTRYyld.Area=='India')] 
        CNTRYha=CNTRYha[~(CNTRYha.Area=='India')]         
        CNTRYp=CNTRYp[~(CNTRYp.Area=='Belgium')]
        CNTRYp=CNTRYp[~(CNTRYp.Area=='Luxembourg')]
 
    #create the pandas dataframe to fill, use a multi index of [year, state]
    # To get all state/year combinations 
    cntries = np.unique(CNTRYp.Area)
    index = pd.MultiIndex.from_product([years,cntries], names=['year', 'country'])
    df = pd.DataFrame(index=index, columns=colNams)    
    
    #linearly detrend each state yields, then 
    for cntNam in cntries:
        cntYield = CNTRYyld[CNTRYyld.Area==cntNam].Value
        if cntYield.size<10:
            #print(cntNam)
            continue

        #now calculate yield anomalies    
        stYrs = cntYield.index
        HA = np.array(CNTRYha[CNTRYha.Area==cntNam].loc[stYrs].Value.dropna())
        Prod = np.array(CNTRYp[CNTRYp.Area==cntNam].loc[stYrs].Value.dropna())
        cntYield = cntYield.dropna();
        HA = HA[Prod!=0];Prod=Prod[Prod!=0];
        cntYield = cntYield*yldFactor; Prod = Prod*prdFactor; HA = HA*HAfactor 

        stYldAbs = np.array(signal.detrend(cntYield,bp=0))
        stYld = np.array(signal.detrend(cntYield,bp=0)/
                (cntYield-signal.detrend(cntYield,bp=0))) #%yield anom = (yld obs - yld exp)/yld exp
        expYld = (cntYield-signal.detrend(cntYield,bp=0))
        stYldFD = np.array(cntYield[1:])-np.array(cntYield[:-1]) 
        stYldSm9 = np.array(cntYield[4:-4]-moving_average(np.array(cntYield),4))/moving_average(np.array(cntYield),4)
        stYldSm5 = np.array(cntYield[2:-2]-moving_average(np.array(cntYield),2))/moving_average(np.array(cntYield),2)   
        stYldSm9Abs= np.array(cntYield[4:-4]-moving_average(np.array(cntYield),4))
        stYldSm5Abs = np.array(cntYield[2:-2]-moving_average(np.array(cntYield),2))
        stYldGau = (cntYield.values-ndimage.filters.gaussian_filter1d(cntYield.values,3))/ndimage.filters.gaussian_filter1d(cntYield.values,3)
        select = list(zip(stYrs,np.repeat(cntNam,stYrs.size)))

        df.loc[select,'country']=np.repeat(np.nan,stYrs.size)
        df.loc[select,'state']=np.repeat(np.nan,stYrs.size)
        df.loc[select,'country name']=np.repeat(cntNam,stYrs.size)
        df.loc[select,'yield']=np.array(cntYield)
        df.loc[select,'Production']= Prod
        df.loc[select,'directProdAnom']= (Prod - ndimage.filters.gaussian_filter1d(Prod,3))
        df.loc[select,'directProdExp']= ndimage.filters.gaussian_filter1d(Prod,3)
        df.loc[select,'ProdAnomLinAbs']= stYldAbs*HA
        df.loc[select,'prodAnomGau']= (cntYield.values-ndimage.filters.gaussian_filter1d(cntYield.values,3))*HA
        df.loc[select[2:-2],'prodAnomSmooth5']= stYldSm5Abs*HA[2:-2]
        df.loc[select[4:-4],'prodAnomSmooth9']= stYldSm9Abs*HA[4:-4]
        df.loc[select,'yldAnom']=stYld
        df.loc[select,'yldAnomGau']=stYldGau
        df.loc[select[2:-2],'yldAnomSmooth5']=stYldSm5
        df.loc[select[4:-4],'yldAnomSmooth9']=stYldSm9
        df.loc[select,'expectedProd']=np.array(expYld)*HA
        df.loc[select,'expectedProdGau']=ndimage.filters.gaussian_filter1d(cntYield.values,3)*HA
        df.loc[select[2:-2],'expectedProdSmooth5']=moving_average(np.array(cntYield),2)*HA[2:-2]
        df.loc[select[4:-4],'expectedProdSmooth9']=moving_average(np.array(cntYield),4)*HA[4:-4]
        df.loc[select,'expectedYld']=np.array(expYld)
        df.loc[select,'expectedYldGau']=ndimage.filters.gaussian_filter1d(cntYield.values,3)
        df.loc[select[2:-2],'expectedYldSmooth5']=moving_average(np.array(cntYield),2)
        df.loc[select[4:-4],'expectedYldSmooth9']=moving_average(np.array(cntYield),4)
        df.loc[select,'yldAnomLinAbs']=stYldAbs
        df.loc[select,'yldAnomGauAbs']=(cntYield.values-ndimage.filters.gaussian_filter1d(cntYield.values,3))
        df.loc[select[2:-2],'yldAnomSmooth5Abs']=stYldSm5Abs
        df.loc[select[4:-4],'yldAnomSmooth9Abs']=stYldSm9Abs
        if np.shape(select)[0]>1:
            df.loc[select[1:],'yldAnomFD']=stYldFD
        df.loc[select,'Harvested_Area']=HA            
    
        #Region
        regions = np.load('/Volumes/Data_Archive/Data/SREX_shapefiles/numpy objects/refRegionsMod5.npy')
        
        #Pull data on crop harvest days
        CCC_file = netCDF4.Dataset('/Volumes/Data_Archive/Data/CropCalendar/Sacks/1.0Deg/'+sacksCropNames[crop]+'.crop.calendar.fill.nc')
        CCC = CCC_file['harvest'][:]
        
        #limit the analysis to land areas, so that when you sum of subnational areas you don't get values over the ocean                   
        MircaCodes = MircaCodes*LSmsk
        MircaCodesTemp = copy.deepcopy(MircaCodes)
        code = adminCodes['UnitCode'][adminCodes['UnitName']==cntNam]
        if code.size==0:
            mismatch.append(cntNam)
            continue
        for num in code:
            MircaCodesTemp[MircaCodes==num]=1
        MircaCodesTemp[MircaCodesTemp!=1]=0


        #Finally, calculate the julian harvest day
        har = np.nanmean(CCC[MircaCodesTemp==1])
        if np.sum(MircaCodesTemp==1)>0:
            reg = stats.mode(regions[MircaCodesTemp==1])[0][0]
            smState = np.squeeze(np.nanmean(sm[:,MircaCodesTemp==1],axis=1))
            sumpState = np.squeeze(np.nanmean(sump[:,MircaCodesTemp==1],axis=1))
            kddState = np.squeeze(np.nanmean(kdd[:,MircaCodesTemp==1],axis=1))
        else: 
            reg=np.nan
            smState = np.repeat(np.nan,smYrs.shape[0])
            sumpState = np.repeat(np.nan,sumpYrs.shape[0])
            kddState = np.repeat(np.nan,kddYrs.shape[0])
            
        
        df.loc[list(zip(smYrs,np.repeat(cntNam,smYrs.size))),'SoilM']=smState
        df.loc[list(zip(sumpYrs,np.repeat(cntNam,sumpYrs.size))),'sumP']=sumpState
        df.loc[list(zip(kddYrs,np.repeat(cntNam,kddYrs.size))),'KDD']=kddState 
        df.loc[select,'flwrSeas']=np.repeat(np.nan,stYrs.size)
        df.loc[select,'harDay']=np.repeat(har,stYrs.size)
        df.loc[select,'region']=np.repeat(reg,stYrs.size)
        ENselect = list(zip(np.array(ENyrs),np.repeat(cntNam,np.size(ENyrs))))
        ENselect_1 = list(zip(np.array(ENyrs)-1,np.repeat(cntNam,np.size(ENyrs))))
        ENselect1 = list(zip(np.array(ENyrs)+1,np.repeat(cntNam,np.size(ENyrs))))
        LNselect = list(zip(np.array(LNyrs),np.repeat(cntNam,np.size(LNyrs))))
        LNselect_1 = list(zip(np.array(LNyrs)-1,np.repeat(cntNam,np.size(LNyrs))))
        LNselect1 = list(zip(np.array(LNyrs)+1,np.repeat(cntNam,np.size(LNyrs))))
        LNselect2 = list(zip(np.array(LNyrs)+2,np.repeat(cntNam,np.size(LNyrs))))
        df.loc[ENselect,'EN']=0;df.loc[ENselect_1,'EN']=-1;df.loc[ENselect1,'EN']=1
        df.loc[LNselect,'LN']=0;df.loc[LNselect_1,'LN']=-1;
        df.loc[LNselect,'LN2']=0;df.loc[LNselect_1,'LN2']=-1;#repeat for when LN2 == LN-1 (so it is overwritten)
    if crop=='Soybeans':crop='Soy'           
    df.to_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/Dataframes/FAO_'+crop+'_'+str(years[0])+'_'+str(years[-1])+'.csv',',')     
elapsed = time.clock()-start
print( elapsed)
     