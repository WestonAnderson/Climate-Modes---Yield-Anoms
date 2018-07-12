#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 15:12:12 2018

@author: weston
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
import pandas as pd
import netCDF4
import numpy as np
import matplotlib
from unidecode import unidecode
import xarray as xr


exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/GADM_dict.py').read())
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/project/Climate Modes - Yield Anoms/ENSOgroups.py').read())

years = [1981,2011]

## IOD ## 
iod = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/SVD/'+
                             'all crops1961_2012_Jan-DecENSOyr0har_iodSVD1_janSVD_allTrops.csv')  ## TAV ## 
iod.reset_index(inplace=True);
iod.set_index(['state','year'],inplace=True,drop=True)
iod.rename(columns={'SVD':'IOD'},inplace=True)
## TAV ## 
tav = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/SVD/'+
                             'all crops1961_2012_Jan-DecENSOyr0har_tavSVD2_janSVD_allTrops.csv')  
tav.reset_index(inplace=True);
tav.set_index(['state','year'],inplace=True,drop=True)
tav.rename(columns={'SVD':'TAV'},inplace=True)
###### ENSO ######
enso = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/SVD/'+
                             'all crops1961_2012_Jan-DecENSOyr0har_SVD2_janSVD_allTrops.csv')   
enso.reset_index(inplace=True);
enso.set_index(['state','year'],inplace=True,drop=True)
enso.rename(columns={'SVD':'ENSO'},inplace=True)
###### NAO ######
nao = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/SVD/'+
                             'DJFwheat1961_2012_Jan-DecENSOyr0har_naoSVD1_janSVD.csv')   
nao.reset_index(inplace=True);
nao.set_index(['state','year'],inplace=True,drop=True)
nao.rename(columns={'SVD':'NAO'},inplace=True)

crDF = pd.concat([enso,nao['NAO'],iod['IOD'],tav['TAV']],axis=1)
crDF.reset_index(inplace=True);
crDF = crDF[(crDF.year>=years[0])&(crDF.year<=years[1])]
#crDF.set_index(['state','year'],inplace=True,drop=True)

Regions = ['Northeast Brazil','West Africa','Southwest Mexico']
color = ['cornflowerblue','maroon','gold','peru']

fig = plt.figure(); ax1 = plt.subplot(111);
for iR in range(len(Regions)):
    regName = Regions[iR]
    mRegs = ['maize_'+x.upper () for x in mRegions[regName]]
    reg = crDF[['Harvested_Area','yldAnomGau','ENSO','NAO','IOD','TAV','year']][crDF['state'].isin(mRegs)]
    reg['ProdAnomGau']=reg['Harvested_Area']*reg['yldAnomGau']
    reg['ensoProd']=reg['Harvested_Area']*reg['ENSO']
    reg['naoProd']=reg['Harvested_Area']*reg['NAO']
    reg['iodProd']=reg['Harvested_Area']*reg['IOD']
    reg['tavProd']=reg['Harvested_Area']*reg['TAV']      
    reg.set_index('year',drop=True,inplace=True)
    
    regYld = reg['ProdAnomGau'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    ensoYld =reg['ensoProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    naoYld =reg['naoProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    tavYld =reg['tavProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    iodYld =reg['iodProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    enYld =(reg['naoProd'].sum(level='year').fillna(0)+reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
    eniYld =(reg['iodProd'].sum(level='year').fillna(0)+reg['naoProd'].sum(level='year').fillna(0)+
              reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
    enitYld =(reg['tavProd'].sum(level='year').fillna(0)+reg['iodProd'].sum(level='year').fillna(0)+reg['naoProd'].sum(level='year').fillna(0)+
              reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
 
    ax1.plot(ensoYld.values,tavYld.values,'o',color=color[iR])
    print(np.corrcoef(ensoYld.values,tavYld.values)[0][1])
    
ax1.vlines(0,-.4,.4); ax1.hlines(0,-.4,.4)
ax1.set_ylim(-.4,.4);ax1.set_xlim(-.4,.4)
fig.set_size_inches(12,12);
plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/basin comparison/maize Atlantic Ocean.png')
plt.close()
    



Regions = ['Southeast Africa','East Africa']

fig = plt.figure(); ax1 = plt.subplot(111);
for iR in range(len(Regions)):
    regName = Regions[iR]
    mRegs = ['maize_'+x.upper () for x in mRegions[regName]]
    reg = crDF[['Harvested_Area','yldAnomGau','ENSO','NAO','IOD','TAV','year']][crDF['state'].isin(mRegs)]
    reg['ProdAnomGau']=reg['Harvested_Area']*reg['yldAnomGau']
    reg['ensoProd']=reg['Harvested_Area']*reg['ENSO']
    reg['naoProd']=reg['Harvested_Area']*reg['NAO']
    reg['iodProd']=reg['Harvested_Area']*reg['IOD']
    reg['tavProd']=reg['Harvested_Area']*reg['TAV']      
    reg.set_index('year',drop=True,inplace=True)
    
    regYld = reg['ProdAnomGau'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    ensoYld =reg['ensoProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    naoYld =reg['naoProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    tavYld =reg['tavProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    iodYld =reg['iodProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    enYld =(reg['naoProd'].sum(level='year').fillna(0)+reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
    eniYld =(reg['iodProd'].sum(level='year').fillna(0)+reg['naoProd'].sum(level='year').fillna(0)+
              reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
    enitYld =(reg['tavProd'].sum(level='year').fillna(0)+reg['iodProd'].sum(level='year').fillna(0)+reg['naoProd'].sum(level='year').fillna(0)+
              reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
 
    ax1.plot(ensoYld.values,iodYld.values,'o',color=color[iR])
    print(np.corrcoef(ensoYld.values,iodYld.values)[0][1])
    
ax1.vlines(0,-.4,.4); ax1.hlines(0,-.4,.4)
ax1.set_ylim(-.4,.4);ax1.set_xlim(-.4,.4)
fig.set_size_inches(12,12);
plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/basin comparison/maize Indian Ocean.png')
plt.close()
        
    
"""    
for regName in ['Northeast Brazil','West Africa','Southwest Mexico']:
    mRegs = ['maize_'+x.upper () for x in mRegions[regName]]
    reg = crDF[['Harvested_Area','yldAnomGau','ENSO','NAO','IOD','TAV','year']][crDF['state'].isin(mRegs)]
    reg['ProdAnomGau']=reg['Harvested_Area']*reg['yldAnomGau']
    reg['ensoProd']=reg['Harvested_Area']*reg['ENSO']
    reg['naoProd']=reg['Harvested_Area']*reg['NAO']
    reg['iodProd']=reg['Harvested_Area']*reg['IOD']
    reg['tavProd']=reg['Harvested_Area']*reg['TAV']      
    reg.set_index('year',drop=True,inplace=True)
    
    regYld = reg['ProdAnomGau'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    ensoYld =reg['ensoProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    naoYld =reg['naoProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    tavYld =reg['tavProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    iodYld =reg['iodProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
    enYld =(reg['naoProd'].sum(level='year').fillna(0)+reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
    eniYld =(reg['iodProd'].sum(level='year').fillna(0)+reg['naoProd'].sum(level='year').fillna(0)+
              reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
    enitYld =(reg['tavProd'].sum(level='year').fillna(0)+reg['iodProd'].sum(level='year').fillna(0)+reg['naoProd'].sum(level='year').fillna(0)+
              reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')

    if i==0:
        ts = np.zeros([regYld.values.shape[0],len(mRegions)])
        ts[:,i] = regYld.values
        i=i+1
    else:
        ts[:,i] = regYld.values
        i=i+1
    
    eVar2 = 100*(np.linalg.norm(ensoYld.fillna(0).values)**2/np.linalg.norm(regYld.values)**2)
    nVar2 = 100*(np.linalg.norm(naoYld.fillna(0).values)**2/np.linalg.norm(regYld.values)**2)
    iVar2 = 100*(np.linalg.norm(iodYld.fillna(0).values)**2/np.linalg.norm(regYld.values)**2)
    tVar2 = 100*(np.linalg.norm(tavYld.fillna(0).values)**2/np.linalg.norm(regYld.values)**2)
    enitVar2 = 100*(np.linalg.norm(enitYld.values)**2/np.linalg.norm(regYld.values)**2)

    eVar = 1-(np.var(regYld-ensoYld)/np.var(regYld))
    nVar = 1-(np.var(regYld-naoYld)/np.var(regYld))
    iVar = 1-(np.var(regYld-iodYld)/np.var(regYld))
    tVar = 1-(np.var(regYld-tavYld)/np.var(regYld))
    enitVar = 1-(np.var(regYld-enitYld)/np.var(regYld))
    yLim = regYld.abs().max()*1.1

    nB =0
    pB =0
    fig = plt.figure(); ax1 = plt.subplot(121); ax2 = plt.subplot(185)
    if eVar>=0:b1 = ax2.bar(1, eVar2, 0.5, bottom=pB,color='skyblue',label='ENSO');pB=pB+eVar2
    else:b1 = ax2.bar(1, eVar2, 0.5, bottom=nB,color='skyblue',label='ENSO');nB=np.nansum([nB,eVar2]) 
    if nVar>0:b2 = ax2.bar(1, nVar2, 0.5, bottom=pB,color='maroon',label='NAO');pB=pB+nVar2
    else:b2 = ax2.bar(1, nVar2, 0.5, bottom=nB,color='maroon',label='NAO');nB=np.nansum([nB,nVar2]) 
    if iVar>0:b2 = ax2.bar(1, iVar2, 0.5, bottom=pB,color='gold',label='IOD');pB=pB+iVar2
    else:b2 = ax2.bar(1, iVar2, 0.5, bottom=nB,color='gold',label='IOD');nB=np.nansum([nB,iVar2]) 
    if tVar>0:b2 = ax2.bar(1, tVar2, 0.5, bottom=pB,color='peru',label='TAV');pB=pB+tVar2
    else:b2 = ax2.bar(1, tVar2, 0.5, bottom=nB,color='peru',label='TAV');nB=np.nansum([nB,tVar2]) 
    ax1.plot(regYld.sort_index(),'k',lw=1.5);
    ax1.plot(ensoYld.sort_index(),'-',color='cornflowerblue',lw=2)
    ax1.plot(naoYld.sort_index(),'-',color='maroon',lw=2)
    ax1.plot(iodYld.sort_index(),'-',color='gold',lw=2)
    ax1.plot(tavYld.sort_index(),'-',color='peru',lw=2)
    ax2.text(0.67,100-7,str("NAO:   "+str(np.round(nVar2,1))+"%"))
    ax2.text(0.67,100-17,str("IOD:    "+str(np.round(iVar2,1))+"%"))
    ax2.text(0.67,100-27,str("TAV:   "+str(np.round(tVar2,1))+"%"))
    ax2.text(0.67,100-37,str("ENSO: "+str(np.round(eVar2,1))+"%"))
    ax2.text(0.79,pB+5,str(str(np.round(enitVar2,1))+"%"),fontsize=14)
    ax1.set_ylim(-yLim,yLim)
    ax2.set_xlim(.6,1.35)
    ax2.set_xticks([])
    ax2.set_ylim(0,100)
    ax2.set_xlabel('')
    ax1.set_title(regName)
    #ax2.legend()
    fig.set_size_inches(12,4);
    plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/basin comparison/maize Indian Ocean.png')
    #plt.close()    
"""