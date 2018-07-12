#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 08:21:10 2018

@author: weston

Read in groups of states and plot the yield anomaly time series and the running
standard deviation, then the trend in the standard deviation

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

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #

years = [1981,2011]

#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/GADM_dict.py').read())
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/project/Climate Modes - Yield Anoms/ENSOgroups.py').read())

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


for regName in mRegions:
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

    eVar = 1-(np.var(regYld-ensoYld)/np.var(regYld))
    nVar = 1-(np.var(regYld-naoYld)/np.var(regYld))
    iVar = 1-(np.var(regYld-iodYld)/np.var(regYld))
    tVar = 1-(np.var(regYld-tavYld)/np.var(regYld))
    nVarMar = 1-(np.var(regYld-enYld)/np.var(regYld))-eVar
    iVarMar = 1-(np.var(regYld-eniYld)/np.var(regYld))-eVar-nVar
    tVarMar = 1-(np.var(regYld-enitYld)/np.var(regYld))-eVar-nVar-iVar
    yLim = regYld.abs().max()*1.1
    
    nB =0
    pB =0
    fig = plt.figure(); ax1 = plt.subplot(121); ax2 = plt.subplot(164)
    if eVar>=0:b1 = ax2.bar(1, eVar, 0.5, bottom=pB,color='b');pB=pB+eVar
    else:b1 = ax2.bar(1, eVar, 0.5, bottom=nB,color='b');nB=nB+eVar    
    if nVar>0:b2 = ax2.bar(1, nVar, 0.5, bottom=pB,color='r');pB=pB+nVar
    else:b2 = ax2.bar(1, nVar, 0.5, bottom=nB,color='r');nB=nB+nVar
    if iVar>0:b2 = ax2.bar(1, iVar, 0.5, bottom=pB,color='g');pB=pB+iVar
    else:b2 = ax2.bar(1, iVar, 0.5, bottom=nB,color='g');nB=nB+iVar
    if tVar>0:b2 = ax2.bar(1, tVar, 0.5, bottom=pB,color='k');pB=pB+tVar
    else:b2 = ax2.bar(1, tVar, 0.5, bottom=nB,color='k');nB=nB+tVar
    ax1.plot(regYld);
    ax1.plot(ensoYld,'--r')
    ax1.plot(enYld,'--g')
    ax1.plot(eniYld,'--b')
    ax1.plot(enitYld,'--k')
    ax1.set_ylim(-yLim,yLim)
    ax2.set_ylim(0,1)
    fig.set_size_inches(8,4);
    
    plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/time series/maize yield_'+regName+'.png')
    plt.close()
    
    
for regName in wRegions:
    wRegs = ['wheat_'+x.upper () for x in wRegions[regName]]
    reg = crDF[['Harvested_Area','yldAnomGau','ENSO','NAO','IOD','TAV','year']][crDF['state'].isin(wRegs)]
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

    eVar = 1-(np.var(regYld-ensoYld)/np.var(regYld))
    nVar = 1-(np.var(regYld-enYld)/np.var(regYld))-eVar
    iVar = 1-(np.var(regYld-eniYld)/np.var(regYld))-eVar-nVar
    tVar = 1-(np.var(regYld-enitYld)/np.var(regYld))-eVar-nVar-iVar
    yLim = regYld.abs().max()*1.1
    
    nB =0
    pB =0
    fig = plt.figure(); ax1 = plt.subplot(121); ax2 = plt.subplot(164)
    if eVar>=0:b1 = ax2.bar(1, eVar, 0.5, bottom=pB,color='b');pB=pB+eVar
    else:b1 = ax2.bar(1, eVar, 0.5, bottom=nB,color='b');nB=nB+eVar    
    if nVar>0:b2 = ax2.bar(1, nVar, 0.5, bottom=pB,color='r');pB=pB+nVar
    else:b2 = ax2.bar(1, nVar, 0.5, bottom=nB,color='r');nB=nB+nVar
    if iVar>0:b2 = ax2.bar(1, iVar, 0.5, bottom=pB,color='g');pB=pB+iVar
    else:b2 = ax2.bar(1, iVar, 0.5, bottom=nB,color='g');nB=nB+iVar
    if tVar>0:b2 = ax2.bar(1, tVar, 0.5, bottom=pB,color='k');pB=pB+tVar
    else:b2 = ax2.bar(1, tVar, 0.5, bottom=nB,color='k');nB=nB+tVar
    ax1.plot(regYld);
    ax1.plot(ensoYld,'--r')
    ax1.plot(enYld,'--g')
    ax1.plot(eniYld,'--b')
    ax1.plot(enitYld,'--k')
    ax1.set_ylim(-yLim,yLim)
    ax2.set_ylim(0,1)
    fig.set_size_inches(8,4);
    
    plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/time series/wheat yield_'+regName+'.png')
    plt.close()
    
mRegs = [x for x in crDF.state.unique() if x[0:5]=='maize']
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

eVar = 1-(np.var(regYld-ensoYld)/np.var(regYld))
nVar = 1-(np.var(regYld-enYld)/np.var(regYld))-eVar
iVar = 1-(np.var(regYld-eniYld)/np.var(regYld))-eVar-nVar
tVar = 1-(np.var(regYld-enitYld)/np.var(regYld))-eVar-nVar-iVar
yLim = regYld.abs().max()*1.1

nB =0
pB =0
fig = plt.figure(); ax1 = plt.subplot(121); ax2 = plt.subplot(164)
if eVar>=0:b1 = ax2.bar(1, eVar, 0.5, bottom=pB,color='b');pB=pB+eVar
else:b1 = ax2.bar(1, eVar, 0.5, bottom=nB,color='b');nB=nB+eVar    
if nVar>0:b2 = ax2.bar(1, nVar, 0.5, bottom=pB,color='r');pB=pB+nVar
else:b2 = ax2.bar(1, nVar, 0.5, bottom=nB,color='r');nB=nB+nVar
if iVar>0:b2 = ax2.bar(1, iVar, 0.5, bottom=pB,color='g');pB=pB+iVar
else:b2 = ax2.bar(1, iVar, 0.5, bottom=nB,color='g');nB=nB+iVar
if tVar>0:b2 = ax2.bar(1, tVar, 0.5, bottom=pB,color='k');pB=pB+tVar
else:b2 = ax2.bar(1, tVar, 0.5, bottom=nB,color='k');nB=nB+tVar
ax1.plot(regYld);
ax1.plot(ensoYld,'--r')
ax1.plot(enYld,'--g')
ax1.plot(eniYld,'--b')
ax1.plot(enitYld,'--k')
ax1.set_ylim(-yLim,yLim)
ax2.set_ylim(0,1)
fig.set_size_inches(8,4);

plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/time series/maize yield_global.png')
plt.close()


mRegs = [x for x in crDF.state.unique() if x[0:5]=='wheat']
reg = crDF[['Harvested_Area','yldAnomGau','ENSO','NAO','IOD','TAV','year']][crDF['state'].isin(mRegs)]
reg['ProdAnomGau']=reg['Harvested_Area']*reg['yldAnomGau']
reg['ensoProd']=reg['Harvested_Area']*reg['ENSO']
reg['naoProd']=reg['Harvested_Area']*reg['NAO']
reg['iodProd']=reg['Harvested_Area']*reg['IOD']
reg['tavProd']=reg['Harvested_Area']*reg['TAV']
reg.set_index('year',drop=True,inplace=True)

regYld = reg['ProdAnomGau'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
ensoYld =reg['ensoProd'].sum(level='year')/reg['Harvested_Area'].sum(level='year')
enYld =(reg['naoProd'].sum(level='year').fillna(0)+reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
eniYld =(reg['iodProd'].sum(level='year').fillna(0)+reg['naoProd'].sum(level='year').fillna(0)+
          reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')
enitYld =(reg['tavProd'].sum(level='year').fillna(0)+reg['iodProd'].sum(level='year').fillna(0)+reg['naoProd'].sum(level='year').fillna(0)+
          reg['ensoProd'].sum(level='year').fillna(0))/reg['Harvested_Area'].sum(level='year')

eVar = 1-(np.var(regYld-ensoYld)/np.var(regYld))
nVar = 1-(np.var(regYld-enYld)/np.var(regYld))-eVar
iVar = 1-(np.var(regYld-eniYld)/np.var(regYld))-eVar-nVar
tVar = 1-(np.var(regYld-enitYld)/np.var(regYld))-eVar-nVar-iVar
yLim = regYld.abs().max()*1.1

nB =0
pB =0
fig = plt.figure(); ax1 = plt.subplot(121); ax2 = plt.subplot(164)
if eVar>=0:b1 = ax2.bar(1, eVar, 0.5, bottom=pB,color='b');pB=pB+eVar
else:b1 = ax2.bar(1, eVar, 0.5, bottom=nB,color='b');nB=nB+eVar    
if nVar>0:b2 = ax2.bar(1, nVar, 0.5, bottom=pB,color='r');pB=pB+nVar
else:b2 = ax2.bar(1, nVar, 0.5, bottom=nB,color='r');nB=nB+nVar
if iVar>0:b2 = ax2.bar(1, iVar, 0.5, bottom=pB,color='g');pB=pB+iVar
else:b2 = ax2.bar(1, iVar, 0.5, bottom=nB,color='g');nB=nB+iVar
if tVar>0:b2 = ax2.bar(1, tVar, 0.5, bottom=pB,color='k');pB=pB+tVar
else:b2 = ax2.bar(1, tVar, 0.5, bottom=nB,color='k');nB=nB+tVar
ax1.plot(regYld);
ax1.plot(ensoYld,'--r')
ax1.plot(enYld,'--g')
ax1.plot(eniYld,'--b')
ax1.plot(enitYld,'--k')
ax1.set_ylim(-yLim,yLim)
ax2.set_ylim(0,1)
fig.set_size_inches(8,4);

plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/time series/wheat yield_global.png')
plt.close()