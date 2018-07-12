#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 11:45:44 2018
@author: weston
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from unidecode import unidecode

eofs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+
        '/all crops/crop_all crops_Jan-DecENSOyr0har_EOFs_1981-2011_janSVD_allTrops_25N25S.npy')
pcs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+
        '/all crops/crop_all crops_Jan-DecENSOyr0har_PCs_1981-2011_janSVD_allTrops_25N25S.npy')
CRdata = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/'+
                    'stats regressions/all cropsVar1961_2012_Jan-DecENSOyr0har_janSVD_allTrops_25N25S.csv')
CRdata['temp'] =1
ha = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/SVD/'+
                           'all crops1961_2012_Jan-DecENSOyr0har_SVD2_janSVD_allTrops_25N25S.csv')
ha = ha[['Harvested_Area','expectedYldGau','country.name']][ha.year==2010]

CRdata = pd.concat([CRdata,ha],axis=1)
CRdata.dropna(subset=['temp'],inplace=True)

CRdata.reset_index(inplace=True)
eof_states=list(CRdata['index']) #find the states  

eof_prod = -0.5*eofs[0,:]*np.std(pcs[:,0])+0.5*eofs[1,:]*np.std(pcs[:,1])

#Read in the crop yield anomaly data
year=1983
crDFyr = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/country data/Dataframes/maize.csv') 
crDFyr = crDFyr[['state','yldAnomGau']].loc[crDFyr.year==year]
crDFyr.state=[unidecode('maize_'+x.upper()) for x in crDFyr.state] #make state names the same
crDFyr=crDFyr.rename(columns={'state':'index'})
CRdata=CRdata.merge(crDFyr,'left')

SVD = pd.DataFrame(data={'state':eof_states,'SVD':eof_prod,'HA':CRdata['Harvested_Area'],'expYld':CRdata['expectedYldGau'],
                        'prod':(eof_prod*CRdata['expectedYldGau'])*CRdata['Harvested_Area'],
                        'country':CRdata['country.name'],'maizeYRprd':CRdata['yldAnomGau']*CRdata['expectedYldGau']*CRdata['Harvested_Area']})
SVD.set_index('state',drop=False,inplace=True)
SVD['country'].loc[SVD.state.isin(['maize_PR','maize_SC','maize_RS',
       'wheat_PR','wheat_SC','wheat_RS','soy_PR','soy_SC','soy_RS'])] = 'S. Brazil'

SVD['maizeYRprd'].loc[np.isnan(SVD.maizeYRprd)]=0

regions = {'southeast South America':['Argentina','Uruguay','Paraguay','S. Brazil'],
           'US':['United States'],'India':['India'],'Mexico':['Mexico'],
           'Brazil':['Brazil'],'W. Africa':["Côte d'Ivoire",'Senegal','Gambia',
            'Guinea-Bissau','Mali','Sierra Leone','Ghana','Togo','Benin','Burkina Faso',
            'Nigeria','Cameroon'],'China':['China'],'Southeast Africa': ['South Africa',
            'Zimbabwe','Mozambique','Zambia','Botswana'],'Australia':['Australia'],
           'Eurasia-IND+CHN':['Spain','Portugal','France','Germany','United Kingdom',
            'Afghanistan','Iran','Sweden','Denmark','Norway','Turkey','Romania','Syria',
            'Pakistan','Iraq','Ireland','Poland','Austria','Belarus','Greece','Italy',
            'Mongolia','Yugoslavia','Czechoslovakia','Belgium-Luxembourg','USSR'],
             'East Africa':['Ethiopia','Kenya','Somalia','Somaliland','Eritrea',
              'Djibouti','Uganda'],'Canada':['Canada']}
print('\n\n')
for crop in ['maize','soy','wheat']:
        num=20
        cStates = [x for x in SVD.state.unique() if x[0:len(crop)]==crop]
        cData = SVD[SVD.index.isin(cStates)].sort_values(by='prod')
        prodSum = round(cData['prod'].values.sum()/1000000000,2)
        rData = cData.copy()
        cData = cData[np.abs(cData['prod'])>np.sort(np.abs(cData['prod']))[-num]]
        plt.figure()
        plt.title(crop)
        plt.bar(range(cData['prod'].values.size),cData['prod'].values)
        plt.xticks(range(cData['prod'].values.size),cData.state,rotation=90)
        plt.tight_layout()
        plt.text(cData['prod'].values.size/2.5,500000000,str(prodSum)+' M tonns')
        plt.show();plt.close()
        print(crop)
        print('global: '+str(prodSum))
        for reg in regions:
            regData = rData[rData.country.isin(regions[reg])]
            regSum = round(regData['prod'].values.sum()/1000000000,2)
            print(reg+' '+str(regSum))
        print('\n\n')
        
crop='maize'
cStates = [x for x in SVD.state.unique() if x[0:len(crop)]==crop]
cData = SVD[SVD.index.isin(cStates)].sort_values(by='maizeYRprd')
prodSum = round(cData['maizeYRprd'].values.sum()/1000000000,2)
rData = cData.copy()
cData = cData[np.abs(cData['maizeYRprd'])>np.sort(np.abs(cData['maizeYRprd']))[-num]]
print("Maize in 1983 [Millions of tonnes]")
print('global: '+str(prodSum))
for reg in regions:
    regData = rData[rData.country.isin(regions[reg])]
    regSum = round(regData['maizeYRprd'].values.sum()/1000000000,2)
    print(reg+' '+str(regSum))
print('\n\n')