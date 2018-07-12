#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 17:55:02 2018

@author: weston
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
import pandas as pd
from scipy import stats
from eofs.standard import Eof
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #

yrMin = 1981
yrMax = 2011
PCs = 3 #python indexing
monCutoff = 1
timePeriod =''# ' 1960-2010' #'' or ' 1960-2010'
rot = '' #'varimax' or ''
notes = '_marSVD'

#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = months[monCutoff-1]+'-'+months[monCutoff-2]+'ENSOyr0har'

#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                 Helper Functions                  #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#Read in moving average: moving_average(a, smooth)
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/movingAverage.py').read())
#Read in moving standard deviation: moving_std(a, smooth)
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/movingStd.py').read())
#create surrogate timeseries
#follows from Ebusuzaki (1997) and Schrieber and Shmitz (2000)
#  function call is   ############ surrogate(input TS): ############                  #
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/surrogateData.py').read())
#Read in the shapefiles
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/GADM_dict.py').read())
#import the varimax rotation
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/varimax.py').read())
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/project/Climate Modes - Yield Anoms/NAO state list.py').read())

# Restrict SSTs to tropical Pacific and Indian Ocean
latMax = 5; latMin = -5; lonMax = 270; lonMin = 160
#latMax = 10; latMin = -10; lonMax = 280; lonMin = 160
#latMax = 5; latMin = -5; lonMax = 240; lonMin = 190
#latMax = 15; latMin = -15; lonMax = 360; lonMin = 0

#Read in the NAO time series to plot against
NAOdf = pd.read_csv('/Volumes/Data_Archive/Data/NAO/nao_station_seasonal.csv',index_col=0)
nao = NAOdf.loc[yrMin:yrMax]['DJF'].values


#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   Read in the data                        #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#         
#crop data first
wDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/'+
                                'stats regressions'+timePeriod+'/wheat1961_2012_'+ENSOyr+'.csv')   
wDF.state=[x.upper() for x in wDF.state] #make state names the same
wDF['state'] = ',wheat_'.join(wDF['state'].values).split(',')
wDF['state'].iloc[0] = str('wheat_'+wDF['state'].iloc[0])
#count up total production before dropping any states for having a short production time series and before adjusting the years for the SVD
wDF=wDF[np.isfinite(wDF.GAM_E)] #drop incomplete entries
wDF.state[wDF.state=='wheat_NEW SOUTH WALES(B)']='wheat_NEW SOUTH WALES'

mDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/'+
                                'stats regressions'+timePeriod+'/maize1961_2012_'+ENSOyr+'.csv')   
mDF.state=[x.upper() for x in mDF.state] #make state names the same
mDF['state'] = ',maize_'.join(mDF['state'].values).split(',')
mDF['state'].iloc[0] = str('maize_'+mDF['state'].iloc[0])
#count up total production before dropping any states for having a short production time series and before adjusting the years for the SVD
mDF=mDF[np.isfinite(mDF.GAM_E)] #drop incomplete entries

sDF = pd.DataFrame.from_csv('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/'+
                                'stats regressions'+timePeriod+'/soy1961_2012_'+ENSOyr+'.csv')   
sDF.state=[x.upper() for x in sDF.state] #make state names the same
sDF['state'] = ',soy_'.join(sDF['state'].values).split(',')
sDF['state'].iloc[0] = str('soy_'+sDF['state'].iloc[0])
#count up total production before dropping any states for having a short production time series and before adjusting the years for the SVD
sDF=sDF[np.isfinite(sDF.GAM_E)] #drop incomplete entries

cropDF =wDF.append(mDF.append(sDF,ignore_index=True),ignore_index=True)
cropDF.set_index(['state','year'],drop=True,inplace=True) #re-index the data
obs = cropDF['yldAnomGau']*cropDF['expectedYldGau']*cropDF['Harvested_Area']
exp = cropDF['expectedYldGau']*cropDF['Harvested_Area']
cropDF.reset_index(inplace=True) #re-index the data to alter year
cropDF.year = cropDF.year - cropDF.yrAdded #adjust years for SVD matrix
CRpctAnom = cropDF.pivot(index='year',columns='state',values='yldAnomGau') #table for SVD
CRpctAnom=CRpctAnom.loc[yrMin:yrMax,:] #crop to the correct years
CRpctAnom=CRpctAnom.dropna(axis=1) #drop incomplete columns
eof_states=list(CRpctAnom.columns) #find the states 

#Duplicate the DF for the NAO analysis with only the NAO states
NAOstates=[x.upper() for x in NAOstates] #make state names the same
wNAOstates = ',wheat_'.join(NAOstates).split(',')
wNAOstates[0] = str('wheat_'+wNAOstates[0])
mNAOstates = ',maize_'.join(NAOstates).split(',')
mNAOstates[0] = str('maize_'+mNAOstates[0])
sNAOstates = ',soy_'.join(NAOstates).split(',')
sNAOstates[0] = str('soy_'+sNAOstates[0])
crNAOstates = wNAOstates+mNAOstates+sNAOstates

#claculate the EOFs for the global analysis
solver = Eof(np.array(CRpctAnom))
eofs = solver.eofsAsCorrelation(neofs=PCs)
if rot == 'varimax': #rotate with varimax
    eofs = varimax(eofs.T);eofs = eofs.T
pcs = solver.pcs(npcs=PCs, pcscaling=1)
eigens = solver.varianceFraction(10)
errors = solver.northTest(10,vfscaled=True) #calculate the errors on EOFs


fig = plt.figure(); ax1 = plt.subplot(111)
ax1.set_title('Percent Variance Explained');
ax1.plot(range(1,11),eigens*100,'ko')
ax1.vlines(range(1,11),(eigens-errors)*100,(eigens+errors)*100)
fig.set_size_inches(4,8);
fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/EOF/all crops_EigenVectors_'+ENSOyr+notes+'.png')
plt.close()

cEOFs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/all crops/crop_all crops'+
              '_'+ENSOyr+'_EOFs_1981-2011'+notes+'.npy') 
cPCs = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD/all crops/crop_all crops'+
              '_'+ENSOyr+'_PCs_1981-2011'+notes+'.npy')

#compare to Niño 3.4 and to the second PC of mode 2
Nino34 = pd.read_csv('/Volumes/Data_Archive/Data/ENSO/NinoIndex34.csv',',',index_col=0,skiprows=1)
enso = Nino34.loc[yrMin:yrMax]['1'].values
plt.figure()
plt.plot(-pcs[:,0]/np.std(pcs[:,0]),'b');
plt.plot(enso/np.std(enso),'k');
plt.plot(cPCs[:,1]/np.std(cPCs[:,1]),'--b');#plt.show()

print(np.corrcoef(pcs[:,0],enso)[0][1])
print(np.corrcoef(cPCs[:,1],enso)[0][1])
print(np.corrcoef(cPCs[:,1],pcs[:,0])[0][1])
#EOF analysis has a combined mode for ENSO. This makes sense because the
# SVD shows expansion coefficients that are nearly identical (or at least
# highly correlated)
print(np.corrcoef(cPCs[:-1,1],cPCs[1:,0])[0][1])

rhoS = []
pS=[]
rhoP = []
for j in range(3):
    for i in range(PCs):
        print('SVD'+str(j+1)+' with EOF'+str(i+1)+': '+str(np.round(stats.spearmanr(cEOFs[j,:],eofs[i,:]),2)))
        rhoS.append(stats.spearmanr(cEOFs[j,:],eofs[i,:])[0])
        pS.append(stats.spearmanr(cEOFs[j,:],eofs[i,:])[1])
        rhoP.append(np.corrcoef(cEOFs[j,:],eofs[i,:])[0][1])


for PC in range(PCs+1):
    continue
    eof_prod = eofs[PC,:]   
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
    #                   Plot the Regions                #
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#         
    for crop in ['wheat','maize','soy']:
        eof_states = eof_states.copy()
        #Correct the names to match the mapping names     
    
        thisCMAP =  'BrBG'# matplotlib.colors.Colormap(Margot2_4.mpl_colormap)
        clrBr = np.arange(-.45,0.45+.05,.05)
        norm = Normalize(vmin=-0.45, vmax=0.45, clip=False)
        mapper=cm.ScalarMappable(norm=norm, cmap=thisCMAP)
        
        fig = plt.figure()
        m = Basemap(projection='kav7',lon_0=0)
        #m.pcolor(sstLons.astype(int),sstLats.astype(int),eof_sst,cmap='RdBu_r',latlon=True,vmin=-0.05,vmax=0.05)
        
        
        SVDcol = pd.DataFrame(data={'state':eof_states,'SVD':eof_prod})
        
        
        #Plot everything at the country level first from FAO data before going into subnational data
        world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                name='cnts', drawbounds=True)
        names = [];colors = {};i=0
        for shape_dict in m.cnts_info:
            names.append(shape_dict['NAME'])
            if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['NAME']]=state_color
            #else: colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.cnts):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
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
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            #else: colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly)
        mapper.set_array(clrBr)
        
        BRshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/BRA_adm_shp/BRA_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if (shape_dict['HASC_1']=='BR.AM'):
                continue#colors[shape_dict['HASC_1']]=0
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['HASC_1'][3:5]))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['HASC_1'][3:5])]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly) 
        mapper.set_array(clrBr)
        
        ARshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/ARG_adm_shp/ARG_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if (shape_dict['HASC_1']=='AR.DF')|(shape_dict['HASC_1']=='AR.SC')|(shape_dict['HASC_1']=='AR.TF'):
                continue#colors[shape_dict['HASC_1']]=0
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
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
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly) 
        mapper.set_array(clrBr)
        
        if crop is 'wheat':
            AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly) 
            mapper.set_array(clrBr) 
    
            CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                            name='states', drawbounds=True)
            # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly) 
        
            INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            mapper.set_array(clrBr)   
        
        if crop is 'maize':
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                            name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            mapper.set_array(clrBr)
        
            INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            mapper.set_array(clrBr) 
        
        
        cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
        cb = ColorbarBase(cax,cmap=thisCMAP,norm=norm,orientation='horizontal')
        mapper.set_array(clrBr);
        fig.set_size_inches(30,30)
        plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/EOF/'+crop+'_EOF'+str(PC+1)+ENSOyr+str(PCs)+'PCs'+rot+notes+'.png')
        plt.close()
        
        
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   NAO analysis.                   #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#   
CRpctAnomNao = CRpctAnom[list(set(CRpctAnom.columns).intersection(crNAOstates))]
WpctAnomNao = CRpctAnom[list(set(CRpctAnom.columns).intersection(wNAOstates))]


#Repeat for the NAO analysis
solverNao = Eof(np.array(CRpctAnomNao))
eofsNao = solverNao.eofsAsCorrelation(neofs=PCs)
if rot == 'varimax': #rotate with varimax
    eofsNao = varimax(eofsNao.T);eofsNao = eofsNao.T
pcsNao = solverNao.pcs(npcs=PCs, pcscaling=1)
eigensNao = solverNao.varianceFraction(10)
errorsNao = solverNao.northTest(10,vfscaled=True) #calculate the errors on EOFs

#Repeat for the NAO analysis with only wheat
solverNaoW = Eof(np.array(WpctAnomNao))
eofsNaoW = solverNaoW.eofsAsCorrelation(neofs=PCs)
if rot == 'varimax': #rotate with varimax
    eofsNaoW = varimax(eofsNaoW.T);eofsNaoW = eofsNaoW.T
pcsNaoW = solverNaoW.pcs(npcs=PCs, pcscaling=1)
eigensNaoW = solverNaoW.varianceFraction(10)
errorsNaoW = solverNaoW.northTest(10,vfscaled=True) #calculate the errors on EOFs

fig = plt.figure(); ax1 = plt.subplot(111)
ax1.set_title('Percent Variance Explained');
ax1.plot(range(1,11),eigensNao*100,'ko')
ax1.vlines(range(1,11),(eigensNao-errorsNao)*100,(eigensNao+errorsNao)*100)
fig.set_size_inches(4,8);#plt.show()
fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/EOF/NAO_all crops_EigenVectors_'+ENSOyr+notes+'.png')
plt.close()

fig = plt.figure(); ax1 = plt.subplot(111)
ax1.set_title('Percent Variance Explained');
ax1.plot(range(1,11),eigensNaoW*100,'ko')
ax1.vlines(range(1,11),(eigensNaoW-errorsNaoW)*100,(eigensNaoW+errorsNaoW)*100)
fig.set_size_inches(4,8);#plt.show()
fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/EOF/NAO_wheat_EigenVectors_'+ENSOyr+notes+'.png')
plt.close()


for iPC in range(3):
    print('\n\n')
    print(iPC+1)
    for seas in NAOdf.columns:
        print(seas)
        print(np.corrcoef(pcsNao[:,iPC],NAOdf.loc[yrMin:yrMax][seas])[0][1])
        print(np.corrcoef(pcsNao[:,iPC],NAOdf.loc[yrMin-1:yrMax-1][seas])[0][1])
        print(np.corrcoef(pcsNao[:,iPC],NAOdf.loc[yrMin+1:yrMax+1][seas])[0][1])
    print('annual avg NAO')
    print(np.corrcoef(pcsNao[:,iPC],np.mean(NAOdf,1).loc[yrMin-1:yrMax-1].values)[0][1])
    print(np.corrcoef(pcsNao[:,iPC],np.mean(NAOdf,1).loc[yrMin:yrMax].values)[0][1])
    print(np.corrcoef(pcsNao[:,iPC],np.mean(NAOdf,1).loc[yrMin+1:yrMax+1].values)[0][1])

plt.figure()
plt.plot(pcsNao[:,0]/np.std(pcsNao[:,0]),'k');
plt.plot(-1*NAOdf.loc[yrMin-1:yrMax-1]['DJF'].values/np.std(NAOdf.loc[yrMin:yrMax]['DJF'].values),'b');#plt.show()
#plt.plot(NAOdf.loc[yrMin:yrMax]['JFM'].values/np.std(NAOdf.loc[yrMin:yrMax]['JFM'].values),'b');#plt.show()
plt.plot(NAOdf.loc[yrMin+1:yrMax+1]['SON'].values/np.std(NAOdf.loc[yrMin:yrMax]['SON'].values),'r');#plt.show()

plt.figure()
plt.plot(pcsNao[:,0]/np.std(pcsNao[:,0]),'k');
plt.plot(np.mean(NAOdf,1).loc[yrMin:yrMax].values/np.std(np.mean(NAOdf,1).loc[yrMin:yrMax].values),'b');#plt.show()


for iPC in range(3):
    print('\n\n')
    print(iPC+1)
    for seas in NAOdf.columns:
        print(seas)
        print(np.corrcoef(pcsNaoW[:,iPC],NAOdf.loc[yrMin:yrMax][seas])[0][1])
        print(np.corrcoef(pcsNaoW[:,iPC],NAOdf.loc[yrMin-1:yrMax-1][seas])[0][1])
        print(np.corrcoef(pcsNaoW[:,iPC],NAOdf.loc[yrMin+1:yrMax+1][seas])[0][1])
    print('annual avg NAO')
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf,1).loc[yrMin-1:yrMax-1].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf,1).loc[yrMin:yrMax].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf,1).loc[yrMin+1:yrMax+1].values)[0][1])


#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
#                   NAO analysis                    #
#               after removing ENSO                 #
#^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#   
ensoYlds = pd.read_pickle('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/'+
                          'SVD/all crops/crop_all crops_Jan-DecENSOyr0har_recon_2EOFs_1981-2011_marSVD.pkl')

CRpctAnomLessENSO = CRpctAnom - ensoYlds
CRpctAnomNao = CRpctAnomLessENSO[list(set(CRpctAnomLessENSO.columns).intersection(crNAOstates))]
WpctAnomNao = CRpctAnomLessENSO[list(set(CRpctAnomLessENSO.columns).intersection(wNAOstates))]


#Repeat for the NAO analysis
solverNao = Eof(np.array(CRpctAnomNao))
eofsNao = solverNao.eofsAsCorrelation(neofs=PCs)
if rot == 'varimax': #rotate with varimax
    eofsNao = varimax(eofsNao.T);eofsNao = eofsNao.T
pcsNao = solverNao.pcs(npcs=PCs, pcscaling=1)
eigensNao = solverNao.varianceFraction(10)
errorsNao = solverNao.northTest(10,vfscaled=True) #calculate the errors on EOFs

#Repeat for the NAO analysis with only wheat
solverNaoW = Eof(np.array(WpctAnomNao))
eofsNaoW = solverNaoW.eofsAsCorrelation(neofs=PCs)
if rot == 'varimax': #rotate with varimax
    eofsNaoW = varimax(eofsNaoW.T);eofsNaoW = eofsNaoW.T
pcsNaoW = solverNaoW.pcs(npcs=PCs, pcscaling=1)
eigensNaoW = solverNaoW.varianceFraction(10)
errorsNaoW = solverNaoW.northTest(10,vfscaled=True) #calculate the errors on EOFs

fig = plt.figure(); ax1 = plt.subplot(111)
ax1.set_title('Percent Variance Explained');
ax1.plot(range(1,11),eigensNao*100,'ko')
ax1.vlines(range(1,11),(eigensNao-errorsNao)*100,(eigensNao+errorsNao)*100)
fig.set_size_inches(4,8);#plt.show()

fig = plt.figure(); ax1 = plt.subplot(111)
ax1.set_title('Percent Variance Explained');
ax1.plot(range(1,11),eigensNaoW*100,'ko')
ax1.vlines(range(1,11),(eigensNaoW-errorsNaoW)*100,(eigensNaoW+errorsNaoW)*100)
fig.set_size_inches(4,8);#plt.show()


for iPC in range(PCs):
    print('\n\n')
    print(iPC+1)
    for seas in NAOdf.columns:
        print(seas)
        print(np.corrcoef(pcsNao[:,iPC],NAOdf.loc[yrMin-1:yrMax-1][seas])[0][1])
        print(np.corrcoef(pcsNao[:,iPC],NAOdf.loc[yrMin:yrMax][seas])[0][1])
        print(np.corrcoef(pcsNao[:,iPC],NAOdf.loc[yrMin+1:yrMax+1][seas])[0][1])
    print('annual avg NAO')
    print(np.corrcoef(pcsNao[:,iPC],np.mean(NAOdf,1).loc[yrMin-1:yrMax-1].values)[0][1])
    print(np.corrcoef(pcsNao[:,iPC],np.mean(NAOdf,1).loc[yrMin:yrMax].values)[0][1])
    print(np.corrcoef(pcsNao[:,iPC],np.mean(NAOdf,1).loc[yrMin+1:yrMax+1].values)[0][1])

plt.figure()
plt.plot(pcsNao[:,0]/np.std(pcsNao[:,0]),'k');
plt.plot(NAOdf.loc[yrMin-1:yrMax-1]['DJF'].values/np.std(NAOdf.loc[yrMin:yrMax]['DJF'].values),'b');#plt.show()
plt.plot(NAOdf.loc[yrMin:yrMax]['JAS'].values/np.std(NAOdf.loc[yrMin:yrMax]['JAS'].values),'r');#plt.show()


plt.figure()
plt.plot(pcsNao[:,1]/np.std(pcsNao[:,1]),'k');
plt.plot(NAOdf.loc[yrMin:yrMax]['DJF'].values/np.std(NAOdf.loc[yrMin:yrMax]['DJF'].values),'b');#plt.show()
plt.plot(NAOdf.loc[yrMin:yrMax]['ASO'].values/np.std(NAOdf.loc[yrMin:yrMax]['ASO'].values),'r');#plt.show()



for iPC in range(PCs):
    print('\n\n')
    print(iPC+1)
    for seas in NAOdf.columns:
        print(seas)
        print(np.corrcoef(pcsNaoW[:,iPC],NAOdf.loc[yrMin-1:yrMax-1][seas])[0][1])
        print(np.corrcoef(pcsNaoW[:,iPC],NAOdf.loc[yrMin:yrMax][seas])[0][1])
        print(np.corrcoef(pcsNaoW[:,iPC],NAOdf.loc[yrMin+1:yrMax+1][seas])[0][1])
    print('annual avg NAO')
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf,1).loc[yrMin-1:yrMax-1].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf,1).loc[yrMin:yrMax].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf,1).loc[yrMin+1:yrMax+1].values)[0][1])

    print('Dec - Sep avg NAO')
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf.iloc[:,:9],1).loc[yrMin-1:yrMax-1].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf.iloc[:,:9],1).loc[yrMin:yrMax].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf.iloc[:,:9],1).loc[yrMin+1:yrMax+1].values)[0][1])

    print('Dec - Mar avg NAO')
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf.iloc[:,:4],1).loc[yrMin-1:yrMax-1].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf.iloc[:,:4],1).loc[yrMin:yrMax].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf.iloc[:,:4],1).loc[yrMin+1:yrMax+1].values)[0][1])

    print('DJF -1 + ASO 0 avg NAO')
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf[['DJF','ASO']],1).loc[yrMin-1:yrMax-1].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf[['DJF','ASO']],1).loc[yrMin:yrMax].values)[0][1])
    print(np.corrcoef(pcsNaoW[:,iPC],np.mean(NAOdf[['DJF','ASO']],1).loc[yrMin+1:yrMax+1].values)[0][1])



plt.figure()
plt.plot(pcsNaoW[:,0]/np.std(pcsNaoW[:,0]),'k');
plt.plot(-1*NAOdf.loc[yrMin-1:yrMax-1]['DJF'].values/np.std(NAOdf.loc[yrMin:yrMax]['DJF'].values),'b');#plt.show()
#plt.plot(NAOdf.loc[yrMin:yrMax]['JFM'].values/np.std(NAOdf.loc[yrMin:yrMax]['JFM'].values),'b');#plt.show()
plt.plot(-1*NAOdf.loc[yrMin:yrMax]['ASO'].values/np.std(NAOdf.loc[yrMin:yrMax]['ASO'].values),'r');#plt.show()


plt.figure()
plt.plot(pcsNaoW[:,1]/np.std(pcsNaoW[:,1]),'k');
plt.plot(NAOdf.loc[yrMin:yrMax]['AMJ'].values/np.std(NAOdf.loc[yrMin:yrMax]['AMJ'].values),'darkgreen');#plt.show()



for PC in range(PCs):
    eof_prod = eofsNao[PC,:]   
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
    #                   Plot the Regions                #
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#         
    for crop in ['wheat','maize','soy']:
        eof_states = list(set(CRpctAnom.columns).intersection(crNAOstates))

    
        thisCMAP =  'BrBG'# matplotlib.colors.Colormap(Margot2_4.mpl_colormap)
        clrBr = np.arange(-.45,0.45+.05,.05)
        norm = Normalize(vmin=-0.45, vmax=0.45, clip=False)
        mapper=cm.ScalarMappable(norm=norm, cmap=thisCMAP)
        
        fig = plt.figure()
        m = Basemap(projection='kav7',lon_0=0)
        #m.pcolor(sstLons.astype(int),sstLats.astype(int),eof_sst,cmap='RdBu_r',latlon=True,vmin=-0.05,vmax=0.05)
        
        
        SVDcol = pd.DataFrame(data={'state':eof_states,'SVD':eof_prod})
        
        
        #Plot everything at the country level first from FAO data before going into subnational data
        world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                name='cnts', drawbounds=True)
        names = [];colors = {};i=0
        for shape_dict in m.cnts_info:
            names.append(shape_dict['NAME'])
            if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['NAME']]=state_color
            #else: colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.cnts):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
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
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            #else: colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly)
        mapper.set_array(clrBr)
        
        BRshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/BRA_adm_shp/BRA_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if (shape_dict['HASC_1']=='BR.AM'):
                continue#colors[shape_dict['HASC_1']]=0
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['HASC_1'][3:5]))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['HASC_1'][3:5])]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly) 
        mapper.set_array(clrBr)
        
        ARshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/ARG_adm_shp/ARG_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if (shape_dict['HASC_1']=='AR.DF')|(shape_dict['HASC_1']=='AR.SC')|(shape_dict['HASC_1']=='AR.TF'):
                continue#colors[shape_dict['HASC_1']]=0
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
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
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly) 
        mapper.set_array(clrBr)
        
        if crop is 'wheat':
            AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly) 
            mapper.set_array(clrBr) 
    
            CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                            name='states', drawbounds=True)
            # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly) 
        
            INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            mapper.set_array(clrBr)   
        
        if crop is 'maize':
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                            name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            mapper.set_array(clrBr)
        
            INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            mapper.set_array(clrBr) 
        
        
        cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
        cb = ColorbarBase(cax,cmap=thisCMAP,norm=norm,orientation='horizontal')
        mapper.set_array(clrBr);
        fig.set_size_inches(30,30)
        plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/EOF/NAO_'+crop+'_EOF'+str(PC+1)+ENSOyr+str(PCs)+'PCs'+rot+notes+'.png')
        plt.close()
        



for PC in range(PCs):
    eof_prod = eofsNaoW[PC,:]   
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#
    #                   Plot the Regions                #
    #^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*~^~*#         
    for crop in ['wheat']:
        eof_states = list(set(CRpctAnomLessENSO.columns).intersection(wNAOstates))
    
        thisCMAP =  'BrBG'# matplotlib.colors.Colormap(Margot2_4.mpl_colormap)
        clrBr = np.arange(-.45,0.45+.05,.05)
        norm = Normalize(vmin=-0.45, vmax=0.45, clip=False)
        mapper=cm.ScalarMappable(norm=norm, cmap=thisCMAP)
        
        fig = plt.figure()
        m = Basemap(projection='kav7',lon_0=0)
        #m.pcolor(sstLons.astype(int),sstLats.astype(int),eof_sst,cmap='RdBu_r',latlon=True,vmin=-0.05,vmax=0.05)
        
        
        SVDcol = pd.DataFrame(data={'state':eof_states,'SVD':eof_prod})
        
        
        #Plot everything at the country level first from FAO data before going into subnational data
        world = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/ne_50m_admin_0_countries/ne_50m_admin_0_countries',
                                name='cnts', drawbounds=True)
        names = [];colors = {};i=0
        for shape_dict in m.cnts_info:
            names.append(shape_dict['NAME'])
            if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['NAME']]=state_color
            #else: colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.cnts):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
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
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            #else: colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly)
        mapper.set_array(clrBr)
        
        BRshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/BRA_adm_shp/BRA_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if (shape_dict['HASC_1']=='BR.AM'):
                continue#colors[shape_dict['HASC_1']]=0
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['HASC_1'][3:5]))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['HASC_1'][3:5])]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly) 
        mapper.set_array(clrBr)
        
        ARshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/ARG_adm_shp/ARG_adm1', 
                        name='states', drawbounds=True)
        # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
        names = [];colors = {};i=0
        for shape_dict in m.states_info:
            names.append(shape_dict['HASC_1'])
            if (shape_dict['HASC_1']=='AR.DF')|(shape_dict['HASC_1']=='AR.SC')|(shape_dict['HASC_1']=='AR.TF'):
                continue#colors[shape_dict['HASC_1']]=0
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
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
            elif len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                colors[shape_dict['HASC_1']]=state_color
            else: continue#colors[shape_dict['HASC_1']]=0
        ax = plt.gca() # get current axes instance
        for nshape, seg in enumerate(m.states):
            if names[nshape] not in colors: continue
            poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
            ax.add_patch(poly) 
        mapper.set_array(clrBr)
        
        if crop is 'wheat':
            AUshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/AUS_adm_shp/AUS_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly) 
            mapper.set_array(clrBr) 
    
            CAshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/CAN_adm_shp/CAN_adm1', 
                            name='states', drawbounds=True)
            # collect the state names from the shapefile attributes so we can look up the shape obect for a state by it's name
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly) 
        
            INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            mapper.set_array(clrBr)   
        
        if crop is 'maize':
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            MXshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/MEX_adm_shp/MEX_adm1', 
                            name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            mapper.set_array(clrBr)
        
            INDshp = m.readshapefile('/Volumes/Data_Archive/Data/adminBoundaries/GADM/IND_adm_shp/IND_adm1', 
                        name='states', drawbounds=True)
            names = [];colors = {};i=0
            for shape_dict in m.states_info:
                names.append(shape_dict['HASC_1'])
                if len(np.where(np.array(eof_states)==str(crop+'_'+shape_dict['NAME_1'].upper()))[0]) != 0:
                    state_color = SVDcol[SVDcol.state==str(crop+'_'+shape_dict['NAME_1'].upper())]['SVD'].values.astype(float)[0]
                    colors[shape_dict['HASC_1']]=state_color
                else: continue#colors[shape_dict['HASC_1']]=0
            ax = plt.gca() # get current axes instance
            for nshape, seg in enumerate(m.states):
                if names[nshape] not in colors: continue
                poly = Polygon(seg,facecolor=mapper.to_rgba(colors[names[nshape]]), edgecolor='k')
                ax.add_patch(poly)
            mapper.set_array(clrBr) 
        
        
        cax = fig.add_axes([0.27, 0.1, 0.5, 0.05]) # posititon
        cb = ColorbarBase(cax,cmap=thisCMAP,norm=norm,orientation='horizontal')
        mapper.set_array(clrBr);
        fig.set_size_inches(30,30)
        plt.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/EOF/NAOw_'+crop+'_EOF'+str(PC+1)+ENSOyr+str(PCs)+'PCs'+rot+notes+'.png')
        plt.close()