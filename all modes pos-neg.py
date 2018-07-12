#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 10:02:28 2018
@author: weston

NOTES: PCs were coded to lag early spring backward. So loading in A-M-M-J-J-A-S-O-N-D 
        in PC-year 1998 will be calendar year 1998, but in J-F-M in 
        PC-year 1998 will be calendar year 1999

EDIT - 1/29/18 - change the data source to be the post-processed data in the data folder
                instead of the data in the project files
               - Include precip, geopotential height, winds and anomalous ascent
"""

import numpy as np
import netCDF4
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import mpl_toolkits
plt.ioff()

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#
#        Script buttons and knobs       #
years = [1981,2011] #from the SVD analysis
saveNotes = ''#'_ensoLC'
modes=['ENSO','IOD','TAV (AMM)','NAO']
thresh = .65
#       End of buttons and knobs        #
#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#

#ENSO life cycle years
#'EN': [1982,1987,1997,2009,2006,1994]#2006,1994, 'LN':[1983,1988,1998,2010,2007,1995]#2007,1995
#EN':[1997-1,1982-1,2009-1,1987-1],'LN':

#strong years
#'EN':[1982,1997,2009,2002,1991,1987], 'LN':[1988,1998,2007,1983,1995,2010]


posModeYrs = {'ENSO':[1982-1,1997-1,2009-1,2002-1],'NAO':[1998-1,1994-1,1983-1,1988-1],
              'IOD':[1992-1,1996-1,1981-1,1984-1],'TAV (AMM)':[2010-1,2005-1,1981-1,1997-1]}

negModeYrs = {'ENSO':[1988-1,1998-1,2007-1,1995-1],'NAO':[2009-1,1997-1,1984-1,1986-1],
              'IOD':[1994-1,1997-1,1982-1,1987-1],'TAV (AMM)':[1991-1,1986-1,1994-1,2002-1]}

corners = {'ENSO':[-60,60,-179,179],'NAO':[20,80,-50,50],'TAV (AMM)':[-60,30,-80,20],'IOD':[-40,30,20,160]}

#Read in agreement mask. agreeMSK(DF, frac, require agree w/ avg? (Y/N) )
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/agreement.py').read())
exec(open('/Users/weston/Desktop/Columbia/Research/Code/PyCode/general/helper_Functions/custom_div_cmap.py').read())

yrList = np.array(range(years[0],years[1]+1))
seasons = ['DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND','NDJ']
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ENSOyr = 'Jan-DecENSOyr0har'


#read in the netCDF file containing the variables
ncGeo = netCDF4.Dataset('/Volumes/Data_Archive/Data/NCEP_monthly_nc/detrended/gph.nc',
                        'r',format='NETCDF4') #read in the netCDF file for the variable
geoVar=ncGeo.variables['hgt'][(ncGeo.variables['year'][:]>=years[0])&(ncGeo.variables['year'][:]<years[1]+2),
                              :,np.where(ncGeo.variables['level'][:]==200)[0][0],:,:]

ncSM = netCDF4.Dataset('/Volumes/Data_Archive/Data/SoilMoisture/GLDAS/Noah10_nc/detrended/SoilM.nc','r',format='NETCDF4') 
sm = np.nansum(ncSM.variables['SoilM'][(ncSM.variables['year'][:]>=years[0])&(ncSM.variables['year'][:]<years[1]+2),:,:3,:,:],2)

ncMaxT = netCDF4.Dataset('//Volumes/Data_Archive/Data/Temp/BerkeleyEarth/monthly/Complete_TMAX_LatLong1.nc','r',format='NETCDF4') #read in the netCDF file for the variable
maxT = ncMaxT.variables['temperature'][(ncMaxT.variables['time'][:]>=years[0])&(ncMaxT.variables['time'][:]<years[1]+2),:,:] 

ncUwnd = netCDF4.Dataset('/Volumes/Data_Archive/Data/Wind/detrend/uwnd.nc','r',format='NETCDF4') #read in the netCDF file for the variable
Uwnd=ncUwnd.variables['uwnd'][(ncUwnd.variables['year'][:]>=years[0])&(ncUwnd.variables['year'][:]<years[1]+2),
                              :,np.where(ncUwnd.variables['level'][:]==925)[0][0],:,:]
ncVwnd = netCDF4.Dataset('/Volumes/Data_Archive/Data/Wind/detrend/vwnd.nc','r',format='NETCDF4') #read in the netCDF file for the variable
Vwnd=ncVwnd.variables['vwnd'][(ncVwnd.variables['year'][:]>=years[0])&(ncVwnd.variables['year'][:]<years[1]+2),
                              :,np.where(ncVwnd.variables['level'][:]==925)[0][0],:,:]



#Read in the SST anomaly data
sstData = netCDF4.Dataset('/Volumes/Data_Archive/Data/Sea_Surface_Temp/ERSSTv3b_anom.cdf')
#T is months since 1960-01-01, and I want to start at 1950 with one year lag
sstLats=sstData.variables['Y'][:]
sstLons=sstData.variables['X'][:]
sstLons,sstLats=np.meshgrid(sstLons,sstLats)

sstAnom=sstData.variables['anom'][(sstData.variables['T'][:]>=-12*(1960-years[0]))&
                        (sstData.variables['T'][:]<=12*(2+years[1]-1960)),...]
sstAnom = np.squeeze(sstAnom);sstAnom=np.ma.filled(sstAnom,0)
sstAnom_M = np.squeeze(sstAnom);sstAnom=np.ma.filled(sstAnom,0)


#Remove the climatology to get anomalies by month
for m in range(0,12):
    geoVar[:,m,...] = geoVar[:,m,...] - np.nanmean(geoVar[:,m,...],axis=0)
    sm[:,m,...] = sm[:,m,...] - np.nanmean(sm[:,m,...],axis=0)
    Uwnd[:,m,...] = Uwnd[:,m,...] - np.nanmean(Uwnd[:,m,...],axis=0)
    Vwnd[:,m,...] = Vwnd[:,m,...] - np.nanmean(Vwnd[:,m,...],axis=0)
maxT[:,...] = maxT[:,...] - np.nanmean(maxT[:,...],axis=0)

#Rearrange to have only one time dimension
geoVar_M = geoVar.reshape([geoVar.shape[0]*geoVar.shape[1],geoVar.shape[2],geoVar.shape[3]])
sm_M = sm.reshape([sm.shape[0]*sm.shape[1],sm.shape[2],sm.shape[3]])
maxT_M = maxT#.reshape([maxT.shape[0]*maxT.shape[1],maxT.shape[2],maxT.shape[3]])
Uwnd_M = Uwnd.reshape([Uwnd.shape[0]*Uwnd.shape[1],Uwnd.shape[2],Uwnd.shape[3]])
Vwnd_M = Vwnd.reshape([Vwnd.shape[0]*Vwnd.shape[1],Vwnd.shape[2],Vwnd.shape[3]])
          
NCyrs =ncGeo.variables['year'][(ncGeo.variables['year'][:]>=years[0])&(ncGeo.variables['year'][:]<years[1]+2)]
gphLats=ncGeo.variables['latitude'][:]
gphLons=ncGeo.variables['longitude'][:]#-180.
gphLons, gphLats = np.meshgrid(gphLons,gphLats)

mLats=ncMaxT.variables['latitude'][:]
mLons=ncMaxT.variables['longitude'][:]#-180.
mLons, mLats = np.meshgrid(mLons,mLats)
smLats=ncSM.variables['latitude'][:]
smLons=ncSM.variables['longitude'][:]#-180.
smLons, smLats = np.meshgrid(smLons,smLats)


for iS in [-5,-3,-1,0,1,3,5,7,11,12,13,15]:
    for mode in modes:
        posCompYrs = posModeYrs[mode]
        negCompYrs = negModeYrs[mode]
        posYrInd = np.array([iY == NCyrs for iY in posCompYrs])
        posYrInd = np.where(posYrInd.sum(0)==1)[0]
        
        negYrInd = np.array([iY == NCyrs for iY in negCompYrs])
        negYrInd = np.where(negYrInd.sum(0)==1)[0]
        yrOffset = 0 + (iS//12)
        posYrInd = np.tile(np.arange((iS%12-1),(iS%12+2)),posYrInd.shape[0])+ \
                  np.repeat(posYrInd+yrOffset+1,3)*12
        negYrInd = np.tile(np.arange((iS%12-1),(iS%12+2)),negYrInd.shape[0])+ \
                  np.repeat(negYrInd+yrOffset+1,3)*12   
        limPos = np.int(posYrInd.shape[0]-np.ceil(np.sum(posYrInd>=sm_M.shape[0])/3)*3)
        limNeg = np.int(negYrInd.shape[0]-np.ceil(np.sum(negYrInd>=sm_M.shape[0])/3)*3)
                  
        #index
        geoVar = geoVar_M[posYrInd,...]
        posGPH = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
        geoVar = geoVar_M[negYrInd,...]
        negGPH = np.nanmean(geoVar.reshape([geoVar.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    
        Uwnd = Uwnd_M[posYrInd,...]
        posUwnd = np.nanmean(Uwnd.reshape([Uwnd.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
        Uwnd = Uwnd_M[negYrInd,...]
        negUwnd = np.nanmean(Uwnd.reshape([Uwnd.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    
        Vwnd = Vwnd_M[posYrInd,...]
        posVwnd = np.nanmean(Vwnd.reshape([Vwnd.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
        Vwnd = Vwnd_M[negYrInd,...]
        negVwnd = np.nanmean(Vwnd.reshape([Vwnd.shape[0]//3,3,gphLons.shape[0],gphLons.shape[1]]),1)
    
        sm = sm_M[posYrInd[:limPos],...]
        posSM = np.nanmean(sm.reshape([sm.shape[0]//3,3,smLons.shape[0],smLons.shape[1]]),1)
        sm = sm_M[negYrInd[:limNeg],...]
        negSM = np.nanmean(sm.reshape([sm.shape[0]//3,3,smLons.shape[0],smLons.shape[1]]),1)
    
        maxT = maxT_M[posYrInd,...]
        posMaxT = np.nanmean(maxT.reshape([maxT.shape[0]//3,3,mLons.shape[0],mLons.shape[1]]),1)
        maxT = maxT_M[negYrInd,...]
        negMaxT = np.nanmean(maxT.reshape([maxT.shape[0]//3,3,mLons.shape[0],mLons.shape[1]]),1)
    
        sstAnom = sstAnom_M[posYrInd,...]
        posSST = np.nanmean(sstAnom.reshape([sstAnom.shape[0]//3,3,sstLons.shape[0],sstLons.shape[1]]),1)
        sstAnom = sstAnom_M[negYrInd,...]
        negSST = np.nanmean(sstAnom.reshape([sstAnom.shape[0]//3,3,sstLons.shape[0],sstLons.shape[1]]),1)
                  
        posGPH_all = np.nanmean(posGPH,0);posGPH[posGPH==0]=np.nan
        negGPH_all = np.nanmean(negGPH,0);negGPH[negGPH==0]=np.nan
        posGPH = np.nanmean(posGPH,0)#*agreeMSK(posGPH,thresh,'N');posGPH[posGPH==0]=np.nan
        negGPH = np.nanmean(negGPH,0)#*agreeMSK(negGPH,thresh,'N');negGPH[negGPH==0]=np.nan   

        posUwnd_all = np.nanmean(posUwnd,0);
        negUwnd_all = np.nanmean(negUwnd,0);
        posUwnd = np.nanmean(posUwnd,0)#*agreeMSK(posUwnd,thresh,'N');posUwnd[posUwnd==0]=np.nan
        negUwnd = np.nanmean(negUwnd,0)#*agreeMSK(negUwnd,thresh,'N');negUwnd[negUwnd==0]=np.nan

        posVwnd_all = np.nanmean(posVwnd,0);
        negVwnd_all = np.nanmean(negVwnd,0);
        posVwnd = np.nanmean(posVwnd,0)#*agreeMSK(posVwnd,thresh,'N');posVwnd[posVwnd==0]=np.nan
        negVwnd = np.nanmean(negVwnd,0)#*agreeMSK(negVwnd,thresh,'N');negVwnd[negVwnd==0]=np.nan

        posSST = np.nanmean(posSST,0)#*agreeMSK(posSST,thresh,'N');posSST[posSST==0]=np.nan
        negSST = np.nanmean(negSST,0)#*agreeMSK(negSST,thresh,'N');negSST[negSST==0]=np.nan        
        posMaxT = np.nanmean(posMaxT,0)#*agreeMSK(posMaxT,thresh,'N');posMaxT[posMaxT==0]=np.nan
        negMaxT = np.nanmean(negMaxT,0)#*agreeMSK(negMaxT,thresh,'N');negMaxT[negMaxT==0]=np.nan        
        posSM = np.nanmean(posSM,0)#*agreeMSK(posSM,thresh,'N');posSM[posSM==0]=np.nan
        negSM = np.nanmean(negSM,0)#*agreeMSK(negSM,thresh,'N');negSM[negSM==0]=np.nan
    
        sstDiff = posSST-negSST; sstDiff = np.append(sstDiff[:,90:],sstDiff[:,:90],1)
        smDiff = posSM-negSM; smDiff = np.append(smDiff[:,90:],smDiff[:,:90],1)
        maxtDiff = posMaxT-negMaxT; maxtDiff = np.append(maxtDiff[:,90:],maxtDiff[:,:90],1)
        
        #define color scales 
        hw=2.5
        gphMax = 200
        gphLvls = np.arange(-gphMax,gphMax+10,15) 
        tScale = 3
        qScale = 10
        latSpc=20
        lonSpc=30
        if (mode == 'IOD')|(mode == 'TAV (AMM)')|(mode=='NAO'):
            tScale = 2
            qScale = 10
            latSpc=10
            lonSpc=20
            hw=4
            if (mode == 'TAV (AMM)'):
                qScale = 5
            elif (mode=='NAO'):
                hw=3
                latSpc=10
                lonSpc=10
                qscale=15
                gphLvls = np.arange(-gphMax,gphMax+10,30) 
            

        pScale = 5
        smScale = 40
        omegaScale = .025
        
        
        fig = plt.figure();
        ax11=plt.subplot(111)
        
        m13 = Basemap(projection='cyl',llcrnrlat=corners[mode][0],urcrnrlat=corners[mode][1],
                      llcrnrlon=corners[mode][2],urcrnrlon=corners[mode][3],resolution='c',ax=ax11)
        m13.drawcoastlines();m13.drawcountries()

        im13=plt.contourf(sstLons-180, sstLats,sstDiff,shading='flat',cmap='RdBu_r',latlon=True,levels=np.arange(-tScale,tScale+tScale/5,tScale/5),extend='both')  
#        if mode is 'NAO':
#            im13=plt.contourf(mLons, mLats,posMaxT-negMaxT,shading='flat',cmap='PuOr_r',latlon=True,levels=np.arange(-tScale,tScale+tScale/5,tScale/5) ,extend='both')  
#        else:
        im13=plt.contourf(smLons, smLats,posSM-negSM,shading='flat',cmap='BrBG',latlon=True,levels=np.arange(-smScale,smScale+smScale/5,smScale/5),extend='both')  
        cb13 = m13.colorbar(im13,"right", size="5%", pad=2);
        cb13.ax.set_yticklabels(cb13.ax.get_yticklabels(), fontsize=24)
        
        im13=m13.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd_all[::2,::2]-negUwnd_all[::2,::2],posVwnd[::2,::2]-negVwnd[::2,::2],
                        latlon=True, color='grey',pivot='middle', units='inches',scale=qScale,headwidth=hw,width=0.03)
        im13=m13.quiver(gphLons[::2,::2],gphLats[::2,::2],posUwnd[::2,::2]-negUwnd[::2,::2],posVwnd[::2,::2]-negVwnd[::2,::2],
                        latlon=True, pivot='middle', units='inches',scale=qScale,headwidth=hw,width=0.03)

        if mode == 'NAO':
            plt.quiverkey(im13,0.5,1.05,5,'5 m/s',fontproperties={'size':24})
        else:
            plt.quiverkey(im13,0.5,1.05,5,'5 m/s',fontproperties={'size':24})
        plt.text(corners[mode][3]-5,corners[mode][1]+5,'GPH lvls: '+str(gphLvls[1]-gphLvls[0]),fontsize=24)


        im13=m13.contour(gphLons, gphLats,posGPH_all-negGPH_all,gphLvls, shading='flat',colors='darkgrey',latlon=True)   
        im13=m13.contour(gphLons, gphLats,posGPH-negGPH,gphLvls, shading='flat',colors='k',latlon=True) 
        m13.drawparallels(np.arange(int(-95),int(95),latSpc),labels=[0,1,0,0],fontsize=24,linewidth=0.0)
        m13.drawmeridians(np.arange(int(-185),int(185),lonSpc),labels=[0,0,0,1], fontsize=24,linewidth=0.0)
        
        fig.set_size_inches(25,10);
        fig.savefig('/Users/weston/Desktop/Columbia/Research/Results/Climate Modes - Yield Anoms/SVD climate/event composites/'+mode+'_'+ENSOyr+'_'+seasons[iS%12]+'_ENSOyrs_yrOffset'+str(yrOffset)+saveNotes+'.png')
        plt.close()        
    
    
