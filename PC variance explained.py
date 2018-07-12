#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 07:57:02 2018

@author: weston
"""
import numpy as np
import matplotlib.pyplot as plt

years = [1981,2011]
timePeriod =''# ' 1960-2010' #'' or ' 1960-2010'
crop = 'maize'
ENSOyr = 'May-AprENSOyr15har'


cPC = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/'+crop+'/crop_'+crop+
              '_'+ENSOyr+'_3PCs_'+str(years[0])+'-'+str(years[1])+'_v2.npy') 
sstPC = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/'+crop+'/sst_'+crop+
              '_'+ENSOyr+'_3PCs_'+str(years[0])+'-'+str(years[1])+'_v2.npy')             

eigens = np.load('/Volumes/Data_Archive/Results/Climate Modes - Yield Anoms/SVD'+timePeriod+'/'+crop+'/sst_'+crop+
              '_'+ENSOyr+'_eigens_'+str(years[0])+'-'+str(years[1])+'_v2.npy')             
