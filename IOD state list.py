#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 14:09:27 2018

@author: weston
"""



wIODstates_fewer = ['Shanghai', 'Jiangsu', 'Zhejiang', 'Anhui', 'Fujian','Shaanxi', 
              'Jiangxi',  'Henan', 'Hubei', 'Hunan', 'Guangdong', 
              'Guangxi', 'Hainan', 'Chongqing', 'Sichuan', 'Guizhou', 'Yunnan',
              'New South Wales(b)', 'Victoria', 'Queensland', 'South Australia',
              'Western Australia', 'Tasmania', 'Northern Territory', 'Australian Capital Territory', 
              'Andhra Pradesh', 'Arunachal Pradesh', 'Assam', 'Bihar', 'Gujarat', 
              'Haryana', 'Himachal Pradesh', 'Jammu and Kashmir', 'Karnataka', 
              'Madhya Pradesh', 'Maharashtra', 'Manipur', 'Meghalaya', 'Nagaland', 
              'Odisha', 'Punjab', 'Rajasthan', 'Sikkim', 'Tamil Nadu', 'Tripura',
              'Uttar Pradesh', 'West Bengal', 'D & N Haveli', 'Delhi',
              'Afghanistan','Bangladesh','Bhutan','Botswana','Cambodia','Eritrea', 
              'Ethiopia', 'Kenya','Madagascar', 'Malawi', 'Mozambique', 'Myanmar', 
              'Nepal', 'New Zealand', 'Somalia', 'South Africa', 
              'South Sudan', 'Sudan', 'Swaziland', 'Tanzania','Zambia', 'Zimbabwe']#,

#mIODstates  = ['none']
#sIODstates  = ['none']

mIODstates_fewer = ['Shanghai', 'Jiangsu', 'Zhejiang', 'Anhui', 'Fujian','Shaanxi',
              'Jiangxi',  'Henan', 'Hubei', 'Hunan', 'Guangdong', 
              'Guangxi', 'Hainan', 'Chongqing', 'Sichuan', 'Guizhou', 'Yunnan',
              'Andhra Pradesh', 'Arunachal Pradesh', 'Assam', 'Gujarat', 'Haryana', 
              'Himachal Pradesh', 'Jammu and Kashmir', 'Karnataka', 'Madhya Pradesh',
              'Maharashtra', 'Manipur', 'Meghalaya', 'Mizoram', 'Nagaland', 'Odisha',
              'Punjab', 'Rajasthan', 'Sikkim', 'Tamil Nadu', 'Tripura', 'Uttar Pradesh', 
              'West Bengal', 'Australia', 'Bangladesh', 'Bhutan', 
              'Botswana','Cambodia', 'Djibouti','Eritrea', 'Ethiopia', 'Indonesia', 
              'Kenya', 'Kuwait', "Lao People's Democratic Republic", 'Madagascar', 'Malawi', 
              'Malaysia',  'Micronesia (Federated States of)', 'Mozambique',
              'Myanmar', 'Nepal', 'New Zealand', 'Philippines','Somalia', 'South Africa', 
              'South Sudan',  'Sri Lanka', 'Sudan', 'Swaziland', 'Tanzania', 
              'Thailand', 'Timor-Leste','Vietnam', 'Zambia', 'Zimbabwe']#,

sIODstates_fewer  = ['Shanghai', 'Jiangsu', 'Zhejiang', 'Anhui', 'Fujian','Shaanxi',
              'Jiangxi',  'Henan', 'Hubei', 'Hunan', 'Guangdong', 
              'Guangxi', 'Hainan', 'Chongqing', 'Sichuan', 'Guizhou', 'Yunnan',
              'Australia', 'Bangladesh', 'Bhutan',  'Cambodia', 'Ethiopia', 'India', 'Indonesia', 
              'Kenya', "Lao People's Democratic Republic",  'Madagascar', 'Malawi', 
              'Malaysia', 'Myanmar', 'Nepal',  'New Zealand',
              'Philippines', 'South Africa',  'Sri Lanka',  'Tanzania',
              'Thailand', 'Timor-Leste', 'Vietnam',  'Zambia', 'Zimbabwe']#,



wIODstates = ['Tianjin', 'Hebei', 'Shanxi', 'Nei Mongol', 'Liaoning', 'Jilin',
              'Heilongjiang', 'Shanghai', 'Jiangsu', 'Zhejiang', 'Anhui', 'Fujian',
              'Jiangxi', 'Shandong', 'Henan', 'Hubei', 'Hunan', 'Guangdong', 
              'Guangxi', 'Hainan', 'Chongqing', 'Sichuan', 'Guizhou', 'Yunnan',
              'Tibet', 'Shaanxi', 'Gansu', 'Qinghai', 'Ningxia Hui', 'Xinjiang',
              'New South Wales', 'Victoria', 'Queensland', 'South Australia',
              'Western Australia', 'Tasmania', 'Northern Territory', 'Australian Capital Territory', 
              'Andhra Pradesh', 'Arunachal Pradesh', 'Assam', 'Bihar', 'Gujarat', 
              'Haryana', 'Himachal Pradesh', 'Jammu and Kashmir', 'Karnataka', 
              'Madhya Pradesh', 'Maharashtra', 'Manipur', 'Meghalaya', 'Nagaland', 
              'Odisha', 'Punjab', 'Rajasthan', 'Sikkim', 'Tamil Nadu', 'Tripura',
              'Uttar Pradesh', 'West Bengal', 'D & N Haveli', 'Delhi',
              'Afghanistan','Bangladesh','Bhutan','Botswana','Cambodia','Eritrea', 
              'Ethiopia','Iran', 'Iraq','Israel','Japan','Jordan', 'Kenya', 'Kuwait', 
              'Lebanon', 'Madagascar', 'Malawi', 'Mozambique', 'Myanmar', 
              'Nepal', 'New Zealand', 'Occupied Palestinian Territory', 
              'Oman', 'Pakistan','Qatar', 'Saudi Arabia', 'Somalia', 'South Africa', 
              'South Sudan', 'Sudan', 'Swaziland', 'Syria', 'Tanzania', 
              'United Arab Emirates', 'Yemen', 'Zambia', 'Zimbabwe']#,

mIODstates  = ['Tianjin', 'Hebei', 'Shanxi', 'Nei Mongol', 'Liaoning', 'Jilin', 
              'Heilongjiang', 'Shanghai', 'Jiangsu', 'Zhejiang', 'Anhui', 'Fujian',
              'Jiangxi', 'Shandong', 'Henan', 'Hubei', 'Hunan', 'Guangdong', 
              'Guangxi', 'Hainan', 'Chongqing', 'Sichuan', 'Guizhou', 'Yunnan',
              'Tibet', 'Shaanxi', 'Gansu', 'Qinghai', 'Ningxia Hui', 'Xinjiang',
              'Andhra Pradesh', 'Arunachal Pradesh', 'Assam', 'Gujarat', 'Haryana', 
              'Himachal Pradesh', 'Jammu and Kashmir', 'Karnataka', 'Madhya Pradesh',
              'Maharashtra', 'Manipur', 'Meghalaya', 'Mizoram', 'Nagaland', 'Odisha',
              'Punjab', 'Rajasthan', 'Sikkim', 'Tamil Nadu', 'Tripura', 'Uttar Pradesh', 
              'West Bengal','Afghanistan', 'Australia', 'Bangladesh', 'Bhutan', 
              'Botswana','Cambodia', 'Djibouti','Eritrea', 'Ethiopia', 'Indonesia', 
              'Iran', 'Iraq', 'Israel', 'Japan', 'Jordan', 'Kenya', 'Kuwait', 
              "Lao People's Democratic Republic", 'Lebanon', 'Madagascar', 'Malawi', 
              'Malaysia',  'Micronesia (Federated States of)', 'Mozambique',
              'Myanmar', 'Nepal', 'New Zealand', 'Oman',  'Pakistan', 
              'Philippines','Qatar','Saudi Arabia','Somalia', 'South Africa', 
              'South Sudan',  'Sri Lanka', 'Sudan', 'Swaziland','Syria', 'Tanzania', 
              'Thailand', 'Timor-Leste','United Arab Emirates','Vietnam', 'Yemen', 'Zambia', 'Zimbabwe']#,

sIODstates = ['Tianjin', 'Hebei', 'Shanxi', 'Nei Mongol', 'Liaoning', 'Jilin', 
              'Heilongjiang', 'Shanghai', 'Jiangsu', 'Zhejiang', 'Anhui', 'Fujian',
              'Jiangxi', 'Shandong', 'Henan', 'Hubei', 'Hunan', 'Guangdong', 
              'Guangxi', 'Hainan', 'Chongqing', 'Sichuan', 'Guizhou', 'Yunnan',
              'Tibet', 'Shaanxi', 'Gansu', 'Qinghai', 'Ningxia Hui', 'Xinjiang', 
              'Australia', 'Bangladesh', 'Bhutan',  'Cambodia', 
              'Ethiopia', 'India', 'Indonesia', 'Iran', 'Iraq', 'Japan', 'Jordan',
              'Kenya', "Lao People's Democratic Republic",  'Madagascar', 'Malawi', 
              'Malaysia', 'Myanmar', 'Nepal',  'New Zealand','Pakistan', 
              'Philippines', 'South Africa',  'Sri Lanka',  'Syria', 'Tanzania',
              'Thailand', 'Timor-Leste', 'Vietnam',  'Zambia', 'Zimbabwe']#,

