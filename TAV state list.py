#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 14:09:27 2018

@author: weston
"""

TAVstates_least = ['Angola',  'Benin','Burkina Faso', 'Cameroon', "Côte d'Ivoire", 
       'Eq. Guinea','Gabon', 'Gambia', 'Ghana', 'Guinea', 'Guinea-Bissau', 'Guyana',
       'Liberia', 'Namibia',  'Nigeria', 'Senegal','Sierra Leone', 'Singapore',
       'Suriname', 'Togo', 'Tonga', 'Trinidad and Tobago','Venezuela', 'W. Sahara',
       'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA','MG','MT','PA',
       'PB','PE','PI','RJ','RN','RO','RR','SE','TO']#,

TAVstates_fewer = ['American Samoa', 'Andorra','Angola', 'Barbados','Belize', 'Benin',
       'Bermuda', 'Burkina Faso', 'Cameroon', 'Central African Rep.', 'Chad',
       'Colombia', 'Congo','Costa Rica',  'Cuba',  "Côte d'Ivoire", 'Dem. Rep. Congo',
       'Dominica', 'Dominican Rep.',  'El Salvador',
       'Eq. Guinea','Gabon', 'Gambia', 'Ghana', 'Guam','Guatemala', 'Guinea', 'Guinea-Bissau', 'Guyana',
       'Haiti', 'Honduras','Jamaica', 'Liberia','Mali','Namibia','Nicaragua',
       'Niger', 'Nigeria', 'Panama', 'Paraguay', 'Puerto Rico','Senegal',
       'Sierra Leone', 'Singapore', 'Suriname', 'Togo', 'Tonga', 'Trinidad and Tobago',
       'Uruguay','Venezuela', 'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA','MG','MS','MT','PA',
       'PB','PE','PI','PR','RJ','RN','RO','RR','RS','SC','SE','SP','TO',
       'Campeche', 'Chiapas',  'Hidalgo','Morelos',  
       'Nuevo León', 'Oaxaca', 'Puebla', 'Querétaro', 'Quintana Roo',
       'San Luis Potosí',  'Tabasco', 'Tlaxcala', 'Veracruz', 'Yucatan']#,

TAVstates_some = ['Angola', 'Belize', 'Benin',
       'Burkina Faso', 'Cameroon', 'Congo',
       'Colombia', 'Costa Rica',  'Cuba',  "Côte d'Ivoire", 
       'Dominica', 'Dominican Rep.', 'El Salvador',
       'Eq. Guinea','Gabon', 'Gambia', 'Ghana', 'Guam','Guatemala', 'Guinea', 'Guinea-Bissau', 'Guyana',
       'Haiti', 'Honduras','Jamaica',  'Liberia', 
       'Mexico','Namibia','Nicaragua',  'Nigeria', 'Mali',
       'Panama', 'Paraguay', 'Puerto Rico','Senegal','Sierra Leone', 'Singapore', 'Suriname', 
       'Togo', 'Tonga', 'Trinidad and Tobago','Venezuela', 'W. Sahara',
       'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA','MG','MS','MT','PA',
       'PB','PE','PI','PR','RJ','RN','RO','RR','RS','SC','SE','SP','TO',
       'Aguascalientes', 'Baja California', 'Baja California Sur',
       'Campeche', 'Chiapas', 'Chihuahua', 'Coahuila', 'Colima',
       'Distrito Federal', 'Durango', 'Guanajuato', 'Guerrero', 'Hidalgo',
       'Jalisco', 'Michoacán', 'Morelos', 'México', 'Nayarit',
       'Nuevo León', 'Oaxaca', 'Puebla', 'Querétaro', 'Quintana Roo',
       'San Luis Potosí', 'Sinaloa', 'Sonora', 'Tabasco', 'Tamaulipas',
       'Tlaxcala', 'Veracruz', 'Yucatán', 'Zacatecas']#,

TAVstates = ['American Samoa', 'Andorra','Angola', 'Barbados','Belize', 'Benin',
       'Bermuda', 'Burkina Faso', 'Cameroon', 'Central African Rep.', 'Chad',
       'Colombia', 'Congo','Costa Rica',  'Cuba',  "Côte d'Ivoire", 'Dem. Rep. Congo',
       'Dominica', 'Dominican Rep.',  'El Salvador',
       'Eq. Guinea','Gabon', 'Gambia', 'Ghana', 'Guam','Guatemala', 'Guinea', 'Guinea-Bissau', 'Guyana',
       'Haiti', 'Honduras','Jamaica', 'Liberia','Mali','Namibia','Nicaragua',
       'Niger', 'Nigeria', 'Panama', 'Paraguay', 'Puerto Rico','Senegal',
       'Sierra Leone', 'Singapore', 'Suriname', 'Togo', 'Tonga', 'Trinidad and Tobago',
       'Uruguay','Venezuela', 'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA','MG','MS','MT','PA',
       'PB','PE','PI','PR','RJ','RN','RO','RR','RS','SC','SE','SP','TO',
       'Buenos Aires', 'Catamarca', 'Chaco', 'Chubut',
       'Ciudad de Buenos Aires', 'Corrientes', 'Cordoba', 'Entre Rios',
       'Formosa', 'Jujuy', 'La Pampa', 'La Rioja', 'Mendoza', 'Misiones',
       'Neuquén', 'Río Negro', 'Salta', 'San Juan', 'San Luis',
       'Santa Cruz', 'Santa Fe', 'Santiago del Estero', 'Tierra del Fuego',
       'Tucumán','Campeche', 'Chiapas',  'Hidalgo','Morelos',  
       'Nuevo León', 'Oaxaca', 'Puebla', 'Querétaro', 'Quintana Roo',
       'San Luis Potosí',  'Tabasco', 'Tlaxcala', 'Veracruz', 'Yucatan']#,

TAVstates_wUS = ['ALABAMA', 'ARKANSAS', 'FLORIDA', 'GEORGIA',  'LOUISIANA', 'MARYLAND', 
             'MISSISSIPPI', 'NORTH CAROLINA','SOUTH CAROLINA', 'TEXAS',
        'American Samoa', 'Andorra','Angola', 'Barbados','Belize', 'Benin',
       'Bermuda', 'Burkina Faso', 'Cameroon', 'Central African Rep.', 'Chad',
       'Colombia', 'Congo','Costa Rica',  'Cuba',  "Côte d'Ivoire", 'Dem. Rep. Congo',
       'Dominica', 'Dominican Rep.',  'El Salvador',
       'Eq. Guinea','Gabon', 'Gambia', 'Ghana', 'Guam','Guatemala', 'Guinea', 'Guinea-Bissau', 'Guyana',
       'Haiti', 'Honduras','Jamaica', 'Liberia', 'Mali','Mexico','Namibia','Nicaragua',
       'Niger', 'Nigeria', 'Panama', 'Paraguay', 'Puerto Rico','Senegal',
       'Sierra Leone', 'Singapore', 'Suriname', 'Togo', 'Tonga', 'Trinidad and Tobago',
       'Uruguay','Venezuela', 'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA','MG','MS','MT','PA',
       'PB','PE','PI','PR','RJ','RN','RO','RR','RS','SC','SE','SP','TO',
       'Buenos Aires', 'Catamarca', 'Chaco', 'Chubut',
       'Ciudad de Buenos Aires', 'Corrientes', 'Cordoba', 'Entre Rios',
       'Formosa', 'Jujuy', 'La Pampa', 'La Rioja', 'Mendoza', 'Misiones',
       'Neuquén', 'Río Negro', 'Salta', 'San Juan', 'San Luis',
       'Santa Cruz', 'Santa Fe', 'Santiago del Estero', 'Tierra del Fuego',
       'Tucumán','Aguascalientes', 'Baja California', 'Baja California Sur',
       'Campeche', 'Chiapas', 'Chihuahua', 'Coahuila', 'Colima',
       'Distrito Federal', 'Durango', 'Guanajuato', 'Guerrero', 'Hidalgo',
       'Jalisco', 'Michoacán', 'Morelos', 'México', 'Nayarit',
       'Nuevo León', 'Oaxaca', 'Puebla', 'Querétaro', 'Quintana Roo',
       'San Luis Potosí', 'Sinaloa', 'Sonora', 'Tabasco', 'Tamaulipas',
       'Tlaxcala', 'Veracruz', 'Yucatán', 'Zacatecas']#,

