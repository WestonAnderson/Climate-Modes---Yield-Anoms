#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 08:23:09 2018

@author: weston
"""
mRegions = {   
'US Midwest': ['IOWA','ILLINOIS','MISSOURI','INDIANA','KENTUCKY',
       'OHIO','SOUTH DAKOTA','NEBRASKA','KANSAS','NORTH DAKOTA','WISCONSIN','MICHIGAN','MINNESOTA'],#

'US Corn Belt': ['IOWA','ILLINOIS','MISSOURI','INDIANA','KENTUCKY','OHIO','WISCONSIN','MICHIGAN','MINNESOTA'],#
               
               
'Northeast Brazil': ['MA','TO','PA','PI','BA','CE','RN','PB','PE','AL','SE'],

'West Africa':["Côte d'Ivoire",'Senegal','Gambia','Guinea','Guinea-Bissau','Mali',
'Sierra Leone','Ghana','Liberia','Togo','Benin','Burkina Faso','Nigeria','Cameroon'],

'Southeast Africa': ['South Africa','Zimbabwe','Malawi','Mozambique','Zambia','Botswana'],#'Tanzania'

'East Africa':['Ethiopia','Kenya','Somalia','Somaliland','Eritrea','Djibouti','Uganda'],

'West India': ['Rajasthan','Gujarat','Madhya Pradesh','Maharashtra','Andhra Pradesh',
        'Karnataka'], #'Odisha'?

'India': ['Andaman and Nicobar', 'Andhra Pradesh', 'Arunachal Pradesh',
       'Assam', 'Bihar', 'Chandigarh', 'Chhattisgarh',
       'Dadra and Nagar Haveli', 'Daman and Diu', 'Goa', 'Gujarat',
       'Haryana', 'Himachal Pradesh', 'Jammu and Kashmir', 'Jharkhand',
       'Karnataka', 'Kerala', 'Lakshadweep', 'Madhya Pradesh',
       'Maharashtra', 'Manipur', 'Meghalaya', 'Mizoram', 'NCT of Delhi',
       'Nagaland', 'Odisha', 'Puducherry', 'Punjab', 'Rajasthan', 'Sikkim',
       'Tamil Nadu', 'Telangana', 'Tripura', 'Uttar Pradesh',
       'Uttarakhand', 'West Bengal'],
               
'Southeast South America':['Uruguay','Buenos Aires','Entre Rios','La Pampa','San Luis','Santiago del Estero',
        'Cordoba','Tucuman','Chaco','Corrientes','PR','SC','RS'],

'Brazil-S.BR':['RO','AC','AM','AP','PA','TO','MA','PI','CE','RN','PB','PE',
        'AL','SE','BA','MT','MS','GO','DF','MG','ES','RJ','SP','PR','SC','RS'],
                           
'Southwest Mexico':['Colima','Jalisco','Guerrero','Oaxaca','Michoacan'],
                    #,'Puebla','Morelos','Mexico','Guanajuato','Queretaro','Chiapas','Puebla'

'United States': ['ALABAMA','ARIZONA','ARKANSAS','CALIFORNIA','COLORADO','DELAWARE','FLORIDA','GEORGIA',
 'IDAHO', 'ILLINOIS', 'INDIANA', 'IOWA', 'KANSAS', 'KENTUCKY', 'LOUISIANA', 'MARYLAND', 'MICHIGAN',
 'MINNESOTA', 'MISSISSIPPI', 'MISSOURI', 'MONTANA', 'NEBRASKA', 'NEW JERSEY', 'NEW MEXICO', 'NEW YORK',
 'NORTH CAROLINA', 'NORTH DAKOTA', 'OHIO', 'OKLAHOMA', 'OREGON', 'PENNSYLVANIA', 'SOUTH CAROLINA', 'SOUTH DAKOTA',
 'TENNESSEE', 'TEXAS', 'UTAH', 'VIRGINIA', 'WASHINGTON', 'WEST VIRGINIA', 'WISCONSIN', 'WYOMING'],

'Americas':['AC', 'AGUASCALIENTES', 'AL', 'ALABAMA', 'AM', 'AP', 'ARIZONA', 'ARKANSAS', 'BA', 'BAJA CALIFORNIA',
            'BELIZE', 'BOLIVIA', 'BUENOS AIRES', 'CABO VERDE', 'CALIFORNIA', 'CANADA', 'CATAMARCA', 'CE', 'CHACO',
 'CHIAPAS', 'CHIHUAHUA', 'CHILE', 'COLOMBIA', 'COLORADO', 'CORDOBA', 'CORRIENTES', 'DELAWARE', 'DF', 'DISTRITO FEDERAL',
 'DURANGO', 'ECUADOR', 'EL SALVADOR', 'ENTRE RIOS', 'ES', 'FLORIDA', 'FRENCH GUIANA', 'GEORGIA', 'GO', 'GUATEMALA',
 'IDAHO', 'ILLINOIS', 'INDIANA', 'IOWA', 'KANSAS', 'KENTUCKY', 'LA PAMPA', 'LOUISIANA', 'MA', 'MARYLAND', 'MEXICO',
 'MG', 'MICHIGAN', 'MINNESOTA', 'MISIONES', 'MISSISSIPPI', 'MISSOURI', 'MONTANA', 'MS', 'MT', 'NEBRASKA', 'NEW JERSEY',
 'NEW MEXICO', 'NEW YORK', 'NORTH CAROLINA', 'NORTH DAKOTA', 'OAXACA', 'OHIO', 'OKLAHOMA', 'OREGON', 'PA', 'PANAMA',
 'PARAGUAY', 'PB', 'PE', 'PENNSYLVANIA', 'PERU', 'PI', 'PR', 'RJ', 'RN', 'RO', 'RR', 'RS', 'SAN LUIS', 'SAN LUIS POTOSI',
 'SANTA FE', 'SANTIAGO DEL ESTERO', 'SC', 'SE', 'SINALOA', 'SONORA', 'SOUTH CAROLINA', 'SOUTH DAKOTA', 'SP', 'TABASCO',
 'TENNESSEE', 'TEXAS', 'TLAXCALA', 'TO', 'TUCUMAN', 'URUGUAY', 'UTAH', 'VENEZUELA', 'VERACRUZ', 'VIRGINIA', 'WASHINGTON', 
 'WEST VIRGINIA', 'WISCONSIN', 'WYOMING', 'YUCATAN', 'ZACATECAS'], 

'North China Plain': ['Liaoning','Shandong','Hebei', 'Shanxi','Anhui','Shaanxi','Henan','Jiangsu']
}

wRegions={
        
'United States': [ 'ALABAMA', 'ARIZONA', 'ARKANSAS', 'CALIFORNIA', 'COLORADO', 'DELAWARE', 'FLORIDA',
 'GEORGIA', 'IDAHO', 'ILLINOIS', 'INDIANA', 'IOWA', 'KANSAS', 'KENTUCKY', 'LOUISIANA', 'MARYLAND', 'MICHIGAN',
 'MINNESOTA', 'MISSISSIPPI', 'MISSOURI', 'MONTANA', 'NEBRASKA', 'NEW JERSEY', 'NEW MEXICO', 'NEW YORK', 'NORTH CAROLINA',
 'NORTH DAKOTA', 'OHIO', 'OKLAHOMA', 'OREGON', 'PENNSYLVANIA', 'SOUTH CAROLINA', 'SOUTH DAKOTA', 'TENNESSEE', 'TEXAS',
 'UTAH', 'VIRGINIA', 'WASHINGTON', 'WEST VIRGINIA', 'WISCONSIN', 'WYOMING'],
        
'US Great Plains': ['KANSAS','OKLAHOMA','TEXAS','NORTH DAKOTA','SOUTH DAKOTA','COLORADO',
       'NEBRASKA','WYOMING','MONTANA','NEW MEXICO'], 
                    
'West India': ['Rajasthan','Gujarat','Madhya Pradesh','Maharashtra','Andhra Pradesh',
        'Karnataka'], #'Odisha'?  

'India': ['Andaman and Nicobar', 'Andhra Pradesh', 'Arunachal Pradesh',
       'Assam', 'Bihar', 'Chandigarh', 'Chhattisgarh',
       'Dadra and Nagar Haveli', 'Daman and Diu', 'Goa', 'Gujarat',
       'Haryana', 'Himachal Pradesh', 'Jammu and Kashmir', 'Jharkhand',
       'Karnataka', 'Kerala', 'Lakshadweep', 'Madhya Pradesh',
       'Maharashtra', 'Manipur', 'Meghalaya', 'Mizoram', 'NCT of Delhi',
       'Nagaland', 'Odisha', 'Puducherry', 'Punjab', 'Rajasthan', 'Sikkim',
       'Tamil Nadu', 'Telangana', 'Tripura', 'Uttar Pradesh',
       'Uttarakhand', 'West Bengal'],
               
'Southeast Australia':['Victoria','New South Wales(b)','Queensland', 'South Australia'],

'Ethiopia':['Ethiopia'],

'East Africa':['Ethiopia','Kenya','Somalia','Somaliland','Eritrea','Djibouti','Uganda'],

'Southeast Africa': ['South Africa','Zimbabwe','Malawi','Mozambique','Zambia','Botswana'],#'Tanzania'

'Southeast South America':['Uruguay','Buenos Aires','Entre Rios','La Pampa','San Luis','Santiago del Estero',
        'Cordoba','Tucuman','Chaco','Corrientes','PR','SC','RS'],

'Americas':[
'ALABAMA', 'ALBERTA','ARIZONA',
 'ARKANSAS','BOLIVIA','BRITISH COLUMBIA',
 'BUENOS AIRES','CALIFORNIA','CHACO','CHILE','COLOMBIA','COLORADO',
 'DELAWARE','DF','ENTRE RIOS','FLORIDA','FORMOSA','GEORGIA','GO','IDAHO','ILLINOIS','INDIANA', 'IOWA',
 'JUJUY', 'KANSAS','LA PAMPA','MARYLAND', 'MEXICO', 'MG','MICHIGAN', 'MINNESOTA', 'MISSISSIPPI','MISSOURI',
 'MONTANA','MS', 'MT','NEBRASKA', 'NEVADA', 'NEW BRUNSWICK',
'NEW JERSEY', 'NEW MEXICO', 'NEW YORK','NORTH CAROLINA','NORTH DAKOTA', 'NOVA SCOTIA', 'OHIO',
 'OKLAHOMA', 'ONTARIO', 'OREGON', 'PARAGUAY', 'PENNSYLVANIA', 'PERU', 'PR',
 'QUEBEC', 'RS', 'SAN LUIS', 'SANTA FE', 'SC', 'SOUTH CAROLINA','SOUTH DAKOTA', 'SP', 'TENNESSEE',
 'TEXAS', 'URUGUAY', 'UTAH', 'VENEZUELA', 'VIRGINIA', 'WASHINGTON', 'WISCONSIN', 'WYOMING',
 'RO','AC','AM','AP','PA','TO','MA','PI','CE','RN','PB','PE',
 'AL','SE','BA','MT','MS','GO','DF','MG','ES','RJ','SP','PR','SC','RS'],                         
                         
'North Africa': ['Morocco','Algeria','Tunisia'],

'North Africa extended': ['Morocco','Algeria','Tunisia','Libya','Egypt'],

'Iberian Peninsula': ['Spain','Portugal'],

'IP + N. Af': ['Morocco','Algeria','Tunisia','Egypt','Spain','Portugal'],

'Europe+CSR-IP':['France','Germany','United Kingdom',
'Sweden','Denmark','Norway','Ireland','Poland','Austria','Greece','Italy',
'Yugoslavia','Czechoslovakia','Belgium-Luxembourg','Romania','Bulgaria','Ukraine',
'Serbia','Lithuania','Latvia','Estonia','USSR'],

             
'Europe':['Spain','Portugal','France','Germany','United Kingdom',
'Sweden','Denmark','Norway','Ireland','Poland','Austria','Greece','Italy',
'Yugoslavia','Czechoslovakia','Belgium-Luxembourg','Romania','Bulgaria','Ukraine',
'Serbia','Lithuania','Latvia','Estonia'],

'East Europe':['Romania','Yugoslavia','Czechoslovakia','USSR'],

               
'Europe+CSR':['Spain','Portugal','France','Germany','United Kingdom',
'Sweden','Denmark','Norway','Ireland','Poland','Austria','Greece','Italy',
'Yugoslavia','Czechoslovakia','Belgium-Luxembourg','Romania','Bulgaria','Ukraine',
'Serbia','Lithuania','Latvia','Estonia','USSR'],

'Eurasia-IND-CHN':[
'Spain','Portugal','France','Germany','United Kingdom','Afghanistan','Iran',
'Sweden','Denmark','Norway','Turkey','Romania','Bulgaria','Ukraine',
'Serbia','Lithuania','Latvia','Estonia','Syria','Pakistan','Iraq',
'Ireland','Poland','Austria','Belarus','Greece','Italy',
'Mongolia','Yugoslavia','Czechoslovakia','Belgium-Luxembourg','USSR'],

           
           
'Brazil-S.BR':['RO','AC','AM','AP','PA','TO','MA','PI','CE','RN','PB','PE',
        'AL','SE','BA','MT','MS','GO','DF','MG','ES','RJ','SP','PR','SC','RS']


}


sRegions = {   
'US Midwest': ['IOWA','ILLINOIS','MISSOURI','INDIANA','KENTUCKY',
       'OHIO','SOUTH DAKOTA','NEBRASKA','KANSAS','NORTH DAKOTA','WISCONSIN','MICHIGAN','MINNESOTA'],#

'United States': [
 'ALABAMA', 'ARIZONA', 'ARKANSAS', 'CALIFORNIA', 'COLORADO', 'DELAWARE', 'FLORIDA', 'GEORGIA', 'IDAHO', 'ILLINOIS',
 'INDIANA', 'IOWA', 'KANSAS', 'KENTUCKY', 'LOUISIANA', 'MARYLAND', 'MICHIGAN', 'MINNESOTA', 'MISSISSIPPI', 'MISSOURI',
 'MONTANA', 'NEBRASKA', 'NEW JERSEY', 'NEW MEXICO', 'NEW YORK', 'NORTH CAROLINA', 'NORTH DAKOTA', 'OHIO',
 'OKLAHOMA', 'OREGON', 'PENNSYLVANIA', 'SOUTH CAROLINA', 'SOUTH DAKOTA', 'TENNESSEE', 'TEXAS', 'UTAH', 'VIRGINIA', 'WASHINGTON',
 'WEST VIRGINIA', 'WISCONSIN', 'WYOMING'],
               
'West Africa':["Côte d'Ivoire",'Senegal','Gambia','Guinea','Guinea-Bissau','Mali',
'Sierra Leone','Ghana','Liberia','Togo','Benin','Burkina Faso','Nigeria','Cameroon'],

'Southeast Africa': ['South Africa','Zimbabwe','Malawi','Mozambique','Zambia','Botswana'],#'Tanzania'

'East Africa':['Ethiopia','Kenya','Somalia','Somaliland','Eritrea','Djibouti','Uganda'],

'Southeast South America':['Uruguay','Buenos Aires','Entre Rios','La Pampa','San Luis','Santiago del Estero',
        'Cordoba','Tucuman','Chaco','Corrientes','PR','SC','RS'],
                           
'Brazil':['RO','AC','AM','AP','PA','TO','MA','PI','CE','RN','PB','PE',
        'AL','SE','BA','MT','MS','GO','DF','MG','ES','RJ','SP','PR','SC','RS'],

                                     
'Brazil-S.BR':['RO','AC','AM','AP','PA','TO','MA','PI','CE','RN','PB','PE',
        'AL','SE','BA','MT','MS','GO','DF','MG','ES','RJ','SP','PR','SC','RS'],
          
'Brazil CW+SE':['MT','MS','GO','DF','MG','ES','RJ','SP'],
          
          
'Brazil_S':['PR','SC','RS'],
          
'North China Plain': ['Liaoning','Shandong','Hebei', 'Shanxi','Anhui','Shaanxi','Henan','Jiangsu']
}


naoRegions={'North Africa': ['Morocco','Algeria','Tunisia']}