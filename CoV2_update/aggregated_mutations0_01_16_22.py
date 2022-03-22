import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from datetime import datetime as dt
from datetime import timedelta
from scipy.signal import lfilter
import pickle as pickle
import sys

# load required files
# CHANGE DATE
with open(f"../Input_files/jh_daily_ir_01_16_22.pkl",'rb') as file: #Contains only country name for each date
    jh_daily_ir = pickle.load(file)
    
# CHANGE DATE
with open(f"../Input_files/jh_data_01_16_22.pkl",'rb') as file:
    jh_data=pickle.load(file)
    
# CHANGE DATE    
with open(f"../Input_files/meta_01_16_22.pkl",'rb') as file:
    meta = pickle.load(file)

# load previous files to diff with

# CHANGE BOTH FILES BELOW TO LATEST DIFF
#in this case, 1-16-22. The below refers to latest update run (1-16-22, and the previous diff was 12-15-21

# date = '12_15_21'
with open(f"../Input_files/jh_cncb_daily_dates_12_15_21.pkl",'rb') as file: #Contains only country name for each date
    aggregated_mutations = pickle.load(file)
    
with open(f"../Input_files/jh_cncb_daily_dates_identifiers_12_15_21.pkl",'rb') as file:
    processed_identifiers = pickle.load(file)
    
# gff3 parsing


# CHANGE TO THE SINK FOLDER YOU WILL USE FOR UPDATE (create the folder in /gpfs/group/balch/data)    
all_names = os.listdir('/gpfs/group/balch/data/gff3_cncb_01_16_22')
column_names = ['variant type','start','end', 'info']
missing = []
date_range = pd.date_range(start=jh_daily_ir.columns[0], end=jh_daily_ir.columns[-1])

#Check if there is a processed_identifier list, contains genome names already processed
if 'processed_identifiers' not in locals():
    processed_identifiers = []
    
#aggregated mutations, dictionary for each date, each sub dictionary contains a 3d list structure. [[],[],[]]
if 'aggregated_mutations' not in locals():
    #this is for country ir and fr for that specific day
    aggregated_mutations = {date:{} for date in date_range}
else:
    for date in date_range:
        if date not in aggregated_mutations.keys():
            aggregated_mutations[date] = {}

#for identifier in tqdm(meta[meta['Country']=='United States'].index):
for identifier in tqdm(meta.index):
    
    #Use prexisting list
    if identifier in processed_identifiers:
        continue

    #Searching for correct identifier
    #--------------------------
    #No alternate name is ' '
    file_name = ''
    #Check if accession id in file names, if not check related ids
    if '2019-nCoV_'+identifier+'_variants.gff3' in all_names:
        file_name = '2019-nCoV_'+identifier+'_variants.gff3'
    # checking alternate names
    elif meta.loc[identifier,'Related ID'] != ' ':
        for alt_identifier in meta.loc[identifier,'Related ID'].replace(' ','').split(','):
            if '2019-nCoV_'+alt_identifier+'_variants.gff3' in all_names:
                file_name = '2019-nCoV_'+alt_identifier+'_variants.gff3'
                break
        #Added in case alternate names are also not found in gffs
        if file_name == '':
            missing.append(identifier)
            continue
    # If file name has not been updated, then there is no matching identifier, move to next index
    elif file_name == '':
        missing.append(identifier)
        continue
    #--------------------------
    
    #Filtering files with no variants
    #--------------------------
    # CHANGE TO THE SINK FOLDER YOU WILL USE FOR UPDATE (create the folder in /gpfs/group/balch/data) 
    with open(f'/gpfs/group/balch/data/gff3_cncb_01_16_22/{file_name}') as text_file:
        lines = text_file.readlines()
        counter = 0
        for l in lines:
            if '#' in l:
                counter += 1
    #Number of info lines should be less than total, if not then there are no mutations
    #--------------------------
    
    #List to keep track of which files used already, save and import this in future to avoid redundant search
    processed_identifiers.append(identifier)
    
    if counter<len(lines) and meta.loc[identifier,'Country'] in set(jh_data['Country_Region']): #country needs to be found in jh data for inf/fata rates
        
        # CHANGE TO THE SINK FOLDER YOU WILL USE FOR UPDATE (create the folder in /gpfs/group/balch/data) 
        gff = pd.read_csv(f'/gpfs/group/balch/data/gff3_cncb_01_16_22/{file_name}',sep='\t',skiprows=counter,usecols=[1,3,4,8],names=column_names)
        info_df = pd.DataFrame(gff['info'].str.split(';').values.tolist(),columns=[0,1,'Ref','Alt','Description']).drop([0,1],axis=1)
        gff = gff.drop(['info'],axis=1)
        gff['Country'] = [meta.loc[identifier,'Country']]*gff.shape[0]
        temp_df = pd.concat([gff,info_df],axis=1)
        
        #Filtering alternate amino acid and reference for missense_variant and synonymous_variant
        missenses_ref = temp_df.loc[temp_df['Description'].str.contains('missense_variant'),'Description'].str.split(',').str[1].str[-3]
        synonymous_ref = temp_df.loc[temp_df['Description'].str.contains('synonymous_variant'),'Description'].str.split(',').str[1].str[-1]
        
        # exception handling for sporadic cases with both synonymous and missense annotation
        idx=missenses_ref.index.intersection(synonymous_ref.index)
        if not idx.empty:
            temp_df=temp_df.drop(idx)
            missenses_ref = temp_df.loc[temp_df['Description'].str.contains('missense_variant'),'Description'].str.split(',').str[1].str[-3]
            synonymous_ref = temp_df.loc[temp_df['Description'].str.contains('synonymous_variant'),'Description'].str.split(',').str[1].str[-1]
            print(file_name,idx)
        
        temp_df['Ref_AA'] = pd.concat([missenses_ref,synonymous_ref])
        temp_df['Alt_AA'] = temp_df.loc[temp_df['Description'].str.contains('missense_variant'),'Description'].str.split(',').str[1].str[-1]

        missenses_str = temp_df.loc[temp_df['Description'].str.contains('missense_variant'),'Description'].str.split(',').str[1].str.split('.').str[-1]
        synonymous_str = temp_df.loc[temp_df['Description'].str.contains('synonymous_variant'),'Description'].str.split(',').str[1].str.split('.').str[-1]
        temp_df['AA'] = pd.concat([missenses_str,synonymous_str])
        
        temp_df.fillna('', inplace=True)
        temp_df['descriptor'] = temp_df['start'].astype(str)+','+temp_df['end'].astype(str)+','+temp_df['Ref']+','+temp_df['Alt']+\
        ','+temp_df['Description'].str.split(',').str[0].str.split('=').str[1]+','+temp_df['variant type']+','+temp_df['Ref_AA']+','+temp_df['Alt_AA']+\
        ','+temp_df['AA']

#         if any(temp_df['Ref'].str.contains('REF=A') & temp_df['Alt'].str.contains('ALT=GGTTC')):
#             break
        
    
        for var in temp_df['descriptor']:
            #Dictionary for each sample date, each df then has its own variant descriptor (duplicates across dates)
            if var not in aggregated_mutations[meta.loc[identifier,'Sample Collection Date']].keys():
                aggregated_mutations[meta.loc[identifier,'Sample Collection Date']][var] = []

            #Use country and sample collection date to parse jh_daily_ir/jh_daily_fr
            #Country
            aggregated_mutations[meta.loc[identifier,'Sample Collection Date']][var].append(temp_df.loc[0,'Country'])

            
            
print(len(missing))
# CHANGE DATE
with open(f"../Input_files/jh_cncb_daily_dates_01_16_22.pkl",'wb') as file:
    pickle.dump(aggregated_mutations,file,protocol=pickle.HIGHEST_PROTOCOL)
# CHANGE DATE
with open(f"../Input_files/jh_cncb_daily_dates_identifiers_01_16_22.pkl",'wb') as file:
    pickle.dump(processed_identifiers,file,protocol=pickle.HIGHEST_PROTOCOL)
# CHANGE DATE    
with open(f"../Input_files/missing_01_16_22.pkl",'wb') as file:
    pickle.dump(missing,file,protocol=pickle.HIGHEST_PROTOCOL)


