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
    
# load aggregated_mutations (run separately through aggregated_mutations0 script)

# CHANGE DATE - output of aggregated_mutations0 script
with open(f"../Input_files/jh_cncb_daily_dates_01_16_22.pkl",'rb') as file: #Contains only country name for each date
    aggregated_mutations = pickle.load(file)
    
# load the above file
# CHANGE DATE - output of sequence country counts script
with open(f"../Input_files/seq_country_counts_01_16_22.pkl",'rb') as file:
    sequence_country_counts = pickle.load(file)
    
    
date_range = jh_daily_ir.columns[9::] #skip the first 9 (feb 1), arbitrarily since not too much data yet
aggregated_mutations_window = {date:{} for date in date_range}

#Gets cumulative country counts up to date

#iterate over dates to get cumlative up to that date
for date in tqdm(date_range):
    #iterate over the preceding dates for the day that you are calculating
    for window_date in pd.date_range(start=date_range[0],end=date):
        for var in aggregated_mutations[window_date].keys():
            if var not in aggregated_mutations_window[date].keys():
                aggregated_mutations_window[date][var] = []
                
            #Use country and date to get same ir/fr for that period
            aggregated_mutations_window[date][var] += aggregated_mutations[window_date][var] #Appending country names

# UPDATE INDEXES IN THE DATE_RANGE BELOW - BEGINNING AND END DATE OF YOUR UPDATE
# To get the integers corresponding to such dates, open a console, load jh_daily_ir, obtain date_range (as per code above) and then look it up. For example, in the below '684' is the index for 12-16-21 (beginning of the update) and '717' is the index for 1-17-21 (end of the update).
aggregated_mutations_window_update={k:aggregated_mutations_window[k] for k in date_range[684:717]}  #12_16_21 - 01_17_22

# CHANGE DATES ACCORDING TO BEGINNING / END OF THE UPDATE
with open(f"../Input_files/aggregated_mutations_window00_12_16_21_01_16_22.pkl",'wb') as file:
    pickle.dump(aggregated_mutations_window_update,file,protocol=pickle.HIGHEST_PROTOCOL)
    

    
# get a vector for the entire time interval, required for computing the cumulative values in the next steps (since the above contains only the update portion)
date_range = jh_daily_ir.columns[9::]
basis = {date:{} for date in date_range}

#Gets cumlative country counts up to date

#iterate over dates to get cumlative up to that date
for date in tqdm(date_range):
    #iterate over the preceding dates for the day that you are calculating
    for window_date in pd.date_range(start=date_range[0],end=date):
        for var in aggregated_mutations[window_date].keys():
            basis[date][var] = []  # slot added
# CHANGE DATE            
with open(f"../Input_files/basis_vector_IR_FR_PP_01_16_22.pkl",'wb') as file:
    pickle.dump(basis,file,protocol=pickle.HIGHEST_PROTOCOL)