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

# load required 
# CHANGE DATE
with open(f"../Input_files/jh_daily_ir_01_16_22.pkl",'rb') as file: #Contains only country name for each date
    jh_daily_ir = pickle.load(file)
# CHANGE DATE  
with open(f"../Input_files/jh_data_01_16_22.pkl",'rb') as file:
    jh_data=pickle.load(file)
# CHANGE DATE    
with open(f"../Input_files/jh_daily_deaths_01_16_22.pkl",'rb') as file:
    jh_daily_deaths = pickle.load(file)    
# CHANGE DATE
with open(f"../Input_files/seq_country_counts_01_16_22.pkl",'rb') as file:
    sequence_country_counts = pickle.load(file)
# CHANGE DATE
with open(f"../Input_files/basis_vector_IR_FR_PP_01_16_22.pkl",'rb') as file:
    lagged_FR1 = pickle.load(file)
# CHANGE DATE
with open(f"../Input_files/aggregated_mutations_window00_12_16_21_01_16_22.pkl",'rb') as file:
    aggregated_mutations_window = pickle.load(file)
    

date_range = jh_daily_ir.columns[9::]    

lag_tab = pd.read_csv("../Input_files/lag_tab_edited.csv")
ctr=lag_tab['country']


#Looping over single days
list_dates=jh_daily_ir.columns[9::]

# For this step, create a few parallel scripts, each cycling through ~15 days. For example, if your update covers 30 days, you would create 2 scripts for each IR, PP and lagged_FR.
# Get the indexes corresponding to days by looking up list_dates as before.

# CHANGE INDEXES
for date in tqdm(list_dates[684:700]):
#Looping over single variants
    for var in aggregated_mutations_window[date].keys():
        #For each variant count countries, 
        #totalCounts=len(aggregated_mutations_window[date][var][0])
        for counted in Counter(aggregated_mutations_window[date][var]).items(): #counted is country and counts
            date_lag=date
            if (ctr==counted[0]).any():
                lag=-(lag_tab.loc[lag_tab['country']==counted[0],"lag"].values[0])
                if (date+timedelta(days=lag.tolist()))<=date_range[-1]:
                    date_lag=date+timedelta(days=lag.tolist())
#            print(date)
#            print(date_lag)
            if sequence_country_counts[date][counted[0]]==0:
                lagged_FR1[date][var].append((counted[1])*jh_daily_deaths.loc[counted[0],date_lag])
            else:    
                lagged_FR1[date][var].append(round((counted[1]/sequence_country_counts[date][counted[0]])*jh_daily_deaths.loc[counted[0],date_lag],2))
      
    
# CHANGE DATE
with open(f"../Input_files/lagged_FR1_011622.pkl",'wb') as file:
    pickle.dump(lagged_FR1,file,protocol=pickle.HIGHEST_PROTOCOL)
