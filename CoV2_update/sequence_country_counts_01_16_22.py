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
    
date_range = jh_daily_ir.columns[9::] #skip the first 9 (feb 1), arbitrarily since not too much data yet
sequence_country_counts = {date:dict.fromkeys(set(meta['Country'].str.lower()),0) for date in date_range} #Dictionary of dates, sub dictionary country names with counts

for date in tqdm(date_range):
    for country in set(meta['Country'].str.lower()):
        sequence_country_counts[date][country] = meta[(meta['Sample Collection Date']<=date) & (meta['Country']==country)].shape[0]
        
        
# CHANGE DATE     
with open(f"../Input_files/seq_country_counts_01_16_22.pkl",'wb') as file:
    pickle.dump(sequence_country_counts,file,protocol=pickle.HIGHEST_PROTOCOL)
    
