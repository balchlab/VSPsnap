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
from itertools import combinations

# import day number
#input1 = sys.argv[1]
#input1 = int(input1)

with open(f"../Input_files/jh_daily_ir_03_20_22.pkl",'rb') as file: #Contains only country name for each date
    jh_daily_ir = pickle.load(file)
    
with open(f"../Input_files/meta_03_20_22.pkl",'rb') as file:
    meta = pickle.load(file)


with open(f"../Input_files/unique_muts_filt_VO.pkl",'rb') as file:
    unique_muts_filt = pickle.load(file)

    

co_occurrence=pd.DataFrame(index=unique_muts_filt.index,columns=unique_muts_filt.index)
co_occurrence.fillna(0,inplace=True)


#sparse
s_co = co_occurrence.astype(pd.SparseDtype("int"))
del co_occurrence


def sp_loc(df, index, columns, val):
    """ Insert data in a DataFrame with SparseDtype format

    Only applicable for pandas version > 0.25

    Args
    ----
    df : DataFrame with series formatted with pd.SparseDtype
    index: str, or list, or slice object
        Same as one would use as first argument of .loc[]
    columns: str, list, or slice
        Same one would normally use as second argument of .loc[]
    val: insert values

    Returns
    -------
    df: DataFrame
        Modified DataFrame

    """

    # Save the original sparse format for reuse later
    spdtypes = df.dtypes[columns]

    # Convert concerned Series to dense format
    df[columns] = df[columns].sparse.to_dense()

    # Do a normal insertion with .loc[]
    df.loc[index, columns] = val

    # Back to the original sparse format
    df[columns] = df[columns].astype(spdtypes)

    return df

date_range = jh_daily_ir.columns[755:765]#788]

all_names = os.listdir('/gpfs/group/balch/data/gff3_cncb_03_20_22')
column_names = ['variant type','start','end', 'info']
missing = []

for date in tqdm(date_range):
    meta_daily=meta[(meta['Sample Collection Date']==date)]
    
    for identifier in tqdm(meta_daily.index):
    #Searching for correct identifier
    #--------------------------
        #No alternate name is ' '
        file_name = ''
        #Check if accession id in file names, if not check related ids
        if '2019-nCoV_'+identifier+'_variants.gff3' in all_names:
            file_name = '2019-nCoV_'+identifier+'_variants.gff3'
        # checking alternate names
        elif meta_daily.loc[identifier,'Related ID'] != ' ':
            for alt_identifier in meta_daily.loc[identifier,'Related ID'].replace(' ','').split(','):
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
        with open(f'/gpfs/group/balch/data/gff3_cncb_03_20_22/{file_name}') as text_file:
            lines = text_file.readlines()
            counter = 0
            for l in lines:
                if '#' in l:
                    counter += 1
        #Number of info lines should be less than total, if not then there are no mutations
        #--------------------------

        #List to keep track of which files used already, save and import this in future to avoid redundant search
        #processed_identifiers.append(identifier)

        if counter<len(lines): #and meta.loc[identifier,'Country'] in set(jh_data['Country_Region']): #country needs to be found in jh data for inf/fata rates

            gff = pd.read_csv(f'/gpfs/group/balch/data/gff3_cncb_03_20_22/{file_name}',sep='\t',skiprows=counter,usecols=[1,3,4,8],names=column_names)
            info_df = pd.DataFrame(gff['info'].str.split(';').values.tolist(),columns=[0,1,'Ref','Alt','Description']).drop([0,1],axis=1)
            gff = gff.drop(['info'],axis=1)
            gff['Country'] = [meta.loc[identifier,'Country']]*gff.shape[0]
            temp_df = pd.concat([gff,info_df],axis=1)

            #Filtering alternate amino acid and reference for missense_variant and synonymous_variant
            missenses_ref = temp_df.loc[temp_df['Description'].str.contains('missense_variant'),'Description'].str.split(',').str[1].str[-3]
            synonymous_ref = temp_df.loc[temp_df['Description'].str.contains('synonymous_variant'),'Description'].str.split(',').str[1].str[-1]
            temp_df['Ref_AA'] = pd.concat([missenses_ref,synonymous_ref])
            temp_df['Alt_AA'] = temp_df.loc[temp_df['Description'].str.contains('missense_variant'),'Description'].str.split(',').str[1].str[-1]

            missenses_str = temp_df.loc[temp_df['Description'].str.contains('missense_variant'),'Description'].str.split(',').str[1].str.split('.').str[-1]
            synonymous_str = temp_df.loc[temp_df['Description'].str.contains('synonymous_variant'),'Description'].str.split(',').str[1].str.split('.').str[-1]
            temp_df['AA'] = pd.concat([missenses_str,synonymous_str])

            temp_df.fillna('', inplace=True)
            temp_df['descriptor'] = temp_df['start'].astype(str)+','+temp_df['end'].astype(str)+','+temp_df['Ref']+','+temp_df['Alt']+\
            ','+temp_df['Description'].str.split(',').str[0].str.split('=').str[1]+','+temp_df['variant type']+','+temp_df['Ref_AA']+','+temp_df['Alt_AA']+\
            ','+temp_df['AA']
            
            if(sum(temp_df['descriptor'].isin(s_co.index))>1): # overlap with mutations in co-occurrence matrix at least 2
                temp_df_filt=temp_df[temp_df['descriptor'].isin(s_co.index)]
            
                pairwise=list(combinations(temp_df_filt['descriptor'],2))   # pairwise combos
                for item in pairwise:
                    val=s_co.loc[item[0],item[1]]+1
                    s_co = sp_loc(s_co,item[0],item[1],val)
                    #co_occurrence.loc[item[0],item[1]]+=1
                    #s_co = co_occurrence.astype(pd.SparseDtype("int"))
                    
    s_co.to_pickle(f"../Output Files/co_occurrence/daily/update_03_20_22_VO/cc_daily_{date.strftime('%m_%d_%y')}.pkl",protocol=pickle.HIGHEST_PROTOCOL)

# a=range(570, 684)
# print(*a)
      
# running:

# test:
# ./foo_spawner_slurm.sh run_co_singleDay_12_15_21.sh 570

# run:
# ./foo_spawner_slurm.sh run_co_singleDay_12_15_21.sh 571 572 573 574 575 576 577 578 579 580 581 582 583 584 585 586 587 588 589 590 591 592 593 594 595 596 597 598 599 600 601 602 603 604 605 606 607 608 609 610 611 612 613 614 615 616 617 618 619 620 621 622 623 624 625 626 627 628 629 630 631 632 633 634 635 636 637 638 639 640 641 642 643 644 645 646 647 648 649 650 651 652 653 654 655 656 657 658 659 660 661 662 663 664 665 666 667 668 669 670 671 672 673 674 675 676 677 678 679 680 681 682 683 684



