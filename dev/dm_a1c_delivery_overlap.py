#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2021-03-11 07:23:35


import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

%matplotlib inline
%config InlineBackend.figure_format='retina'


DATE = datetime.now().strftime('%Y-%m-%d')

ALL_EGA_FILE="/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/combined_cohorts_07_12_2018/full_dataset_characterization/expanded_ega/date_of_conception_w_ega.tsv"
# EGA_FILE = "/dors/capra_lab/users/abraha1/prelim_studies/crp_gest_age_assoc/data/EGA_w-in_3days_of_delivery.tsv"
EGA_FILE="/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/combined_cohorts_07_12_2018/full_dataset_characterization/updated_ega_2020_06_16/earliest_delivery_with_date_of_conception_w_ega_updated_EGA_2020_06_17.tsv"
DELIVERY_FILE = "/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/combined_cohorts_07_12_2018/full_dataset_characterization/est_delivery_date_at_least_one_icd_cpt_ega.tsv"
CLIN_FILE = "/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/filtered_labs/filtered_top100_by_num_grids_and_min_4_vals_per_grid.tsv"
DM_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_03_11_from_annika_t1_t2_dates")
t1_date_file = DM_DIR.joinpath("earliest_date_T1D_ICD.txt")
t2_date_file = DM_DIR.joinpath("earliest_date_T2D_ICD.txt")


OUTPUT_DIR = Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_03_dm_a1c_delivery_counts")
LABS_DIR = Path("/dors/capra_lab/users/abraha1/data/davis_clinical_labs/selected_lab_summaries")
A1C_FILES = [LABS_DIR.joinpath('HgbA1C', "20200108_HgbA1C_4sd_outliers_removed.txt"), LABS_DIR.joinpath('HGB_A1C', "20200108_HGB_A1C_4sd_outliers_removed.txt")]

# relevant labs
# C-Gluc/
# F-Gluc/
# Glucose_HH/
# GLUCOSE LEVEL/
# GLUCOSE_PD_FLUID_(TIMED)/
# GLU-EC/
# GlutAc/
# Glutmn/
# GluWB/
# GT2Hr/
# GTFst/
# GTT1Hr/
# GTT2Hr/
# GTT3Hr/
# GTT50g/
# GTTFst/

# HgbA1C/
# HGB_A1C/




# -----------
# function
# -----------

def t1dm_date_to_delivery(x):

    t1_ = x['earliest_t1_date']
    delivery_ = x['consensus_delivery']
    t1_diff_delivery = t1_ - delivery_

    if (t1_diff_delivery > -1*np.timedelta64(9, 'M')) & (t1_diff_delivery < np.timedelta64(0, 'M')):
        return 'within9mo_delivery'
    elif (t1_diff_delivery < -1*np.timedelta64(9, 'M')):
        return 'before_psuedo_pregnancy'
    elif t1_diff_delivery > np.timedelta64(0, 'M'):
        return 't1_after_delivery'
    else:
        return np.nan

def t2dm_date_to_delivery(x):

    t2_ = x['earliest_t2_date']
    delivery_ = x['consensus_delivery']
    t2_diff_delivery = t2_ - delivery_

    if (t2_diff_delivery > -1*np.timedelta64(9, 'M')) & (t2_diff_delivery < np.timedelta64(0, 'M')):
        return 'within9mo_delivery'
    elif (t2_diff_delivery < -1*np.timedelta64(9, 'M')):
        return 'before_psuedo_pregnancy'
    elif t2_diff_delivery > np.timedelta64(0, 'M'):
        return 't2_after_delivery'
    else:
        return np.nan

def merge_a1c_t1_t2(a1c_df, t1_preg_df=None, t2_preg_df=None, no_dm_preg_df=None ):
    # take a t1, t2, or no dm dataframe and tkae the mean a1c value during pregnancy
    if isinstance(t1_preg_df, pd.DataFrame):
        dm_preg_df = t1_preg_df.copy()
        dm_date_name = 'earliest_t1_date'
    elif isinstance(t2_preg_df, pd.DataFrame):
        dm_preg_df = t2_preg_df.copy()
        dm_date_name = 'earliest_t2_date'
    elif isinstance(no_dm_preg_df, pd.DataFrame):
        dm_preg_df = no_dm_preg_df.copy()
        dm_date_name = False

    # merge a1 with specific cohort
    acols=['PERSON_ID','LAB_DATE','LAB_VALUE']
    dm_a1c_preg_df = a1c_df.loc[a1c_df['PERSON_ID'].isin(dm_preg_df['GRID']), acols].copy()


    # conver to date col
    if dm_date_name:
        datecols = ['LAB_DATE', 'delivery_date','date_of_conception',dm_date_name]
    else:
        datecols = ['LAB_DATE', 'delivery_date','date_of_conception']
    dm_preg_a1c_merged_df = pd.merge(dm_a1c_preg_df, dm_preg_df, left_on="PERSON_ID", right_on="GRID", how='inner')
    for dc in datecols:
        dm_preg_a1c_merged_df[dc] = pd.to_datetime(dm_preg_a1c_merged_df[dc])

    # when was the a1c measures?
    within_preg_bool=(dm_preg_a1c_merged_df['LAB_DATE'] >= dm_preg_a1c_merged_df['date_of_conception']) & (dm_preg_a1c_merged_df['LAB_DATE'] <= dm_preg_a1c_merged_df['delivery_date'])
    dm_preg_a1c_merged_df['context'] = np.nan
    dm_preg_a1c_merged_df.loc[within_preg_bool, 'context'] = 'in_preg'
    dm_preg_a1c_merged_df.loc[~within_preg_bool, 'context'] = 'outside_preg'

    # check if dx happened after diabetes
    if dm_date_name:
        dm_preg_a1c_merged_df['delivery_after_dm'] = dm_preg_a1c_merged_df['delivery_date'] > dm_preg_a1c_merged_df[dm_date_name]
    else:
        dm_preg_a1c_merged_df['delivery_after_dm'] = True

    in_preg_df = dm_preg_a1c_merged_df.loc[(dm_preg_a1c_merged_df['context']=='in_preg') & (dm_preg_a1c_merged_df['delivery_after_dm']==True)].copy()
    a1c_mean_df = in_preg_df.groupby('PERSON_ID').agg({'LAB_VALUE':'mean', 'consensus_label':[np.unique]})
    a1c_mean_df.columns = a1c_mean_df.columns.droplevel()

    return a1c_mean_df


# %%
# -----------
# MAIN
# -----------

#load a1c
a1c_df = pd.DataFrame()
for a1c_file in A1C_FILES:
    this_df = pd.read_csv(a1c_file, sep="\t")
    a1c_df = a1c_df.append(this_df)

a1c_df['PERSON_ID'].nunique()

# %%

# %%
# load deliveries
delivery_df = pd.read_csv(DELIVERY_FILE, sep="\t")
delivery_df = delivery_df.loc[:, ['GRID','consensus_delivery','consensus_label', 'delivery_id']].copy()
delivery_df['consensus_delivery'] = pd.to_datetime(delivery_df['consensus_delivery'])
first_delivery_df = delivery_df[~delivery_df.duplicated(['GRID'],keep='first')].copy()
first_delivery_df.shape


t1_df = pd.read_csv( t1_date_file, sep=" ")
t2_df = pd.read_csv( t2_date_file, sep=" ")
t1_df['min_targ_ICD_date'] = pd.to_datetime(t1_df['min_targ_ICD_date'])
t2_df['min_targ_ICD_date'] = pd.to_datetime(t2_df['min_targ_ICD_date'])

# are there duplicated GRIDS?
t1_df.duplicated('GRID', keep=False).sum()
t2_df.duplicated('GRID', keep=False).sum()
t1_df.columns=['GRID','earliest_t1_date']
t2_df.columns=['GRID','earliest_t2_date']

t1_df.shape
t2_df.shape

t2_df['GRID'].nunique()
# t1 and t2 overlap
t1_and_t2_grids = set(t2_df['GRID']).intersection(set(t1_df['GRID']))
len(t1_and_t2_grids)


# %%
# ega
all_ega_df = pd.read_csv(ALL_EGA_FILE, sep="\t")
all_ega_df.shape
all_ega_df['GRID'].nunique()
all_ega_df['delivery_date'] = pd.to_datetime(all_ega_df['delivery_date'])
all_ega_df.sort_values(['GRID','delivery_date'], inplace=True, ascending=True)
all_ega_df = all_ega_df[~all_ega_df['GRID'].duplicated(keep='first')].copy()
all_ega_df.shape
all_ega_df['GRID'].nunique()


ega_df = pd.read_csv(EGA_FILE, sep="\t")
ega_df = ega_df.loc[ (ega_df['more_than_315d']== False) &  (ega_df['less_than_0d']== False)].copy()
ega_df['delivery_date'] = pd.to_datetime(ega_df['delivery_date'] )
# 27,861 "GRIDs" w/ date of conception
ega_df['GRID'].nunique()


# try merging but not requiring exact delivery date match
test_df = pd.merge(all_ega_df, first_delivery_df, on=['GRID'], how='inner', suffixes=('_ega','_delivery'))
test_df['closest_ega_DATE'] = pd.to_datetime(test_df['closest_ega_DATE'])
test_df['delivery_date'] = pd.to_datetime(test_df['delivery_date'])

# within 3 days
test_df['diff'] = test_df['closest_ega_DATE'] - test_df['delivery_date']
within_3d_df = test_df.loc[ (test_df['diff'] < np.timedelta64(3, 'D')) & (test_df['diff'] >= np.timedelta64(0, 'D')) , ['GRID','delivery_date','consensus_label', 'closest_ega_DATE','date_of_conception']]
within_3d_df["GRID"].nunique()
within_3d_df.rename(columns={'delivery_date_delivery':'delivery_date' }, inplace=True)

### merge with t1 and t2
temp_df = pd.merge(within_3d_df, t1_df, on=['GRID'], how='outer')
concep_dm_df = pd.merge(temp_df, t2_df, on=['GRID'], how='outer')

# mark those with both t1 and t2 dx
concep_dm_df['t1_and_t2'] = concep_dm_df['GRID'].isin(t1_and_t2_grids)
concep_dm_df['GRID'].isin(t1_and_t2_grids).sum()

# 98 women w/ delivery type and conception + t1 diagnoses
# 115 women if you don't require an exact match on dlivery date but within 3 days
concep_dm_df[~concep_dm_df['delivery_date'].isna() & ~concep_dm_df['earliest_t1_date'].isna() & (concep_dm_df['t1_and_t2']==False)]
# 652 women w/ delivery type and conception + t1 diagnoses
# 761 women if you don't require an exact match on delivery date
concep_dm_df[~concep_dm_df['delivery_date'].isna() & ~concep_dm_df['earliest_t2_date'].isna() & (concep_dm_df['t1_and_t2']==False)]

# %%
kcols=['GRID','delivery_date','consensus_label','date_of_conception','earliest_t1_date']
t1_preg_df = concep_dm_df.loc[~concep_dm_df['delivery_date'].isna() & ~concep_dm_df['earliest_t1_date'].isna() & (concep_dm_df['t1_and_t2']==False), kcols].copy()
t1_preg_df.shape

kcols=['GRID','delivery_date','consensus_label','date_of_conception','earliest_t2_date']
t2_preg_df = concep_dm_df.loc[~concep_dm_df['delivery_date'].isna() & ~concep_dm_df['earliest_t2_date'].isna() & (concep_dm_df['t1_and_t2']==False), kcols].copy()
t2_preg_df.shape

kcols=['GRID','delivery_date','consensus_label','date_of_conception']
no_t1_t2_preg_df = concep_dm_df.loc[~concep_dm_df['delivery_date'].isna() & concep_dm_df['earliest_t2_date'].isna() & concep_dm_df['earliest_t2_date'].isna() & (concep_dm_df['t1_and_t2']==False), kcols].copy()
no_t1_t2_preg_df.shape

t1_a1c_mean_df = merge_a1c_t1_t2(a1c_df, t1_preg_df=t1_preg_df, )
t2_a1c_mean_df = merge_a1c_t1_t2(a1c_df, t2_preg_df=t2_preg_df  )
nodm_a1c_mean_df = merge_a1c_t1_t2(a1c_df, no_dm_preg_df=no_t1_t2_preg_df  )

term_preterm_only_t2_a1c_mean_df = t2_a1c_mean_df[t2_a1c_mean_df['unique'].isin(['preterm', 'term'])].copy()
term_preterm_only_nodm_a1c_mean_df = nodm_a1c_mean_df[nodm_a1c_mean_df['unique'].isin(['preterm', 'term'])].copy()

pd.value_counts(t2_a1c_mean_df['unique'])
pd.value_counts(nodm_a1c_mean_df['unique'])



# %%
fig,axs=plt.subplots(ncols=3, sharey=True, figsize=(8,4))

sns.stripplot(y=t1_a1c_mean_df["mean"],x=t1_a1c_mean_df["unique"],ax=axs[0], color='black')
sns.boxplot(y=t1_a1c_mean_df["mean"],x=t1_a1c_mean_df["unique"], ax=axs[0])
axs[0].set_title('T1D')
sns.stripplot(y=term_preterm_only_t2_a1c_mean_df["mean"],x=term_preterm_only_t2_a1c_mean_df["unique"],ax=axs[1], color='black')
sns.boxplot(y=term_preterm_only_t2_a1c_mean_df["mean"],x=term_preterm_only_t2_a1c_mean_df["unique"], ax=axs[1])
axs[1].set_title('T2D')

sns.stripplot(y=term_preterm_only_nodm_a1c_mean_df["mean"],x=term_preterm_only_nodm_a1c_mean_df["unique"],ax=axs[2], color='black')
sns.boxplot(y=term_preterm_only_nodm_a1c_mean_df["mean"],x=term_preterm_only_nodm_a1c_mean_df["unique"], ax=axs[2])
axs[2].set_title('No DM')



