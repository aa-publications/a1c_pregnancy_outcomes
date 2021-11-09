#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 'now'


import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path



import matplotlib.pyplot as plt
import seaborn as sns


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

%matplotlib inline
# %config InlineBackend.figure_format='retina'


DATE = datetime.now().strftime('%Y-%m-%d')


# %%
# -----------
# PATHS
# -----------

# get ICD codes for all of EHR
BILLING_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_09_28_sd_billing_data_from_adi")
icd9_file = BILLING_DIR.joinpath("data.V01/tidy_long/long_UTI.ICD9.V01.Feb2021.tsv.gz")
icd10_file = BILLING_DIR.joinpath("data.V01/tidy_long/long_UTI.ICD10.V01.Feb2021.tsv.gz")

query_file=BILLING_DIR.joinpath('query_icd9_10_codes.txt')


DM_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_03_11_from_annika_t1_t2_dates")
t1_date_file = DM_DIR.joinpath("earliest_date_T1D_ICD.txt")
t2_date_file = DM_DIR.joinpath("earliest_date_T2D_ICD.txt")

LABS_DIR = Path("/dors/capra_lab/users/abraha1/data/davis_clinical_labs/selected_lab_summaries")
A1C_FILES = [LABS_DIR.joinpath('HgbA1C', "20200108_HgbA1C_4sd_outliers_removed.txt"), LABS_DIR.joinpath('HGB_A1C', "20200108_HGB_A1C_4sd_outliers_removed.txt")]

PHE_DIR = Path("/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/phecodes/2021_09_27_phecodev2")
PREG_PHE_FILE = Path("/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/phecodes/2021_09_27_phecodev2/phecode_strings_V2_preg_chapter.csv")


# %%
# -----------
# LOAD DATA
# -----------

# phe
preg_phe_df = pd.read_csv( PREG_PHE_FILE, sep=",", names = ["phecode","phecode_string","phecode_category","sex","ICD10_only","phecode_top","leaf"])
phe_to_label_dict = dict(zip(preg_phe_df['phecode'], preg_phe_df['phecode_string']))
icd9_to_phe_df = pd.read_csv(  PHE_DIR.joinpath("ICD9_to_phecode_V2.csv"), sep=",")
icd10_to_phe_df = pd.read_csv( PHE_DIR.joinpath("ICD10_to_phecode_V2.csv"), sep=",")

# set up df w/ billing codes
icd9_preg_df = icd9_to_phe_df.loc[icd9_to_phe_df['phecode'].isin(preg_phe_df['phecode'])].copy()
icd10_preg_df = icd10_to_phe_df.loc[icd10_to_phe_df['phecode'].isin(preg_phe_df['phecode'])].copy()
icd9_to_phe_dict = dict(zip(icd9_preg_df['icd9'], icd9_preg_df['phecode']))
icd10_to_phe_dict = dict(zip(icd10_preg_df['icd10'], icd10_preg_df['phecode']))

# load icd9 and icd10 codes
i9_df = pd.read_csv( icd9_file, sep="\t", names=['GRID','DATE','CODE'])
i10_df = pd.read_csv( icd10_file, sep="\t", names=['GRID','DATE','CODE'])
i9_df['GRID'].nunique()
i10_df['GRID'].nunique()

# load supervision and high-risk pregnancy codes
query_df = pd.read_csv( query_file, sep="\t", names=['code'])
supervision_codes = [x.split('-')[0] for x in query_df['code'].values.tolist()]

# select individuals with supervision or high-risk pregnancy codes
selected_i9_df = i9_df[i9_df['CODE'].isin(supervision_codes)].copy()
selected_i10_df = i10_df[i10_df['CODE'].isin(supervision_codes)].copy()
selected_i9_df['DATE'] = pd.to_datetime(selected_i9_df['DATE'])

# pick the earliest supervision code in each chart
i9_10_df = pd.concat([selected_i9_df, selected_i10_df])
i9_10_df.sort_values(['GRID','DATE'], inplace=True)
cohort_df = i9_10_df[~i9_10_df.duplicated(['GRID'], keep='first')].copy()
cohort_df['GRID'].nunique()
cohort_df['DATE'] = pd.to_datetime(cohort_df['DATE'])


# ADD a 9 month forward window from teh date of supervision code
cohort_df['DATE_forward'] = cohort_df['DATE'] + np.timedelta64(9, 'M')
cohort_df.head()


# load dm datasets
t1_df = pd.read_csv( t1_date_file, sep=" ")
t2_df = pd.read_csv( t2_date_file, sep=" ")
a1c_df = pd.concat([pd.read_csv(A1C_FILES[0], sep="\t"), pd.read_csv(A1C_FILES[1], sep="\t")])
a1c_df = a1c_df.rename(columns={'PERSON_ID':'GRID'}, inplace=False)


# %%
# -----------
# ANALYSIS
# -----------

# select t1dm cases and not1dm_controls
dm_cohort_df = cohort_df.loc[cohort_df['GRID'].isin(t1_df['GRID'])].copy()
no_dm_cohort_df = cohort_df.loc[~cohort_df['GRID'].isin(t1_df['GRID'])].copy()

# combine w/ a1c data
dm_cohort_date_dict =  dict(zip(dm_cohort_df['GRID'], dm_cohort_df['DATE_forward']))
dm_a1c_df = a1c_df.loc[a1c_df['GRID'].isin(dm_cohort_df['GRID'])].copy()
dm_a1c_df['cohort_date'] = dm_a1c_df['GRID'].map(dm_cohort_date_dict)
dm_a1c_df['LAB_DATE'] = pd.to_datetime(dm_a1c_df['LAB_DATE'])

# keep only a1c values occuring before the earliest supervision date
dm_a1c_df['earliest_diff_lab_date'] = (dm_a1c_df['cohort_date'] -  dm_a1c_df['LAB_DATE'])/ np.timedelta64(1, 'Y')
dm_a1c_before_supervision_df = dm_a1c_df.loc[dm_a1c_df['earliest_diff_lab_date'] < 0].copy()
dm_a1c_before_supervision_df['GRID'].nunique()


# create df w/ pregnancy complications phecodes
icd9_filt_df = i9_df.loc[i9_df['CODE'].isin(icd9_preg_df['icd9'])].copy()
icd10_filt_df = i10_df.loc[i10_df['CODE'].isin(icd10_preg_df['icd10'])].copy()

# map to phecodes;
icd9_filt_df['Phecode'] = icd9_filt_df["CODE"].map(icd9_to_phe_dict)
icd10_filt_df['Phecode'] = icd10_filt_df["CODE"].map(icd10_to_phe_dict)
icd_9_10_filt_df = pd.concat([icd9_filt_df, icd10_filt_df])
icd_9_10_filt_df['DATE'] = pd.to_datetime(icd_9_10_filt_df['DATE'])

# keep only the earliest instance of each type of code
icd_9_10_filt_df.sort_values(['GRID', 'Phecode', 'DATE'], inplace=True)
_bool = ~icd_9_10_filt_df.duplicated(subset=['GRID','Phecode'], keep='first')
earliest_phe_filt_df = icd_9_10_filt_df.loc[_bool].copy()


# map one a1c to each person
earliest_phe_filt_in_dm_a1c_df = earliest_phe_filt_df[earliest_phe_filt_df['GRID'].isin(dm_a1c_before_supervision_df['GRID'])].copy()
earliest_phe_filt_in_dm_a1c_df.shape

#check that phecode occured before the cohort date
earliest_phe_filt_in_dm_a1c_df['cohort_date'] = earliest_phe_filt_in_dm_a1c_df['GRID'].map(dm_cohort_date_dict)
earliest_phe_filt_in_dm_a1c_df['cohort_diff_phe_date'] = (earliest_phe_filt_in_dm_a1c_df['cohort_date'] -  earliest_phe_filt_in_dm_a1c_df['DATE'])/ np.timedelta64(1, 'Y')

phe_before_cohort_date_df = earliest_phe_filt_in_dm_a1c_df.loc[earliest_phe_filt_in_dm_a1c_df['cohort_diff_phe_date'] > 0].copy()
phe_before_cohort_date_df.head()
phe_before_cohort_date_df['GRID'].nunique()


# append A1c closest to the phecode date


# %%

# loop through each phecode
for this_phe in phe_before_cohort_date_df['Phecode'].unique():

    this_phe_df = phe_before_cohort_date_df.query("Phecode == @this_phe")
    this_phe_df["GRID"].nunique()
    # pick a1 closest to this
    break
    this_phe_df['a1c_bin'] = pd.cut(this_phe_df['LAB_VALUE'],bins=[3,4,5,6,7,8,9,10,11])


    # dm_cohort_df