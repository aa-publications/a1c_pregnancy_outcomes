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


# %%
# -----------
# MAIN
# -----------

#load a1c
a1c_df = pd.DataFrame()
for a1c_file in A1C_FILES:
    this_df = pd.read_csv(a1c_file, sep="\t")
    a1c_df = a1c_df.append(this_df)


# a1c_df.columns= [x[1:-1] for x in a1c_df.columns]
# a1c_df['GRID']= a1c_df['GRID'].apply(lambda x: x[1:-1])

a1c_df['PERSON_ID'].nunique()

# %%
# load deliveries
delivery_df = pd.read_csv(DELIVERY_FILE, sep="\t")
delivery_df = delivery_df.loc[:, ['GRID','consensus_delivery','consensus_label', 'delivery_id']].copy()
delivery_df['consensus_delivery'] = pd.to_datetime(delivery_df['consensus_delivery'])
first_delivery_df = delivery_df[~delivery_df.duplicated(['GRID'],keep='first')].copy()



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

# t1 and t2 overlap
t1_and_t2_grids = set(t2_df['GRID']).intersection(set(t1_df['GRID']))
len(t1_and_t2_grids)
# %%
# merge delivery datat with diabetes set
first_delivery_t1_df = pd.merge(first_delivery_df, t1_df, on='GRID', how='left')
first_delivery_dm_df = pd.merge(first_delivery_t1_df, t2_df, on='GRID', how='left')

# estimate when the diabetes diagnosis was entered into teh chart
first_delivery_dm_df['t1_date_to_delivery'] = first_delivery_dm_df.apply(lambda x: t1dm_date_to_delivery(x), axis=1)
first_delivery_dm_df['t2_date_to_delivery'] = first_delivery_dm_df.apply(lambda x: t2dm_date_to_delivery(x), axis=1)

# tabulate when teh diabetes diagnosis entered teh chart
t2_count_df = pd.value_counts(first_delivery_dm_df['t2_date_to_delivery']).reset_index()
t1_count_df = pd.value_counts(first_delivery_dm_df['t1_date_to_delivery']).reset_index()


# %%
sns.set(style="ticks",  font_scale=1.2, rc={"figure.figsize": (5, 3)})
fig ,ax = plt.subplots()
sns.barplot(y='index', x='t1_date_to_delivery', data=t1_count_df)
ax.set_xlabel('Number of Women')
ax.set_ylabel('T1DM Dx relative to delivery')
ax.set_title("Women with ≥1 delivery + T1 DM")
# plt.savefig(OUTPUT_DIR.joinpath(f"{DATE}_counts_of_t1d_dx_relative_to_first_delivery.pdf"))

# %%
sns.set(style="ticks",  font_scale=1.2, rc={"figure.figsize": (5, 3)})
fig ,ax = plt.subplots()
sns.barplot(y='index', x='t2_date_to_delivery', data=t2_count_df)
ax.set_xlabel('Number of Women')
ax.set_ylabel('T2DM Dx relative to delivery')
ax.set_title("Women with ≥1 delivery + T2 DM")
# plt.savefig(OUTPUT_DIR.joinpath(f"{DATE}_counts_of_t2d_dx_relative_to_first_delivery.pdf"))

# %%
sec