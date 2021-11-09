#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2021-03-11 08:55:05
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')


LABS_DIR = Path("/dors/capra_lab/users/abraha1/data/davis_clinical_labs/biovu_202006_sdwide_pull/selected_lab_summaries_10042020")
OUTPUT_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_03_11_cleaner_dm_labs")

# DELIVERY_FILE = "/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/combined_cohorts_07_12_2018/full_dataset_characterization/est_delivery_date_at_least_one_icd_cpt_ega.tsv"
DM_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_03_11_from_annika_t1_t2_dates")
t1_date_file = DM_DIR.joinpath("earliest_date_T1D_ICD.txt")
t2_date_file = DM_DIR.joinpath("earliest_date_T2D_ICD.txt")


DELIVERY_FILES = "/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/combined_cohorts_07_12_2018/full_dataset_characterization/updated_ega_2020_06_16/earliest_delivery_with_date_of_conception_w_ega_updated_EGA_2020_06_17.tsv"

# %%
# -----------
# main
# -----------

delivery_df = pd.read_csv(DELIVERY_FILEs, sep="\t")
delivery_df = delivery_df.loc[:, ['GRID','consensus_delivery','consensus_label', 'delivery_id']].copy()
delivery_df['consensus_delivery'] = pd.to_datetime(delivery_df['consensus_delivery'])
delivery_grids = delivery_df['GRID'].unique()

t1_df = pd.read_csv( t1_date_file, sep="\s")
t2_df = pd.read_csv( t2_date_file, sep="\s")

t1_merge = set(delivery_df['GRID'].unique()).intersection(set(t1_df['GRID']))
t2_merge = set(delivery_df['GRID'].unique()).intersection(set(t2_df['GRID']))

grids_to_keep= t1_merge.union(t2_merge).union(set(delivery_df['GRID'].unique()))
len(grids_to_keep)
# %%

for lab_file in LABS_DIR.glob("*A1C*.txt"):

    lab_df = pd.read_csv( lab_file, sep="\t")
    filt_df = lab_df[lab_df['PERSON_ID'].isin(delivery_grids)].copy()

    # assign labs based on delivery date
    filt_df

    lab_df.shape
    filt_df.shape
    break
# %%
lab_df.shape
lab_df.head()

filt_df