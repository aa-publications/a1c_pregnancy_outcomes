#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2021-09-27 21:06:12

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path



import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

%matplotlib inline
%config InlineBackend.figure_format='retina'


DATE = datetime.now().strftime('%Y-%m-%d')


### PATH
ICD_DIR=Path("/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/combined_cohorts_07_12_2018/full_dataset")
icd9_file = ICD_DIR.joinpath('full_ICD9_cohort.tsv')
icd10_file = ICD_DIR.joinpath('full_ICD10_cohort.tsv')


DM_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_03_11_from_annika_t1_t2_dates")
t1_date_file = DM_DIR.joinpath("earliest_date_T1D_ICD.txt")
t2_date_file = DM_DIR.joinpath("earliest_date_T2D_ICD.txt")


LABS_DIR = Path("/dors/capra_lab/users/abraha1/data/davis_clinical_labs/selected_lab_summaries")
A1C_FILES = [LABS_DIR.joinpath('HgbA1C', "20200108_HgbA1C_4sd_outliers_removed.txt"), LABS_DIR.joinpath('HGB_A1C', "20200108_HGB_A1C_4sd_outliers_removed.txt")]

PHE_DIR = Path("/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/phecodes/2021_09_27_phecodev2")
PREG_PHE_FILE = Path("/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/phecodes/2021_09_27_phecodev2/phecode_strings_V2_preg_chapter.csv")

OUTPUT_DIR= Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_09_27_preg_complications_by_a1c")

# %%
# -----------
# MAIN
# -----------

icd9_df = pd.read_csv( icd9_file, sep="\t")
icd10_df = pd.read_csv( icd10_file, sep="\t")
t1_df = pd.read_csv( t1_date_file, sep=" ")
t2_df = pd.read_csv( t2_date_file, sep=" ")
a1c_df = pd.concat([pd.read_csv(A1C_FILES[0], sep="\t"), pd.read_csv(A1C_FILES[1], sep="\t")])
a1c_df = a1c_df.rename(columns={'PERSON_ID':'GRID'}, inplace=False)

preg_phe_df = pd.read_csv( PREG_PHE_FILE, sep=",", names = ["phecode","phecode_string","phecode_category","sex","ICD10_only","phecode_top","leaf"])
phe_to_label_dict = dict(zip(preg_phe_df['phecode'], preg_phe_df['phecode_string']))
icd9_to_phe_df = pd.read_csv(  PHE_DIR.joinpath("ICD9_to_phecode_V2.csv"), sep=",")
icd10_to_phe_df = pd.read_csv( PHE_DIR.joinpath("ICD10_to_phecode_V2.csv"), sep=",")

# identify icd9 or icd10 pregnancy codes of interest
icd9_preg_df = icd9_to_phe_df.loc[icd9_to_phe_df['phecode'].isin(preg_phe_df['phecode'])].copy()
icd10_preg_df = icd10_to_phe_df.loc[icd10_to_phe_df['phecode'].isin(preg_phe_df['phecode'])].copy()



icd9_to_phe_dict = dict(zip(icd9_preg_df['icd9'], icd9_preg_df['phecode']))
icd10_to_phe_dict = dict(zip(icd10_preg_df['icd10'], icd10_preg_df['phecode']))

# keep only preganncy codes in icd9 and icd10 datasets
icd9_filt_df = icd9_df.loc[icd9_df['ICD'].isin(icd9_preg_df['icd9'])].copy()
icd10_filt_df = icd10_df.loc[icd10_df['ICD'].isin(icd10_preg_df['icd10'])].copy()

# keep only the earliest instance of each type of code
icd9_filt_df['Phecode'] = icd9_filt_df["ICD"].map(icd9_to_phe_dict)
icd10_filt_df['Phecode'] = icd10_filt_df["ICD"].map(icd10_to_phe_dict)

icd_9_10_filt_df = pd.concat([icd9_filt_df, icd10_filt_df])
icd_9_10_filt_df['Date'] = pd.to_datetime(icd_9_10_filt_df['Date'])

icd_9_10_filt_df.sort_values(['GRID', 'Phecode', 'Date'], inplace=True)
_bool = ~icd_9_10_filt_df.duplicated(subset=['GRID','Phecode'], keep='first')
earliest_phe_filt_df = icd_9_10_filt_df.loc[_bool].copy()
earliest_phe_filt_df.shape

# %%
earliest_phe_filt_df['Phecode'].nunique()
for i, this_phe in enumerate(earliest_phe_filt_df['Phecode'].unique()):

    this_phe_df = earliest_phe_filt_df.query("Phecode == @this_phe")
    grid_to_earliest_date_dict = dict(zip(this_phe_df['GRID'], this_phe_df['Date']))

    cases_filt_a1c_df = a1c_df.loc[a1c_df['GRID'].isin(this_phe_df['GRID'])].copy()
    controls_filt_a1c_df = a1c_df.loc[~a1c_df['GRID'].isin(this_phe_df['GRID'])].copy()

    # keep a1c values present only before phecode date
    cases_filt_a1c_df['earliest_date'] = cases_filt_a1c_df['GRID'].map(grid_to_earliest_date_dict)
    cases_filt_a1c_df['LAB_DATE'] = pd.to_datetime(cases_filt_a1c_df['LAB_DATE'])
    cases_filt_a1c_df['earliest_date'] = pd.to_datetime(cases_filt_a1c_df['earliest_date'])
    cases_filt_a1c_df['earliest_diff_lab_date'] = (cases_filt_a1c_df['earliest_date'] -  cases_filt_a1c_df['LAB_DATE'])/ np.timedelta64(1, 'Y')
    before_case_a1c_df = cases_filt_a1c_df.loc[cases_filt_a1c_df['earliest_diff_lab_date'] < 0].copy()
    before_case_a1c_df.sort_values(['GRID','earliest_diff_lab_date'], inplace=True)
    cases_closest_a1c_df= before_case_a1c_df[~before_case_a1c_df.duplicated(subset=['GRID'],keep='last')].copy()

    cases_closest_a1c_df
    # cases_closest_a1c_df['a1c_bin'] = pd.cut(cases_closest_a1c_df['LAB_VALUE'],bins=12)
    cases_closest_a1c_df['a1c_bin'] = pd.cut(cases_closest_a1c_df['LAB_VALUE'],bins=[3,4,5,6,7,8,9,10,11])

    sns.set(style="ticks",  font_scale=1.0, rc={"figure.figsize": (6, 4)})
    fig, ax = plt.subplots()
    sns.countplot(x="a1c_bin", data=cases_closest_a1c_df, ax=ax, color='gray')
    figtitle=f"{this_phe} {phe_to_label_dict[this_phe]}"
    ax.set_title(f"{this_phe}_{phe_to_label_dict[this_phe]}")
    ax.set_ylabel("number of individuals")
    sns.despine(ax=ax, top=True)
    ax.set_xticklabels([aa.get_text() for aa in ax.get_xticklabels()], rotation = 330, ha="left")
    for p in ax.patches:
        ht = p.get_height()
        if np.isnan(ht):
            ht = 0
        ax.text(p.get_x() + p.get_width()/2., ht, '%d' % int(ht),
                fontsize=12, color='black', ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR.joinpath(f'{DATE}_{figtitle.replace(" ", "_")}.pdf'))
    plt.clf()
    # %%


