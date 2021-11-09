#!/bin/python
# This script will ...
#   * loads cohort dataset
#   * loads icd9, icd10 and convert to phecode; keep only adverse pregnancy outcomes
#   * keep only adverse outcomes occuring 9months before delivery date and 3 mo. after delivery date
#   * also include specific CPT codes as outcomes
#
#
# Abin Abraham
# created on: 2021-10-22 22:06:54

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

DATE = datetime.now().strftime('%Y-%m-%d')



import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

%matplotlib inline

# %%

BILLING_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_09_28_sd_billing_data_from_adi")
icd9_file = BILLING_DIR.joinpath("data.V01/tidy_long/long_UTI.ICD9.V01.Feb2021.tsv.gz")
icd10_file = BILLING_DIR.joinpath("data.V01/tidy_long/long_UTI.ICD10.V01.Feb2021.tsv.gz")
cpt_file = BILLING_DIR.joinpath("data.V01/tidy_long/long_UTI.CPT.V01.Feb2021.tsv.gz")


# master_file="/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_10_19_master_dataset/2021-10-19_master_data.tsv"
master_file="/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_10_19_master_dataset/2021-10-28_master_data.tsv"

PHE_DIR = Path("/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/phecodes/2021_09_27_phecodev2")
PREG_PHE_FILE = Path("/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/phecodes/2021_09_27_phecodev2/phecode_strings_V2_preg_chapter.csv")
OUTCOME_ANNO_FILE = Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/mcnew_outcomes.xlsx")

OUTPUT_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_10_19_master_dataset")

# -----------
# functions
# -----------
def find_delivery_outcomes(filt_icd9_df, filt_icd10_df, filt_cpt_df, icd9_to_phe_dict, icd10_to_phe_dict, preg_phe_df):
    filt_icd9_df['phecode']= filt_icd9_df['CODE'].map(icd9_to_phe_dict)
    nona_icd9_df = filt_icd9_df.loc[~filt_icd9_df['phecode'].isna()].copy() # remove NaN
    filt_icd10_df['phecode']= filt_icd10_df['CODE'].map(icd10_to_phe_dict)
    nona_icd10_df = filt_icd10_df.loc[~filt_icd10_df['phecode'].isna()].copy() # remove NaN


    all_phe_only_df = pd.concat([nona_icd9_df, nona_icd10_df])
    phe_only_df = all_phe_only_df[all_phe_only_df['phecode'].isin(preg_phe_df['phecode'].unique())].copy()

    # add specific cpt codes too
    filt_cpt_df['phecode'] = filt_cpt_df['CODE']
    phe_cpt_only_df = pd.concat([filt_cpt_df, phe_only_df])

    phe_cpt_only_df['DATE_delivery'] = phe_cpt_only_df['GRID'].map(cohort2delivery_dict)
    phe_cpt_only_df['delivery_plus_3mo'] = phe_cpt_only_df['DATE_delivery'] + np.timedelta64(3, 'M')
    phe_cpt_only_df['delivery_minus_9mo'] = phe_cpt_only_df['DATE_delivery'] - np.timedelta64(9, 'M')
    phe_cpt_only_df['diff_delivery_phe'] = phe_cpt_only_df['delivery_plus_3mo'] - phe_cpt_only_df['DATE']

    delivery_outcomes_df = phe_cpt_only_df[(phe_cpt_only_df['DATE'] < phe_cpt_only_df['delivery_plus_3mo']) & (phe_cpt_only_df['DATE'] > phe_cpt_only_df['delivery_minus_9mo'])].copy()

    return delivery_outcomes_df

# %%
# -----------
# main
# -----------

# load cohort data
cohort_df = pd.read_csv( master_file, sep="\t", dtype={'GRID':str}, parse_dates=['DATE_delivery'])
cohort2delivery_dict = dict(zip(cohort_df['GRID'], cohort_df['DATE_delivery']))


# load phecodes
cohort_anno_df = pd.read_excel(OUTCOME_ANNO_FILE, dtype={'phecode':str})
cohort_anno_df.rename(columns={'Vascular ':'Vascular'}, inplace=True)
preg_phe_df = pd.read_csv( PREG_PHE_FILE, sep=",", names = ["phecode","phecode_string","phecode_category","sex","ICD10_only","phecode_top","leaf"], dtype={'phecode':'str'})
phe_to_label_dict = dict(zip(preg_phe_df['phecode'], preg_phe_df['phecode_string']))
phe_icd9_preg_df = pd.read_csv(  PHE_DIR.joinpath("ICD9_to_phecode_V2.csv"), sep=",", dtype='str')
phe_icd10_preg_df = pd.read_csv( PHE_DIR.joinpath("ICD10_to_phecode_V2.csv"), sep=",", dtype='str')


icd9_to_phe_dict = dict(zip(phe_icd9_preg_df['icd9'], phe_icd10_preg_df['phecode']))
icd10_to_phe_dict = dict(zip(phe_icd10_preg_df['icd10'], phe_icd10_preg_df['phecode']))

# load billing code data
icd9_df = pd.read_csv( icd9_file, sep="\t", names=['GRID','DATE','CODE'], dtype={'GRID':str,'CODE':str}, parse_dates=['DATE'])
icd10_df = pd.read_csv( icd10_file, sep="\t", names=['GRID','DATE','CODE'], dtype={'GRID':str,'CODE':str}, parse_dates=['DATE'])
cpt_df = pd.read_csv( cpt_file, sep="\t", names=['GRID','DATE','CODE'], dtype={'GRID':str,'CODE':str}, parse_dates=['DATE'])

# keep grids in the cohort
filt_grids = lambda df, cohort_df: df.loc[df['GRID'].isin(cohort_df['GRID'].unique())]
filt_icd9_df = filt_grids(icd9_df, cohort_df).copy()
filt_icd10_df = filt_grids(icd10_df, cohort_df).copy()
filt_cpt_df = filt_grids(cpt_df, cohort_df).copy()
filt_icd9_df['code']= 'icd9'
filt_icd10_df['code']= 'icd10'
filt_cpt_df['code']= 'cpt'


# keep only specific cpt codes
# 76811 = anatomy u/s
# 90714 = TdAP (preservative included)
# 90715 = TdAP (preservative free)
filt_cpt_df = filt_cpt_df[filt_cpt_df['CODE'].isin(['90714', '90715', '76811'])].copy()
filt_cpt_df.shape
filt_cpt_df['CODE'].unique()

# keep only
delivery_outcomes_df = find_delivery_outcomes(filt_icd9_df, filt_icd10_df, filt_cpt_df, icd9_to_phe_dict, icd10_to_phe_dict, preg_phe_df)
delivery_outcomes_df['phecode'] = delivery_outcomes_df['phecode'].map(str)

assert delivery_outcomes_df.isna().sum().sum() == 0, 'na is present in the df'
grids_to_add = set(cohort_df['GRID'].unique())- set(delivery_outcomes_df['GRID'].unique())
len(grids_to_add)

# all_grids_df = pd.concat([pd.DataFrame({'GRID':list(grids_to_add)}), delivery_outcomes_df], axis=0)
outcomes_wide_df = delivery_outcomes_df.groupby(['GRID','phecode']).size().reset_index().pivot(index='GRID', columns="phecode", values=0).fillna(0)
outcomes_wide_df.reset_index(inplace=True)
cohort_outcomes_df = pd.concat([pd.DataFrame({'GRID':list(grids_to_add)}), outcomes_wide_df]).fillna(0)

cohort_outcomes_df.loc[:, cohort_outcomes_df.columns != 'GRID'] = cohort_outcomes_df.loc[:, cohort_outcomes_df.columns != 'GRID'].applymap(int)
# cohort_outcomes_df.to_csv(OUTPUT_DIR.joinpath(f"{DATE}_cohort_outcomes.tsv"), sep="\t", index=False)
cohort_outcomes_df.drop(['76811'], axis=1, inplace=True)

# %%
# -----------
# combine with predictors
# -----------

# create composite outcomes
comp_outcomes = ['GESTATION','DELIVERY','Infections', 'Fetal_anomalies','Vascular']
for outcome_set in comp_outcomes:
    set_ = cohort_anno_df.loc[cohort_anno_df[outcome_set] == "T", 'phecode'].values
    cohort_outcomes_df[outcome_set]  = cohort_outcomes_df.loc[:, cohort_outcomes_df.columns.isin(set_)].any(1)
    # print(pd.value_counts(cohort_outcomes_df[outcome_set]))


# tidy up predictors
cohort_df['birthDate'] = pd.to_datetime(cohort_df['birthDate'])
cohort_df['DATE_supervision'] = pd.to_datetime(cohort_df['DATE_supervision'])
cohort_df['age_at_supervision_years'] = (cohort_df['DATE_supervision'] - cohort_df['birthDate'])/np.timedelta64(1, 'Y')
cohort_df['dm_outcome'] = 'no_t1dm'
cohort_df.loc[cohort_df['dm_before_supervision_code']=='t1dm', 'dm_outcome'] = 't1d'


keep_pred_cols = ['GRID','dm_outcome', 'a1c_closest_to_delivery', 'a1c_closest_to_supervision', 'age_at_supervision_years', 'Race','DEP_INDEX', 'DATE_supervision','DATE_delivery', 'DATE_t1d_dx', 'DATE_a1c_delivery',  'DATE_a1c_supervision']
tidy_pred_df = cohort_df.loc[:, keep_pred_cols].copy()

# tidy up outcomes
outcomes_cols = comp_outcomes + list(cohort_outcomes_df.columns[~cohort_outcomes_df.columns.isin( comp_outcomes +['GRID'])])
final_df = pd.merge(tidy_pred_df, cohort_outcomes_df.loc[:, ['GRID']+outcomes_cols], on='GRID', how='inner', validate='one_to_one')
final_df.to_csv(OUTPUT_DIR.joinpath(f'{DATE}_predictors_outcomes.tsv'), sep="\t", index=False, na_rep="NaN")


final_df['DEP_INDEX'].unique()
final_df['dm_outcome'].unique()
print('\n'.join(list(final_df.columns)))
# %%
phe_to_label_dict = dict(zip(preg_phe_df['phecode'], preg_phe_df['phecode_string']))


#
prev = lambda x: (outcomes_wide_df[x].sum()/outcomes_wide_df.shape[0])*100
prev_val = [np.round(prev(x),2) for x in outcomes_wide_df.columns  if not x.startswith('GRID')]
prev_code = [x for x in outcomes_wide_df.columns if not x.startswith('GRID')]
prev_df = pd.DataFrame({'phecode':prev_code, 'prev_percent': prev_val})


prev_df['label'] = prev_df['phecode'].map(phe_to_label_dict)
prev_df.sort_values('prev_percent', ascending=False, inplace=True)
prev_df.to_csv(OUTPUT_DIR.joinpath(f"{DATE}_prevalence_phe_outcome_labels.tsv"),sep="\t", index=False)

