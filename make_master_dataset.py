#!/bin/python
# This script will
#
#
#
# Abin Abraham
# created on: 2021-10-17 20:51:53
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
# -----------
# paths
# -----------

BILLING_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_09_28_sd_billing_data_from_adi")
icd9_file = BILLING_DIR.joinpath("data.V01/tidy_long/long_UTI.ICD9.V01.Feb2021.tsv.gz")
icd10_file = BILLING_DIR.joinpath("data.V01/tidy_long/long_UTI.ICD10.V01.Feb2021.tsv.gz")
cpt_file = BILLING_DIR.joinpath("data.V01/tidy_long/long_UTI.CPT.V01.Feb2021.tsv.gz")
demo_file = BILLING_DIR.joinpath("data.V01/UTI.Demographics.V01.Feb2021.csv.gz")
socio_file = BILLING_DIR.joinpath("data.V01/UTI.SOCIO-ECONOMIC.V01.Feb2021.csv.gz")

supervision_icd9_10_codes_file = BILLING_DIR.joinpath("supervision_icd9_10_codes.txt")
delivery_cpt_codes_file = BILLING_DIR.joinpath("delivery_cpt_codes.txt")
delivery_icd10_codes_file = BILLING_DIR.joinpath("delivery_icd10_codes.txt")
delivery_icd9_codes_file = BILLING_DIR.joinpath("delivery_icd9_codes.txt")


DM_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_03_11_from_annika_t1_t2_dates")
t1_date_file = DM_DIR.joinpath("earliest_date_T1D_ICD.txt")
t1_cases_file = DM_DIR.joinpath("t1d_cases.txt")
t2_date_file = DM_DIR.joinpath("earliest_date_T2D_ICD.txt")
t2_cases_file = DM_DIR.joinpath("t2dm_phekb/phekb2012_t2dm_cases.txt")


LABS_DIR = Path("/dors/capra_lab/users/abraha1/data/davis_clinical_labs/selected_lab_summaries")
A1C_FILES = [LABS_DIR.joinpath('HgbA1C', "20200108_HgbA1C_4sd_outliers_removed.txt"),
             LABS_DIR.joinpath('HGB_A1C', "20200108_HGB_A1C_4sd_outliers_removed.txt")]

ALL_EGA_FILE="/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/combined_cohorts_07_12_2018/full_dataset_characterization/expanded_ega/date_of_conception_w_ega.tsv"

OUTPUT_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_10_19_master_dataset")



# %%
# -----------
# functions
# -----------
def load_inclusion_codes():
    inclusion_codes = dict()
    for key, file in zip(['supervision','delivery_icd9','delivery_icd10', 'delivery_cpt'],
                         [supervision_icd9_10_codes_file, delivery_icd9_codes_file,
                          delivery_icd10_codes_file, delivery_cpt_codes_file]):
        df_ = pd.read_csv( file, sep="\t", names=['name'], dtype=str)
        df_['code'] = df_['name'].apply(lambda x: str(x.split('-')[0]))
        if key != 'delivery_cpt':
            df_['description'] = df_['name'].apply(lambda x: str(x.split('-')[1]))
        inclusion_codes[key] = df_

    return inclusion_codes

def get_earliest_supervised_preg_date(icd9_df, icd10_df, inclusion_codes):
    # use supervision icd9/10 codes and select earliest instance of code for each patient

    icd_df = pd.DataFrame()
    for df in [icd9_df, icd10_df]:
        supervision_df = df.loc[df['CODE'].isin(inclusion_codes['supervision']['code'].values)].copy()
        supervision_df.sort_values(['GRID','DATE'], inplace=True, ascending=True)
        supervision_df.drop_duplicates(subset=['GRID'], keep='first', inplace=True)
        icd_df = icd_df.append(supervision_df)
    icd_df.sort_values(['GRID', 'DATE'], ascending=True, inplace=True)
    icd_df.drop_duplicates(subset=['GRID'], keep='first', inplace=True)
    return icd_df

def get_earliest_delivery_date(icd9_df, icd10_df, cpt_df, inclusion_codes):

    delivery_df = pd.DataFrame()
    for bill_df, key_name  in zip([icd9_df, icd10_df, cpt_df], ['delivery_icd9', 'delivery_icd10', 'delivery_cpt']):


        this_delivery_df = bill_df.loc[bill_df['CODE'].isin(inclusion_codes[key_name]['code'].values)].copy()
        this_delivery_df.sort_values(['GRID','DATE'], inplace=True, ascending=True)
        this_delivery_df.drop_duplicates(subset=['GRID'], keep='first', inplace=True)
        this_delivery_df['code_type']=key_name
        delivery_df = delivery_df.append(this_delivery_df)

    delivery_df.sort_values(['GRID', 'DATE'], ascending=True, inplace=True)
    delivery_df.drop_duplicates(subset=['GRID'], keep='first', inplace=True)

    return delivery_df

def link_supervision_and_delivery_dates(icd_df, delivery_df):

    preg_df = pd.merge(icd_df, delivery_df, on='GRID', how='inner',suffixes=['_supervision', '_delivery'] )
    preg_df['diff_date_delivery_supervision']  = preg_df['DATE_delivery'] - preg_df['DATE_supervision']
    within_45wk_bool = (preg_df['diff_date_delivery_supervision'] < np.timedelta64(45, 'W')) & (preg_df['diff_date_delivery_supervision'] > np.timedelta64(0, 'W'))
    preg_df['valid_preg_dates'] = within_45wk_bool
    valid_preg_df = preg_df.loc[preg_df['valid_preg_dates']==True,:].copy()
    valid_preg_df.drop(columns={'valid_preg_dates'}, inplace=True)
    return valid_preg_df

def load_a1c_and_format(A1C_FILES, valid_preg_delivery_df_):

    a1c_df = pd.concat([pd.read_csv(A1C_FILES[0], sep="\t"), pd.read_csv(A1C_FILES[1], sep="\t")])
    a1c_df = a1c_df.rename(columns={'PERSON_ID':'GRID'}, inplace=False)
    a1c_df['LAB_DATE'] = pd.to_datetime(a1c_df['LAB_DATE'])


    a1c_merged_df = pd.merge(a1c_df, valid_preg_delivery_df_.loc[:, ['GRID','DATE_delivery', 'DATE_supervision',]], how='left', on=['GRID'], validate='many_to_one')
    a1c_merged_df['diff_delivery_lab'] = a1c_merged_df['LAB_DATE'] - a1c_merged_df['DATE_delivery']
    a1c_merged_df['diff_supervision_lab'] = a1c_merged_df['LAB_DATE'] - a1c_merged_df['DATE_supervision']

    # a1c before delivery
    col_=['GRID','LAB_DATE', 'LAB_VALUE','DATE_delivery','DATE_supervision','diff_delivery_lab','diff_supervision_lab']
    a1c_before_delivery_df = a1c_merged_df.loc[a1c_merged_df['diff_supervision_lab'] <= np.timedelta64(0, 'D'),col_ ].copy()
    a1c_before_delivery_df.sort_values(['GRID','diff_delivery_lab'], ascending=False,inplace=True)
    closest_a1c_to_delivery_df = a1c_before_delivery_df[~a1c_before_delivery_df.duplicated(subset=['GRID'],keep='first')].copy() # pick the a1c value occuring closest to delivery date
    closest_a1c_to_delivery_df.rename(columns={'LAB_DATE':'LAB_DATE_close_to_delivery', 'LAB_VALUE':'LAB_VALUE_close_to_delivery'},inplace=True)
    flat_a1c_df  = closest_a1c_to_delivery_df.copy()

    # repeat with a1c occuring within 3 months of supervision code
    bool_ = (a1c_merged_df['diff_supervision_lab'] < np.timedelta64('3', "M")) & (a1c_merged_df['diff_supervision_lab'] > np.timedelta64('-3', "M"))
    # a1c_merged_df.loc[bool_, ['GRID','LAB_VALUE', 'LAB_DATE', 'DATE_supervision', 'diff_supervision_lab']] # check a1c distribution
    a1c_around_supervision_df = a1c_merged_df.loc[bool_, col_].copy()
    # take the absolute value of the diff_sueprvision lab
    a1c_around_supervision_df['diff_supervision_lab'] = (a1c_around_supervision_df['diff_supervision_lab']/np.timedelta64(1, 'M')).abs()
    a1c_around_supervision_df.sort_values(['GRID','diff_supervision_lab'], ascending=False,inplace=True)
    closest_a1c_to_supervision_df = a1c_around_supervision_df[~a1c_around_supervision_df.duplicated(subset=['GRID'],keep='first')].copy() # pick the a1c value occuring closest to delivery date
    closest_a1c_to_supervision_df.rename(columns={'LAB_DATE':'LAB_DATE_close_to_supervision', 'LAB_VALUE':'LAB_VALUE_close_to_supervision'},inplace=True)

    flat_a1c_all_df = pd.merge(flat_a1c_df, closest_a1c_to_supervision_df.loc[:, ['GRID','DATE_delivery','DATE_supervision','LAB_VALUE_close_to_supervision', 'LAB_DATE_close_to_supervision']], on=['GRID','DATE_delivery','DATE_supervision'])


    # add t1dm column
    grid2t1dx = dict(zip(valid_preg_delivery_df_['GRID'], valid_preg_delivery_df_['t1dm_before_supervision_preg_bool']))
    flat_a1c_all_df['t1dm_before_supervision_preg_bool'] = flat_a1c_all_df['GRID'].map(grid2t1dx)

    return flat_a1c_all_df

def plot_a1c_time_to_delivery():
    to_plot_df = flat_a1c_df.copy()
    to_plot_df['diff_delivery_lab_years']=to_plot_df['diff_delivery_lab']/np.timedelta64(1, 'Y')


    sns.set(style="ticks",  font_scale=1.0, rc={"figure.figsize": (8, 2)})
    fig, axs = plt.subplots(nrows=1, ncols=3, sharex=False, sharey=False)



    sns.boxplot(x=to_plot_df['diff_delivery_lab_years'], fliersize=1, ax=axs[0])
    axs[0].set_title("all data")
    axs[0].set_xlabel("")

    sns.boxplot(x=to_plot_df['diff_delivery_lab_years'], fliersize=1, ax=axs[1])
    axs[1].set_title("zoomed in")
    axs[1].set_xlabel("Time from A1C value to delivery\n(Years)")
    axs[1].set_xlim(-5,0)

    sns.boxplot(x=to_plot_df['diff_delivery_lab_years'], fliersize=1, ax=axs[2])
    axs[2].set_title("zoomed in")
    axs[2].set_xlabel("")
    axs[2].set_xlim(-1,0)

def make_tidy_a1c(flat_a1c_df):
    # keep a1c closest to delivery
    # keep a1c within 9 months of delivery

    cols_=['GRID', 'DATE_delivery', 'DATE_supervision', 'LAB_DATE_close_to_delivery', 'LAB_VALUE_close_to_delivery', 't1dm_before_supervision_preg_bool', 'LAB_VALUE_close_to_supervision', 'LAB_DATE_close_to_supervision']
    tidy_a1c_df = flat_a1c_df.loc[:, cols_].copy()
    tidy_a1c_df['a1c_within_9mo_of_delivery_bool'] = (tidy_a1c_df['DATE_delivery']-tidy_a1c_df['LAB_DATE_close_to_delivery']) <= np.timedelta64(9, 'M')
    tidy_a1c_df.rename(columns={'LAB_VALUE_close_to_delivery':'a1c_closest_to_delivery', 'LAB_DATE_close_to_delivery':'DATE_a1c_delivery',
                                'LAB_VALUE_close_to_supervision':'a1c_closest_to_supervision', 'LAB_DATE_close_to_supervision':'DATE_a1c_supervision',    },inplace=True)



    return tidy_a1c_df

def intersect_t1_and_t2_dm(valid_preg_df, t1_date_file, t2_date_file):
    # intersect t1dm before pregnancy
    # merge valid_preg_df with t1dm data
    # make column "t1dm_before_supervision_preg_bool" --> indicates if T1DM diagnosis before deliery


    # these are just time stamps for the diagnosis
    t1_df = pd.read_csv( t1_date_file, sep=" ", parse_dates=['min_targ_ICD_date'])
    t2_df = pd.read_csv( t2_date_file, sep=" ", parse_dates=['min_targ_ICD_date'])

    # these are t1 t2 labels used to resolve discrepancy if a pt has both dx
    t1_cases_df = pd.read_csv(t1_cases_file , sep=" ", names=['GRID','status'])
    t2_cases_df = pd.read_csv(t2_cases_file , sep=" ", names=['GRID','status'])
    assert len(set(t1_cases_df['GRID'].unique()).intersection(set(t2_cases_df['GRID'].unique()))) == 0, 't1 and t2 sets have overlapping grids'

    # remove grids w/ a t1dm diagnosis in the t2dm dataset
    t1_t2_dm_grids = set(t1_df['GRID'].unique()).intersection(set(t2_df['GRID'].unique()))
    t1_confirmed_grids = set(t1_cases_df.loc[t1_cases_df['GRID'].isin(t1_t2_dm_grids), 'GRID'].unique())
    clean_t2_df = t2_df.loc[ ~t2_df['GRID'].isin(t1_confirmed_grids), :].copy()


    # merge diabetes diagnosis w/ preg
    valid_preg_delivery_df_ = pd.merge(valid_preg_df, t1_df, on='GRID', how='left' )
    valid_preg_delivery_df_['diff_supervision_t1dm'] = valid_preg_delivery_df_['DATE_supervision'] - valid_preg_delivery_df_['min_targ_ICD_date']
    valid_preg_delivery_df_['t1dm_before_supervision_preg_bool'] = valid_preg_delivery_df_['diff_supervision_t1dm'] >= np.timedelta64(0,'D')
    valid_preg_delivery_df_.rename(columns={'min_targ_ICD_date':'DATE_t1d_dx'},inplace=True)


    #t1_t2_dm_df = pd.merge(valid_preg_delivery_df_, clean_t2_df, on='GRID', how='left')
    t1_t2_dm_df = pd.merge(valid_preg_delivery_df_, clean_t2_df, on='GRID', how='left')
    t1_t2_dm_df['diff_supervision_t2dm'] = t1_t2_dm_df['DATE_supervision'] - t1_t2_dm_df['min_targ_ICD_date']
    t1_t2_dm_df['t2dm_before_supervision_preg_bool'] = t1_t2_dm_df['diff_supervision_t2dm'] >= np.timedelta64(0,'D')
    t1_t2_dm_df.rename(columns={'min_targ_ICD_date':'DATE_t2d_dx'},inplace=True)


    t1_t2_dm_df['dm_before_supervision_code']='no_dm'
    t1_t2_dm_df.loc[(t1_t2_dm_df['t1dm_before_supervision_preg_bool']==True) & (t1_t2_dm_df['t2dm_before_supervision_preg_bool']==True) , 'dm_before_supervision_code'] = 't1dm_t2dm'
    t1_t2_dm_df.loc[(t1_t2_dm_df['t1dm_before_supervision_preg_bool']==True) & (t1_t2_dm_df['t2dm_before_supervision_preg_bool']==False), 'dm_before_supervision_code'] = 't1dm'
    t1_t2_dm_df.loc[(t1_t2_dm_df['t1dm_before_supervision_preg_bool']==False) & (t1_t2_dm_df['t2dm_before_supervision_preg_bool']==True ), 'dm_before_supervision_code'] = 't2dm'
    t1_t2_dm_df.drop(['diff_supervision_t1dm','diff_supervision_t2dm'], axis=1, inplace=True)


    # censor dx dates if they occured after supervision
    t1_t2_dm_df.loc[t1_t2_dm_df['dm_before_supervision_code']=='no_dm','DATE_t1d_dx'] = np.nan
    t1_t2_dm_df.loc[t1_t2_dm_df['dm_before_supervision_code']=='no_dm','DATE_t2d_dx'] = np.nan
    t1_t2_dm_df.loc[t1_t2_dm_df['dm_before_supervision_code']=='t1dm','DATE_t2d_dx'] = np.nan
    t1_t2_dm_df.loc[t1_t2_dm_df['dm_before_supervision_code']=='t2dm','DATE_t1d_dx'] = np.nan


    return t1_t2_dm_df

def merge_w_ega(ALL_EGA_FILE, delivery_a1c_t1d_df):

    all_ega_df = pd.read_csv(ALL_EGA_FILE, sep="\t")

    all_ega_df['DATE_delivery']= all_ega_df['GRID'].map(dict(zip(delivery_a1c_t1d_df['GRID'], delivery_a1c_t1d_df['DATE_delivery'])))
    all_ega_df['DATE_delivery'] = pd.to_datetime(all_ega_df['DATE_delivery'])
    all_ega_df['delivery_date'] = pd.to_datetime(all_ega_df['delivery_date'])
    all_ega_df['diff_delivery'] = ((all_ega_df['delivery_date'] - all_ega_df['DATE_delivery'])/np.timedelta64(1,'D')).abs()


    linked_ega_df = all_ega_df[(all_ega_df['diff_delivery'] < 7)].copy()
    linked_ega_df.sort_values(['GRID','diff_delivery'], ascending=True, inplace=True)
    linked_ega_df.drop_duplicates(subset=['GRID'], keep='first',inplace=True)


    delivery_a1c_t1d_ega_df = pd.merge(delivery_a1c_t1d_df, linked_ega_df.loc[:, ['GRID','closest_ega_DATE','closest_ega','DATE_delivery']], on=['GRID','DATE_delivery'], how='left', validate='one_to_one')

    return delivery_a1c_t1d_ega_df
# %%
# -----------
# main
# -----------

###
### OVERVIEW
###

# 1) Select women with supervision codes
# 2) Select women w/ evidence of delivery
# 3) Select earliest pregnancy based on the supervision code
# 4) Syn supervision code and delivery code temporally
#       * from the earliest supervision code, identify the closest delivery code
#       * keep only deliveries that occured within 9 months of the supervision code
# 5) intersect with T1D diagnose occuring before 'earliest supervision' code
# 6) intersect with A1C occuring before the delivery date
# 7) intersect gestational age data
# 8) intersect race, DOB,


# %%

# load billing code data
icd9_df = pd.read_csv( icd9_file, sep="\t", names=['GRID','DATE','CODE'], dtype={'GRID':str,'CODE':str}, parse_dates=['DATE'])
icd10_df = pd.read_csv( icd10_file, sep="\t", names=['GRID','DATE','CODE'], dtype={'GRID':str,'CODE':str}, parse_dates=['DATE'])
cpt_df = pd.read_csv( cpt_file, sep="\t", names=['GRID','DATE','CODE'], dtype={'GRID':str,'CODE':str}, parse_dates=['DATE'])


# load supervision and delivery codes
inclusion_codes = load_inclusion_codes()
inclusion_codes.keys()

# get earliest supervision date and delivery date in each EHR
icd_df = get_earliest_supervised_preg_date(icd9_df, icd10_df, inclusion_codes)
delivery_df = get_earliest_delivery_date(icd9_df, icd10_df, cpt_df, inclusion_codes)
del icd9_df, icd10_df

# link supervision to dlievery date temporally
valid_preg_df = link_supervision_and_delivery_dates(icd_df, delivery_df)
valid_preg_df.shape


# intersect  t1dm data
valid_preg_delivery_df_ = intersect_t1_and_t2_dm(valid_preg_df, t1_date_file, t2_date_file)
pd.value_counts(valid_preg_delivery_df_['dm_before_supervision_code'])

# intersect a1c data
flat_a1c_df = load_a1c_and_format(A1C_FILES, valid_preg_delivery_df_)
tidy_a1c_df = make_tidy_a1c(flat_a1c_df)
tidy_a1c_df.shape




# %%
# -----------
# merge preganncy w/ other variables
# -----------
# merge a1c data with preganncy data
delivery_a1c_t1d_df = pd.merge(valid_preg_delivery_df_.loc[:, ['GRID','DATE_supervision','DATE_delivery','DATE_t1d_dx', 'DATE_t2d_dx', 'dm_before_supervision_code']],
                              tidy_a1c_df.drop('t1dm_before_supervision_preg_bool', axis=1), on=['GRID','DATE_supervision','DATE_delivery'], how='left', validate='one_to_one')
delivery_a1c_t1d_df['a1c_within_9mo_of_delivery_bool'].fillna(False, inplace=True)


# %%
### merge ega data
delivery_a1c_t1d_ega_df = merge_w_ega(ALL_EGA_FILE, delivery_a1c_t1d_df)

# %%
demo_df = pd.read_csv( demo_file, sep=",")
ses_df = pd.read_csv( socio_file, sep=",")

# %%
delivery_a1c_t1d_ega_demo_df = pd.merge(delivery_a1c_t1d_ega_df, demo_df.loc[:, ['GRID','Race','Ethnicity','birthDate']], on=['GRID'], how='left', validate='one_to_one')
delivery_a1c_t1d_ega_demo_ses_df = pd.merge(delivery_a1c_t1d_ega_demo_df, ses_df.loc[:, ['GRID','DEP_INDEX', 'MEDIAN_INCOME']], on=['GRID'], how='left', validate='one_to_one')
delivery_a1c_t1d_ega_demo_ses_df['closest_ega_DATE'] = pd.to_datetime(delivery_a1c_t1d_ega_demo_ses_df['closest_ega_DATE'])

delivery_a1c_t1d_ega_demo_ses_df.dtypes

# %%
delivery_a1c_t1d_ega_demo_ses_df.head()
delivery_a1c_t1d_ega_demo_ses_df.to_csv(OUTPUT_DIR.joinpath(f'{DATE}_master_data.tsv'), sep="\t", index=False)

# %%
# -----------
# internal checks on the data set
# -----------
# check that supervision occured before delivery
((delivery_a1c_t1d_ega_demo_ses_df['DATE_supervision'] - delivery_a1c_t1d_ega_demo_ses_df['DATE_delivery']) < np.timedelta64(0, 'D')).all()

# make sure a1c is before DATE of delivery
bool_a1c = ~delivery_a1c_t1d_ega_demo_ses_df['DATE_a1c'].isna()
((delivery_a1c_t1d_ega_demo_ses_df.loc[bool_a1c, 'DATE_a1c'] - delivery_a1c_t1d_ega_demo_ses_df.loc[bool_a1c, 'DATE_delivery']) <= np.timedelta64(0,'D')).all()

# check that t1d and t2d dx occured before supervision
bool_t1d = ~delivery_a1c_t1d_ega_demo_ses_df['DATE_t1d_dx'].isna()
delivery_a1c_t1d_ega_demo_ses_df['temp_diff_d1'] = (delivery_a1c_t1d_ega_demo_ses_df.loc[bool_t1d, 'DATE_t1d_dx'] - delivery_a1c_t1d_ega_demo_ses_df.loc[bool_t1d, 'DATE_delivery'])
(delivery_a1c_t1d_ega_demo_ses_df.loc[bool_t1d, 'temp_diff_d1']<= np.timedelta64(0, 'D')).all()

bool_t2d = ~delivery_a1c_t1d_ega_demo_ses_df['DATE_t2d_dx'].isna()
delivery_a1c_t1d_ega_demo_ses_df['temp_diff_d2'] = (delivery_a1c_t1d_ega_demo_ses_df.loc[bool_t2d, 'DATE_t2d_dx'] - delivery_a1c_t1d_ega_demo_ses_df.loc[bool_t2d, 'DATE_delivery'])
(delivery_a1c_t1d_ega_demo_ses_df.loc[bool_t2d, 'temp_diff_d2']<= np.timedelta64(0, 'D')).all()
delivery_a1c_t1d_ega_demo_ses_df.drop(['temp_diff_d1, temp_diff_d2'], axis=1, inplace=True)
