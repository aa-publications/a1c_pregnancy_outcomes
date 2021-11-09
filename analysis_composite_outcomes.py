#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2021-10-23 11:28:23


import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
# np.random.seed(123)
DATE = datetime.now().strftime('%Y-%m-%d')


DATA_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_10_19_master_dataset")
predictors_file = DATA_DIR.joinpath('2021-10-19_master_data.tsv')
outcomes_file = DATA_DIR.joinpath(f"2021-10-23_cohort_outcomes.tsv")
composite_file = DATA_DIR.joinpath("outcomes_w_include-exclude.xlsx")


PREG_PHE_FILE = Path("/dors/capra_lab/users/abraha1/projects/PTB_phenotyping/data/phecodes/2021_09_27_phecodev2/phecode_strings_V2_preg_chapter.csv")

OUTPUT_DIR=Path("/dors/capra_lab/users/abraha1/prelim_studies/2020_04_clin_labs/data/a1c_in_preg/2021_10_19_master_dataset/regression")

# %%
# -----------
# functions
# -----------


# %%
# -----------
# main
# -----------
pred_df = pd.read_csv( predictors_file, sep="\t", parse_dates=['DATE_supervision','birthDate'])
outcome_df = pd.read_csv( outcomes_file, sep="\t")
composite_df = pd.read_excel(composite_file, dtype={'phecode':'str'})

keep_outcomes = composite_df.loc[composite_df['Annotate']!='Exclude', 'phecode'].values
keep_mat_outcomes = composite_df.loc[composite_df['Annotate']=='maternal', 'phecode'].values
keep_fetal_outcomes = composite_df.loc[composite_df['Annotate']=='fetal', 'phecode'].values

# create proxy labels compatible with regression function
outcomes_var_og  = outcome_df.columns[outcome_df.columns != 'GRID'].values
proxy_outcome_labels = [f"V{x}" for x in np.arange(1,len(outcomes_var_og)+1)]
og_2_mod_vars = dict(zip(outcomes_var_og, proxy_outcome_labels))
outcome_df.rename(columns=og_2_mod_vars, inplace=True)


preg_phe_df = pd.read_csv( PREG_PHE_FILE, sep=",", names = ["phecode","phecode_string","phecode_category","sex","ICD10_only","phecode_top","leaf"], dtype={'phecode':'str'})
preg_phe_df['outcome_label'] = preg_phe_df['phecode'].map(og_2_mod_vars)


phe_to_label_dict = dict(zip(preg_phe_df['phecode'], preg_phe_df['phecode_string']))
outcome2phe_string= dict(zip(preg_phe_df['outcome_label'], preg_phe_df['phecode_string']))
outcome2phecode= dict(zip(preg_phe_df['outcome_label'], preg_phe_df['phecode']))
phe2outcome= dict(zip(preg_phe_df['phecode'], preg_phe_df['outcome_label']))


# %%
###
###    get predictors
###

# format predictors
clean_pred_df = pred_df.copy()

# keep only A1C within 9 mo
clean_pred_df['a1c'] = np.nan
clean_pred_df.loc[clean_pred_df['a1c_within_9mo_of_delivery_bool']==True, 'a1c'] = clean_pred_df.loc[clean_pred_df['a1c_within_9mo_of_delivery_bool']==True, 'a1c_closest_to_delivery']

# set a1c to 0 if no diabetes dx
clean_pred_df.loc[clean_pred_df['dm_before_supervision_code']=='no_dm', 'a1c'] =0

# remove people w/ both t1 and t2
clean_pred_df = clean_pred_df.loc[clean_pred_df['dm_before_supervision_code'] != "t1dm_t2dm", :].copy()
pd.value_counts(clean_pred_df['dm_before_supervision_code'])

clean_pred_df['dm_outcome'] = 0
clean_pred_df.loc[ (clean_pred_df['dm_before_supervision_code'] == "t1dm"), 'dm_outcome']=1

# calc age at supervision
clean_pred_df['age_at_supervision'] = (clean_pred_df['DATE_supervision'] - clean_pred_df['birthDate'])/np.timedelta64(1, 'Y')

# keep only columns of interest
pred_cols = ['GRID','dm_outcome','a1c','closest_ega','Race','DEP_INDEX','age_at_supervision']
clean_pred_df = clean_pred_df[pred_cols].copy()


###
###    get outcomes
###

# convert outcomes to binary
outcome_df.loc[:, outcome_df.columns!='GRID'] = (outcome_df.loc[:, outcome_df.columns!='GRID'] >0)*1

# keep only outcomes of interest
outcome_df['composite'] = ((outcome_df.loc[:, outcome_df.columns.isin([phe2outcome[x] for x in keep_outcomes])]==1).any(1))*1
outcome_df['mat_composite'] = ((outcome_df.loc[:, outcome_df.columns.isin([phe2outcome[x] for x in keep_mat_outcomes])]==1).any(1))*1
outcome_df['fetal_composite'] = ((outcome_df.loc[:, outcome_df.columns.isin([phe2outcome[x] for x in keep_fetal_outcomes])]==1).any(1))*1

for_impute_df = pd.merge(outcome_df, clean_pred_df, on='GRID', how='inner')


###
###    quantify missingness
###
miss_df = for_impute_df.loc[:, for_impute_df.columns != "GRID"].isna().sum(0).reset_index().sort_values(0, ascending=False)
miss_df.columns = ['variable','count_missing']
miss_df = miss_df[miss_df['count_missing']!=0].copy()
miss_df['total_count'] = for_impute_df['GRID'].nunique()
miss_df['frac_missing'] = (miss_df['count_missing']/miss_df['total_count'])*100

# %%
# covert categorical to dummy variables
# for_impute_df= pd.get_dummies(for_impute_df, columns=['Race', 'dm_before_supervision_code'],drop_first=True)
for_impute_df= pd.get_dummies(for_impute_df, columns=['Race'],drop_first=True)
for_impute_df.set_index('GRID', inplace=True)



# %%

# sklearn impute
# import numpy as np
# from sklearn.experimental import enable_iterative_imputer
# from sklearn.impute import IterativeImputer
# imp_mean = IterativeImputer(random_state=0)
# data = for_impute_df.values
# data_tranform = imp_mean.fit_transform(data)
# imputed_df = pd.DataFrame(data_tranform, columns=for_impute_df.columns)
# imputed_df.columns


# %%
from statsmodels.imputation import mice
import statsmodels.api as sm

# for outcome in proxy_outcome_labels:
# for outcome in ['V37']:
log_df = pd.DataFrame({'outcome':[], 'status':[], 'exception':[]})


all_outcomes = dict()
for outcome in [ 'mat_composite', 'fetal_composite']:

    # outcome='composite'
    imp = mice.MICEData(for_impute_df)
    # imp = mice.MICEData(imputed_df)
    # fml = f"{outcome} ~ {' + '.join(covars)}"
    # fml = f"{outcome} ~ dm_before_supervision_code_t1dm:a1c + dm_before_supervision_code_t2dm:a1c + age_at_supervision + Race_B + Race_I + Race_U + Race_W + DEP_INDEX"
    # fml = f"{outcome} ~ dm_before_supervision_code_t1dm:a1c + dm_before_supervision_code_t2dm:a1c + age_at_supervision + Race_B + Race_I + Race_U + Race_W + DEP_INDEX"
    # fml = f"{outcome} ~ dm_before_supervision_code_t1dm + dm_before_supervision_code_t2dm + a1c + age_at_supervision + Race_B + Race_I + Race_U + Race_W + DEP_INDEX"
    fml = f"{outcome} ~  dm_outcome*a1c + age_at_supervision + Race_B + Race_I + Race_U + Race_W + DEP_INDEX"
    print(fml)
    mice_ = mice.MICE(fml, sm.Logit, imp, fit_kwds={'maxiter':2000, 'method':'newton'})


    try:
        results = mice_.fit(10, 30)
    except Exception as e:
        row = pd.DataFrame({'outcome':[outcome], 'status':['error'], 'exception':[Exception]})
        log_df = log_df.append(row)
        print("CAUGHT ERROR!")



    model = results.model
    converge_status = [model_res.mle_retvals['converged'] for model_res in model.results_list]


    # model.get_margeff(at = “overall”, method = “dyed”, atexog = c(’t1d’, ‘a1c’))


    summary_df = results.summary().tables[1]
    summary_df['models_convereged'] = np.sum(converge_status)
    summary_df['total_models'] = len(converge_status)
    summary_df['P<0.05'] = summary_df['P>|t|'] <0.05
    all_outcomes[outcome]=summary_df

# summary_df.to_csv(OUTPUT_DIR.joinpath(f'{DATE}_{outcome2phecode[outcome]}_summary.tsv'), sep="\t")


# %%
all_outcomes['mat_composite']

# %%
import pickle
pickle.dump(all_outcomes, open(OUTPUT_DIR.joinpath(f'{DATE}_outcomes_df_w_composites.pickle'),'wb'))

# %%
fig = sm.graphics.plot_partregress_grid(crime_model)
model.get_margeff(at = “overall”, method = “dyed”, atexog = c(’t1d’, ‘a1c’))