#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2021-10-15 21:32:24

import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

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
# MAIN
# -----------