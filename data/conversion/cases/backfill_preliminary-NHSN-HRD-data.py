"""
This script loads the archived preliminary NHSN HRD data and estimates the relative change two weeks post-reporting due to backfill.
It then takes the latest available preliminary dataset in `~/data/interim/cases/NHSN-HRD_archive/preliminary`,
backfills it and saves it in the `~/data/interim/cases/NHSN-HRD_archive/preliminary_backfilled` folder.
"""

__author__      = "T.W. Alleman"
__copyright__   = "Copyright (c) 2025 by T.W. Alleman, IDD Group (JHUBSPH) & Bento Lab (Cornell CVM). All Rights Reserved."

##################
## Dependencies ##
##################

import os
import glob
import numpy as np
import pandas as pd

# Define all paths reletive to this file
abs_dir = os.path.dirname(__file__)


########################################
## Retrieve preliminary NHSN HRD data ##
########################################

# Length of rolling backfill window
N = 4

# Find all preliminary .parquet files and read them into a list
parquet_files = sorted(glob.glob(os.path.join(abs_dir, "../../interim/cases/NHSN-HRD_archive/preliminary/*.gzip")))
dfs = []
for file in parquet_files:
    df = pd.read_parquet(file)
    dfs.append(df)

# Loop over the files in threes
abs_backfill_collect = []
for i in range(len(dfs) - 3 + 1):
    # get data from files 0, 1, 2, then 1, 2, 3, then 2, 3, 4 etc.
    data = []
    for j in range(3):
        # get file
        df = dfs[i+j]
        # if first file in sequence of three --> save last date
        if j == 0:
            focal_date = max(df['date'])
        # extract data
        d = df.loc[df['date'] == focal_date][['date', 'name_state', 'fips_state', 'influenza admissions']]
        # rename data column
        d = d.rename(columns={'influenza admissions': 'influenza_admissions_0'})
        data.append(d)

    # differentiate and send to one dataframe
    ## pre-allocate output
    abs_backfill = data[0].copy()
    ## fill
    for k in range(len(data)-1):
        abs_backfill[f'influenza_admissions_{k+1}'] = data[k+1]['influenza_admissions_0']
    abs_backfill_collect.append(abs_backfill)

######################################################
## Generalized Dirichlet–Multinomial Backfill model ##
######################################################

# Based on the work of https://academic.oup.com/biometrics/article/76/3/789/7429141#biom13188-sec-0050

# Generalized Dirichlet priors (sequential hazards)
a0_prior, b0_prior = 45, 5     # immediate reporting in week 0 (E[X] = 0.9)
a1_prior, b1_prior = 40, 10    # fraction of not immediately reported in week 0, reported in week 1 (E[X] = 0.8 --> 2% remaining after 1 week)

# Aggregate counts per state (same as before)
sum_df = pd.concat(abs_backfill_collect[-N:], ignore_index=True)

posterior = (
    sum_df
    .groupby('name_state', as_index=False)
    .agg(
        y0_sum=("influenza_admissions_0", "sum"),
        y1_sum=("influenza_admissions_1", "sum"),
        y2_sum=("influenza_admissions_2", "sum")
    )
)

# Convert cumulative counts into increments
posterior["z0"] = posterior["y0_sum"]
posterior["z1"] = posterior["y1_sum"] - posterior["y0_sum"]
posterior["z2"] = posterior["y2_sum"] - posterior["y1_sum"]

# Final totals
posterior["n"] = posterior["z0"] + posterior["z1"] + posterior["z2"]

# Posterior updates (conjugate)
posterior["a0_post"] = a0_prior + posterior["z0"]
posterior["b0_post"] = b0_prior + (posterior["n"] - posterior["z0"])
posterior["a1_post"] = a1_prior + posterior["z1"]
posterior["b1_post"] = b1_prior + (posterior["n"] - posterior["z0"] - posterior["z1"])

# Posterior means of hazards
posterior["v0_mean"] = posterior["a0_post"] / (posterior["a0_post"] + posterior["b0_post"])
posterior["v1_mean"] = posterior["a1_post"] / (posterior["a1_post"] + posterior["b1_post"])

# Completeness fractions (analytic)
posterior["p_02_mean"] = np.clip(np.round(posterior["v0_mean"], 3), 0, 1)
posterior["p_12_mean"] = np.clip(np.round(1.0 - (1.0 - posterior["v0_mean"]) * (1.0 - posterior["v1_mean"]), 3), 0, 1)


###############################################
## Backfill latest preliminary NHSN HRD data ##
###############################################

# Get the latest dataset and date
latest_df = dfs[-1]
latest_date = max(latest_df['date'])

# backfill the most recent week --> shoot forward to two weeks of backfilling total
latest_df = latest_df.merge(posterior[['name_state', 'p_02_mean', 'p_12_mean']], on='name_state')
latest_df.loc[latest_df['date'] == latest_date, 'influenza admissions'] *= 1/latest_df.loc[latest_df['date'] == latest_date, 'p_02_mean'].values 

# backfill last week's data --> shoot forward to two weeks of backfilling total
latest_minus1_date = sorted(latest_df['date'].unique().tolist())[-2]
latest_df.loc[latest_df['date'] == latest_minus1_date, 'influenza admissions'] *= 1/latest_df.loc[latest_df['date'] == latest_minus1_date, 'p_12_mean'].values 

# remove the p_02_mean and p_12_mean columns
latest_df = latest_df.drop(columns=['p_02_mean', 'p_12_mean'])

# put fips_state back in and sort
posterior = posterior.merge(abs_backfill_collect[-1][['fips_state', 'name_state']], on='name_state')
posterior = posterior.sort_values(by='fips_state')

# add in first and last date used in sliding window
posterior['start_backfill_window'] = sum_df['date'].unique()[0]
posterior['end_backfill_window'] = sum_df['date'].unique()[-1]

# slice out relevant columns
posterior = posterior[['fips_state', 'name_state', 'start_backfill_window', 'end_backfill_window', 'p_02_mean', 'p_12_mean']]

# Save beta-binomial estimates and the data
parquet_filenames = [os.path.basename(f) for f in parquet_files]
posterior.to_csv(os.path.join(abs_dir, '../../interim/cases/NHSN-HRD_archive/preliminary_backfilled/'+parquet_filenames[-1][:-13]+'_backfill_beta-binomial-estimates.csv'), index=False)
latest_df.to_parquet(os.path.join(abs_dir, '../../interim/cases/NHSN-HRD_archive/preliminary_backfilled/'+parquet_filenames[-1]), compression='gzip', index=False)

