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
from scipy.stats import beta as beta_dist

# Define all paths reletive to this file
abs_dir = os.path.dirname(__file__)

########################################
## Retrieve preliminary NHSN HRD data ##
########################################

# Find all preliminary .parquet files and read them into a list
parquet_files = sorted(glob.glob(f"{abs_dir}/*.gzip"))
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
        d = df.loc[df['date'] == focal_date][['date', 'name_state', 'influenza admissions']]
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

###############################################################################
## Fit a beta-binomial model to estimate average missigness of reported data ##
###############################################################################

# For each state i and week X we observe:
#
#   y_{i,X}^{(0)}  = number of influenza admissions reported in week X
#   y_{i,X}^{(2)}  = number of influenza admissions reported in week X+2
#
# where y_{i,X}^{(2)} is treated as the "final" count after reporting delays.
# We assume monotone backfill: early reports are a thinning of the final count.
#
#   y_{i,X}^{(0)} | y_{i,X}^{(2)}, p_i  ~  Binomial(y_{i,X}^{(2)}, p_i)
#
# To stabilize estimation in states with sparse data, we use a Beta prior on p:
#
#   p_i ~ Beta(alpha, beta)
#
# Summing counts of state i over the available weeks, the posterior distribution can be computed analytically:
#
#   p_i | data ~ Beta(
#       alpha + sum_y0_i,
#       beta  + (sum_y2_i - sum_y0_i)
#   )
# Which for alpha = beta = 1 and sum_y0_i >> 1, sum_y2_i >> 1 is simply the ratio of the sum of the backfilled to non-filled cases over the horizon. 

# Beta distribution priors
alpha_02 = 20
beta_02 = 1    # E[p_i] = 0.95 (5% missing on release)
alpha_12 = 50
beta_12 = 1    # E[p_i] = 0.95 (5% missing on release)

# Aggregate dounts per state
sum_df = pd.concat(abs_backfill_collect, ignore_index=True)

# Aggregate evidence per state
posterior = (
    sum_df
    .groupby('name_state', as_index=False)
    .agg(
        y0_sum=("influenza_admissions_0", "sum"),
        y1_sum=("influenza_admissions_1", "sum"),
        y2_sum=("influenza_admissions_2", "sum")
    )
)

# Compute posterior moments of p (0 --> 2)
alpha_post = alpha_02 + posterior["y0_sum"]
beta_post = beta_02 + (posterior["y2_sum"] - posterior["y0_sum"])
posterior["p_02_mean"] = alpha_post / (alpha_post + beta_post)
posterior["p_02_median"] = beta_dist.ppf(
    0.5, alpha_post, beta_post
)
posterior["p_02_low_90"] = beta_dist.ppf(
    0.05, alpha_post, beta_post
)
posterior["p_02_high_90"] = beta_dist.ppf(
    0.95, alpha_post, beta_post
)

# Compute posterior moments of p (1 --> 2)
alpha_post = alpha_12 + posterior["y1_sum"]
beta_post = beta_12 + (posterior["y2_sum"] - posterior["y1_sum"])
posterior["p_12_mean"] = alpha_post / (alpha_post + beta_post)
posterior["p_12_median"] = beta_dist.ppf(
    0.5, alpha_post, beta_post
)
posterior["p_12_low_90"] = beta_dist.ppf(
    0.05, alpha_post, beta_post
)
posterior["p_12_high_90"] = beta_dist.ppf(
    0.95, alpha_post, beta_post
)

# Cap all cases where p > 1 to p = 1 (assumption that all retracting cases represents errors)
posterior = posterior.fillna(1)

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

# deliberately chose not to round values --> we use a poisson likelihood function to fit the model which is continuous

# Save backfill amounts and data
parquet_filenames = [os.path.basename(f) for f in parquet_files]
posterior.to_csv(os.path.join(abs_dir, '../preliminary_backfilled/'+parquet_filenames[-1][:-13]+'_backfill_beta-binomial-estimates.csv'))
latest_df.to_parquet(os.path.join(abs_dir, '../preliminary_backfilled/'+parquet_filenames[-1]), compression='gzip', index=False)

