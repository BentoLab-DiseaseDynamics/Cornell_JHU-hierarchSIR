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
        data.append(df.loc[df['date'] == focal_date][['date', 'name_state', 'influenza admissions']])

    # differentiate and send to one dataframe
    ## pre-allocate output
    abs_backfill = data[0].copy()
    ## change column name
    abs_backfill.rename()
    ## fill
    for k in range(len(data)-1):
        abs_backfill[f'influenza_admissions_{k+1}'] = data[k+1]['influenza admissions']
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
alpha = 20
beta = 1    # E[p_i] = 0.95 (5% missing on release)

# Aggregate dounts per state
sum_df = pd.concat(dfs, ignore_index=True)

# Aggregate evidence per state
posterior = (
    sum_df
    .groupby(state_col, as_index=False)
    .agg(
        y0_sum=("influenza admissions", "sum"),
        y2_sum=("influenza_admissions_2", "sum")
    )
)

# Compute posterior moments of p
alpha_post = alpha + posterior["y0_sum"]
beta_post = beta + (posterior["y2_sum"] - posterior["y0_sum"])
posterior["p_mean"] = alpha_post / (alpha_post + beta_post)
posterior["p_median"] = beta_dist.ppf(
    0.5, alpha_post, beta_post
)
posterior["p_low_90"] = beta_dist.ppf(
    0.05, alpha_post, beta_post
)
posterior["p_high_90"] = beta_dist.ppf(
    0.95, alpha_post, beta_post
)

print(posterior.head(50))
