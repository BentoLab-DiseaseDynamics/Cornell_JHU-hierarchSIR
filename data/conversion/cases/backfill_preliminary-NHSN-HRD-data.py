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

# Length of rolling backfill window
N = 6

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


#################################################################
# Dirichlet–Multinomial Backfill Model (Mathematical Description)
#################################################################
#
# Goal
# ----
# Estimate reporting completeness of influenza hospital admissions in the
# most recent two preliminary NHSN HRD releases by modeling reporting delays.
#
# Data structure
# --------------
# For a given epidemiological week X, let:
#
#   y_{X,X}     = number of admissions in week X reported in week X (on release)
#   y_{X,X+1}   = number of admissions in week X reported in week X+1 (after 1 week of backfill)
#   y_{X,X+2}   = number of admissions in week X reported in week X+2 (after 2 weeks of backfill)
#
# We assume that by week X+2 the data are effectively final. We define reporting increments:
#
#   z_0 = y_{X,X}
#   z_1 = y_{X,X+1} - y_{X,X}
#   z_2 = y_{X,X+2} - y_{X,X+1}
#   z_0 + z_1 + z_2 = y_{X,X+2}
#
# Model
# -----
# Conditional on the reporting-delay proportions \vec{pi} = (pi_0, pi_1, pi_2),
# the incremental completion of the data follow a multinomial allocation of the final count,
#
#   (z_0, z_1, z_2) | π
#       ~ Multinomial(
#             n = y_{X,X+2},
#             p = (pi_0, pi_1, pi_2)
#         )
#
# where:
#   pi_0 = fraction of cases reported immediately (week 0),
#   pi_1 = fraction reported with a 1-week delay,
#   pi_2 = fraction reported with a 2-week delay.
#
# Prior
# -----
# We place a Dirichlet prior on π:
#
#   (π_0, π_1, π_2) ~ Dirichlet(α_0, α_1, α_2)
#
# with parameters chosen to encode prior beliefs:
#
#   E[π_0]           ≈ 0.90   (≈90% completeness at release)
#   E[π_0 + π_1]     ≈ 0.98   (≈98% completeness after 1 week)
#
# Achieved by setting:
#
#   α_k = kappa * m_k
#
# where:
#   m = (0.95, 0.04, 0.01) is the prior mean vector,
#   kappa > 0 controls prior strength (effective sample size).
#
# Posterior
# ---------
# By conjugacy, given aggregated increments z = (z_0, z_1, z_2),
#
#   (π_0, π_1, π_2) | z
#       ~ Dirichlet(
#             α_0 + z_0,
#             α_1 + z_1,
#             α_2 + z_2
#         )
#
# Posterior means
# ---------------
# The posterior mean of each component is:
#
#   E[π_k | z]
#       = (α_k + z_k) / (α_0 + α_1 + α_2 + z_0 + z_1 + z_2),
#       for k = 0, 1, 2.
#
# Reporting completeness factors
# ------------------------------
# The quantities used for backfilling are cumulative completeness fractions:
#
#   p_02 = π_0
#       = fraction of final cases observed at initial release,
#
#   p_12 = π_0 + π_1
#       = fraction of final cases observed after 1 week.
#
# Their posterior means are therefore:
#
#   E[p_02 | z] = E[π_0 | z]
#
#   E[p_12 | z] = E[π_0 + π_1 | z]
#               = E[π_0 | z] + E[π_1 | z]
#
# Backfill correction
# -------------------
# Let y_{X,X} and y_{X-1,X} denote the most recent and second-most-recent
# preliminary observations available at release week X.
#
# These are backfilled to approximate the final (week X+2) counts by:
#
#   y*_{X,X}     = y_{X,X}     / E[p_02 | z]
#   y*_{X-1,X}   = y_{X-1,X}   / E[p_12 | z]
#
# The model guarantees:
#   0 ≤ p_02 ≤ p_12 ≤ 1,
####################################################################

# Dirichlet priors
kappa = 50  # regularization strength 
alpha0 = 0.9 * kappa    # immediate
alpha1 = 0.08 * kappa   # week +1
alpha2 = 0.02 * kappa   # week +2

# Aggregate dounts per state
sum_df = pd.concat(abs_backfill_collect[-N:], ignore_index=True)

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

# Convert cumulative counts into increments
posterior["z0"] = posterior["y0_sum"]
posterior["z1"] = posterior["y1_sum"] - posterior["y0_sum"]
posterior["z2"] = posterior["y2_sum"] - posterior["y1_sum"]

# Total final counts
posterior["z_sum"] = posterior["z0"] + posterior["z1"] + posterior["z2"] # equal to y2_sum

# Posterior means of Dirichlet components
denom = alpha0 + alpha1 + alpha2 + posterior["z_sum"]
posterior["pi0_mean"] = np.round((alpha0 + posterior["z0"]) / denom, 2)
posterior["pi1_mean"] = np.round((alpha1 + posterior["z1"]) / denom, 2)
posterior["pi2_mean"] = np.round((alpha2 + posterior["z2"]) / denom, 2)

# Backfill completeness factors (drop-in replacements)
posterior["p_02_mean"] = np.round(posterior["pi0_mean"], 2)
posterior["p_12_mean"] = np.round(posterior["pi0_mean"] + posterior["pi1_mean"], 2)

# Cap all cases where p > 1 to p = 1 (assumption that all retraction of cases represents an error)
posterior["p_02_mean"] = posterior["p_02_mean"].clip(upper=1)
posterior["p_12_mean"] = posterior["p_12_mean"].clip(upper=1)


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

