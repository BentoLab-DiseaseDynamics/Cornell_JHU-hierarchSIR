"""
This script downloads, formats and archives the NHSN HRD dataset
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
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# Define relevant global  variables
abs_dir = os.path.dirname(__file__)

###################
## Retrieve data ##
###################

# Find all .parquet files
parquet_files = sorted(glob.glob(f"{abs_dir}/*.gzip"))

# Read all files into a list of DataFrames
dfs = []
for file in parquet_files:
    df = pd.read_parquet(file)
    dfs.append(df)

# Loop over the files in threes
rel_backfill_collect = [] 
for i in range(len(dfs) - 3 + 1):
    # get data from files 0, 1, 2, then 1, 2, 3, etc.
    data = []
    for j in range(3):
        # get file
        df = dfs[i+j]
        # if first file --> save last date
        if j == 0:
            focal_date = max(df['date'])
        # extract data
        data.append(df.loc[df['date'] == focal_date][['date', 'name_state', 'influenza admissions']])

    # differentiate and send to one dataframe
    ## pre-allocate output
    rel_backfill = data[0]
    ## fill
    for k in range(len(data)-1):
        rel_backfill[f'rel_adj_week_{k+1}'] = np.round((data[k+1]['influenza admissions'] - data[0]['influenza admissions']) / data[0]['influenza admissions'] * 100, 1)
    rel_backfill_collect.append(rel_backfill)