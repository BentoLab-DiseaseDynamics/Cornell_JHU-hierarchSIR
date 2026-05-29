# Cornell_JHU-hierarchSIR

A hybrid SIR - Bayesian hierarchical discrepancy model for infectious disease forecasting.

https://doi.org/10.64898/2026.05.19.26353597

## Installation (local)

Available platforms: macOS and Linux.

Note: Will not work on Windows because the OS-dependent path to the C++ Boost libraries is not included in `setup.py`.

### Setup and activate a conda environment

Update conda to make sure your version is up-to-date,

```
conda update conda
```

Setup/update the `environment`: All dependencies needed to run the scripts are collected in the conda `hierarchSIR_env.yml` file. To set up the environment,

```
conda env create -f BENTOLAB-HIERARCHSIR_conda-env.yml
conda activate HBENTOLAB-HIERARCHSIR
```

or alternatively, to update the environment (needed after adding a dependency),

```
conda activate BENTOLAB-HIERARCHSIR
conda env update -f BENTOLAB-HIERARCHSIR_conda-env.yml --prune
```

### Install the Boost libraries 

Install the C++ Boost libraries needed to integrate the multi-strain SIR model, for Linux users,

```
sudo apt-get update && sudo apt-get install -y libboost-all-dev
```

Mac users **must** install Boost through Homebrew,

```
brew install boost
```

Note: Boost is a C++ library and is not installed "inside" the conda environment but rather on your local machine. In `setup.py` the software is pointed to the location of the Boost library.

### Install the `hierarchSIR` package

Install the `hierarchSIR` Python package inside the conda environment using,

```
conda activate BENTOLAB-HIERARCHSIR
pip install -e . --force-reinstall
```

Note: The installation script requires the use of `pybind11` to "bind" the multi-strain SIR model in C++ to a Python function. This is the purpose of `pyproject.toml`.
Note: If you make any changes to the C++ files you need to reinstall `hierarchSIR`.

### Model training and forecasting

#### Training (execute once at season start)

The following procedure is performed on a state-by-state basis (for loop over states):

1. Activate the conda environment

    ```
    conda activate BENTOLAB-HIERARCHSIR 
    cd ~/scripts/operational/
    ```

2. Starting from an initial guess, optimize the hierarchSIR parameters for every training season.

    ```
    python optimize-initial_guesses.py
    ```

3. Convert the optimized model parameters per training season into an initial guess of the across-season hyperdistributions of these parameters.

    ```
    python prepare-hyperparameters.py
    ```

3. Train the model to find the across-season hyperdistributions of the parameters.

    ```
    python hierarchical_training.py
    ```

#### Forecast (performed automatically using GH actions)

1. Use the across-season hyperdistributions as priors to forecast 4-weeks ahead during the current season.
    
    ```
    python forecast.py
    ```

## Training on the JHU Rockfish cluster

See `JHU-ROCKFISH_README.md`.

## Workflows

We automatically generate and submit forecasts through the use of Github Actions. The chain of workflows is orchestrated from the master workflow file `forecast_pipeline.yml`, which triggers automatically on Wednesday at 16:00 UTC (corresponds to 12:00 EDT and 11:00 EST). It chains the following workflows together:

1. `fetch-preliminary_NHSN_HRD.yml`: Fetches a preliminary version of the NHSN HRD dataset (released Wednesday around noon), which contains the hospital admission counts per U.S. state in the epiweek ending on the previous Saturday. 

2. `backfill-preliminary_NHSN_HRD.yml`: We then corrected for systematic underreporting in these data using a generalized Dirichlet—multinomial reporting model formulated as a sequential beta-binomial survival process with coefficients estimated on a four-week rolling basis and posterior means expressed in closed form through prior conjugacy. See section 'B.7.1 Underreporting model' in https://doi.org/10.64898/2026.05.19.26353597.

3. `check-in_season.yml`: Forecasts must be submitted from mid October to the end of May. This workflow terminates the pipeline if the date is not in season. If the workflow was manually triggered by the user the pipeline continues.

4. `run-forecasts.yml`: Deploys the forecast model on a local runner (owned by @twallema) because GH servers are to slow.

5. `submit-forecasts.yml`: Submits the forecasts by committing the forecast file to the `BentoLab-DiseaseDynamics/Flusight-forecast-hub` fork of `cdcepi/Flusight-forecast-hub` (uses `SUBMISSION_PUSH_FORECAST_TO_CHALLENGE_REPO_FORK_TOKEN`; fine-grained personal access token, access on the BentoLab-DiseaseDynamics organization, access to repository `FluSight-forecast-hub`**, approval for workflows 'read metadata' and ' Read and Write access to actions, code, and pull requests', valid until May 29 2027), and then opening a PR from `BentoLab-DiseaseDynamics/Flusight-forecast-hub` to `cdcepi/Flusight-forecast-hub` (uses `SUBMISSION_PR_CREATION_TOKEN`; classic personal access token, 'repo' and 'read:org' permissions, valid until June 2027).

Post `run-forecasts.yml` a side branch is run to visualise the forecast and post it in the Bento Lab Slack (`visualise-forecasts.yml` and `post-visualisations_slack.yml`).

The workflows `fetch-preliminary_NHSN_HRD.yml`, `fetch-consolidated_NHSN_HRD.yml`, `backfill-preliminary_NHSN_HRD.yml`, `check-in_season.yml`, `run-forecasts.yml` and `visualise-forecasts.yml` use the automatically generated `GITHUB_TOKEN` with (contents: write) and (pull-requests: write) permission to create branches, commit files and open PRs in the `BentoLab-DiseaseDynamics/Cornell_JHU-hierarchSIR` repository. Tokens are created on @twallema's personal GH account (> settings > Develepor settings > Personal access tokens) and then stored under BentoLab-DiseaseDynamics/Cornell_JHU-hierarchSIR > Settings > Secrets and Variables > Actions.

** The current token from @twallema has access to the forecast repository `Cornell_JHU-hierarchSIR` and `FluSight-forecast-hub`, but only access to `FluSight-forecast-hub` is strictly needed for the submission to work.

### Workflow to-do's

- Re-generate both submission tokens and add them to the repository secretes of both `Cornell_JHU-hierarchSIR` and the new `Cornell_JHU-SCARCHhierarchSIR`.
- Add 'workflow' permission to the fine-grained personal access token `SUBMISSION_PUSH_FORECAST_TO_CHALLENGE_REPO_FORK_TOKEN` and remove access to `BentoLab-DiseaseDynamics/Cornell_JHU-hierarchSIR`.


