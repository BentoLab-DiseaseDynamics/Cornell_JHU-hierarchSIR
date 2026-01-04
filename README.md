# Cornell_JHU-hierarchSIR

A hybrid SIR - Bayesian hierarchical discrepancy model for infectious disease forecasting.

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

- Automated forecast submission requires the use of two tokens PR_CREATION_TOKEN and SUBMIT_FORECAST_TOKEN, saved in the repository secrets (GH repo page > Settings > Secrets and Variables > Actions). Both tokens were generated from @twallema's personal GH account. PR_CREATION_TOKEN is a classic personal access token with 'repo' and 'read:org' permissions and is valid until the end of the challenge (May 31, 2026). This token is used to open a PR from BentoLab-DiseaseDynamics/Flusight-forecast-hub to cdcepi/Flusight-forecast-hub in @twallema's name. SUBMIT_FORECAST_TOKEN is a fine-grained access token, which grants @twallema access to the BentoLab-DiseaseDynamics/Cornell_JHU-hierarchSIR repository with approval for workflows 'read metadata' and ' Read and Write access to actions, code, and pull requests'. The token is valid until Nov 2026. Both tokens will have to be renewed during the summer of 2026 for the 2026-2027 challenge.
