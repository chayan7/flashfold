# FlashFold

A tool for protein and protein complex structure prediction.


## Installation


### Step 1: Install Conda
First, we need to make sure that we have Conda installed. It can be downloaded and installed from [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Or, one can use the following command:
```sh

#For Linux
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash ./Mambaforge-Linux-x86_64.sh

#For MacOS
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash ./Miniforge3-MacOSX-arm64.sh

```


### Step 2: Clone the Repository
```sh
git clone https://github.com/chayan7/flashfold.git
cd flashfold
```

### Step 3: Create and Activate the Conda Environment

#### Linux
```sh
conda env update -f envs/linux-environment.yml
```
#### MacOS
```sh
conda env update -f envs/mac-environment.yml
```
#### Activate the environment
```sh
conda activate flashfold
```

### Step 4: Install the Python Package
```sh
pip install -e .
```

[![Linux Python Package using Conda](https://github.com/chayan7/flashfold/actions/workflows/linux-python-package-conda.yml/badge.svg?event=push)](https://github.com/chayan7/flashfold/actions/workflows/linux-python-package-conda.yml)
