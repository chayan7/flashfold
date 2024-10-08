# FlashFold

A tool for protein and protein complex structure prediction.


## Installation


### Step 1: Install Conda
First, we need to make sure that we have Conda installed. It can be downloaded and installed from 
[Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Or, one can use the following commands (recommended) to install Miniforge or Mambaforge (a minimal installer for conda) 
for Linux and MacOS (Silicon and Intel) platforms:
```sh
#For Linux
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash ./Mambaforge-Linux-x86_64.sh -b

#For MacOS (Silicon)
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash ./Miniforge3-MacOSX-arm64.sh -b

#For MacOS (Intel)
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh
bash ./Mambaforge-MacOSX-x86_64.sh -b
```

### Step 2: Clone the Repository
```sh
git clone https://github.com/chayan7/flashfold.git
cd flashfold
```

### Step 3: Install dependencies under Conda Environment

#### Create conda environment: flashfold
```sh
bash install.sh
```
#### Activate the environment
```sh
conda activate flashfold
```

### Step 4: Install the Python Package
```sh
poetry install
```

[![Linux Python Package using Conda](https://github.com/chayan7/flashfold/actions/workflows/linux-python-package-conda.yml/badge.svg?event=push)](https://github.com/chayan7/flashfold/actions/workflows/linux-python-package-conda.yml)
