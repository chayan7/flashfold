# FlashFold

A tool for protein and protein complex structure prediction.


## Installation

[![Python Package using Conda](https://github.com/chayan7/flashfold/actions/workflows/python-package-conda.yml/badge.svg?branch=main&event=push)](https://github.com/chayan7/flashfold/actions/workflows/python-package-conda.yml)

### Step 1: Install Conda
Make sure you have Conda installed. You can download and install it from [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

### Step 2: Clone the Repository
```sh
git clone https://github.com/chayan7/flashfold.git
cd flashfold
```

### Step 3: Create and Activate the Conda Environment
#### Create the environment
```sh
conda env update -f envs/linux-environment.yml # For Linux

```
or
```sh
conda env update -f envs/mac-environment.yml # For MacOS
```
#### Activate the environment
```sh
conda activate flashfold
```

### Step 4: Install the Python Package
```sh
pip install .
```

### Step 5: Verify the installation
```sh
python -c "import flashfold; print('FlashFold installed successfully!')"
```
