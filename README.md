# FlashFold

A tool for protein and protein complex structure prediction.

# FlashFold

## Installation

### Step 1: Install Conda
Make sure you have Conda installed. You can download and install it from [Anaconda](https://www.anaconda.com/products/distribution) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

### Step 2: Clone the Repository
```sh
git clone https://github.com/chayan7/flashfold.git
cd flashfold
```

### Step 3: Create and Activate the Conda Environment
```sh
conda env create -f environment.yml
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