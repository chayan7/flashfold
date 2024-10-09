# FlashFold: a command-line tool for faster protein and protein complex structure prediction


## Introduction

Proteins are vital to cellular functions and their tertiary structure is key to understanding their biological roles. 
FlashFold predicts the structure of proteins and complexes from amino acid sequences, using AlphaFold2 models with 
a focus on speed. It also provides a table of quality metrics for the predicted structures.

- **Licence**: FlashFold is licensed under the MIT license
- **Language**: Python3 ( > 3.9 )
- **OS**: Linux, MacOS
- **OS-level Dependencies**: 
  - [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold)
  - [HMMER Suite](http://eddylab.org/software/hmmer)

## Installation

⚠️ Note:
Due to the absence of Nvidia GPU/CUDA drivers on macOS, structure predictions are **5-10 times slower on macOS** 
compared to Linux with a GPU

💡 **Tips (for Linux users):** Check prior to installation
- Ensure that the CUDA compiler driver is of version **11.8 or later** (version 12.4 is recommended). 
*You can skip this step if you do not have or plan to use a GPU*.
<pre>$ nvcc --version
nvcc: NVIDIA (R) Cuda compiler driver
Copyright (c) 2005-2022 NVIDIA Corporation
Built on Wed_Sep_21_10:33:58_PDT_2022
Cuda compilation tools, release 11.8, V11.8.89
Build cuda_11.8.r11.8/compiler.31833905_0
</pre> 
DO NOT use `nvidia-smi` to check the version.<br> See 
[NVIDIA CUDA Installation Guide for Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html) if 
you haven't installed it.

- Make sure your GNU compiler version is **9.0 or later** because `GLIBCXX_3.4.26` is required for openmm:
<pre>$ gcc --version
gcc (Ubuntu 9.3.0-17ubuntu1~20.04) 9.3.0
Copyright (C) 2019 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
</pre>
If the version is 8.5.0 or older (e.g. CentOS 7, Rocky/Almalinux 8, etc.), install a new one and add `PATH` to it.
<br>

#### 📦  FlashFold can be installed using the following methods:

<details>
<summary> 1. Python Package Index (PyPI) </summary>

FlashFold can be directly installed through [PyPI](https://pypi.org/) using the following command:
<br>

```sh
pip install -i https://test.pypi.org/simple/ flashfold==1.0.0
```
</details>

<br>
<details>
<summary> 2. Development version </summary>

>FoldFlash internally uses `LocalColabFold` (local version of ColabFold) for structure prediction. 
The installation instructions for LocalColabFold can be found [here](https://github.com/YoshitakaMo/localcolabfold). It 
is definitely useful (but **not necessary**) to install LocalColabFold before installing FlashFold.

Following steps to install:

#### Step 1: Install Conda 

Conda is a package manager that helps to install and manage dependencies. It can be downloaded and installed from:
*(Skip this step if conda is installed)*
- [Anaconda](https://www.anaconda.com/products/distribution) 
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
- [Conda-forge](https://conda-forge.org/) (recommended)

#### Step 2: Clone the git repository
```sh
git clone https://github.com/chayan7/flashfold.git
cd flashfold
```

#### Step 3: Install dependencies under conda environment

```sh
bash install.sh              # Install dependencies
conda activate flashfold     # Activate the environment
```

#### Step 4: Install the Python Package
```sh
poetry install
```
</details>



## Usage


## Acknowledgements

FlashFold utilizes and/or references the following separate libraries
and packages:

*   [Abseil](https://github.com/abseil/abseil-py)
*   [Alphafold](https://github.com/google-deepmind/alphafold)
*   [Colabfold](https://github.com/sokrypton/ColabFold)
*   [HH Suite](https://github.com/soedinglab/hh-suite)*
*   [MMseqs2](https://github.com/soedinglab/mmseqs2)*
*   [HMMER Suite](http://eddylab.org/software/hmmer)
*   [Biopython](https://biopython.org)
*   [Chex](https://github.com/deepmind/chex)
*   [Haiku](https://github.com/deepmind/dm-haiku)
*   [Immutabledict](https://github.com/corenting/immutabledict)
*   [JAX](https://github.com/google/jax/)
*   [Kalign](https://msa.sbc.su.se/cgi-bin/msa.cgi)
*   [matplotlib](https://matplotlib.org/)
*   [ML Collections](https://github.com/google/ml_collections)
*   [NumPy](https://numpy.org)
*   [OpenMM](https://github.com/openmm/openmm)
*   [OpenStructure](https://openstructure.org)
*   [pandas](https://pandas.pydata.org/)
*   [pymol3d](https://github.com/avirshup/py3dmol)
*   [SciPy](https://scipy.org)
*   [Sonnet](https://github.com/deepmind/sonnet)
*   [TensorFlow](https://github.com/tensorflow/tensorflow)
*   [tqdm](https://github.com/tqdm/tqdm)
*   [ujson](https://github.com/ultrajson/ultrajson)


[![Linux Python Package using Conda](https://github.com/chayan7/flashfold/actions/workflows/linux-python-package-conda.yml/badge.svg?event=push)](https://github.com/chayan7/flashfold/actions/workflows/linux-python-package-conda.yml)