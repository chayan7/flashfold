<p align="center">
  <img src="https://github.com/chayan7/flashfold/blob/main/logo/flashfold.png" alt="Logo" width="100">
</p>

# FlashFold: a command-line tool for faster protein and protein complex structure prediction


[![Testing on Linux](https://github.com/chayan7/flashfold/actions/workflows/linux-python-package-conda.yml/badge.svg?event=push)](https://github.com/chayan7/flashfold/actions/workflows/linux-python-package-conda.yml)

## Introduction

Proteins are vital to cellular functions and their tertiary structure is key to understanding their biological roles. 
FlashFold predicts the structure of proteins and complexes from amino acid sequences, using AlphaFold2 models with 
a focus on speed. It also provides a table of quality metrics for the predicted structures.

- License: FlashFold is licensed under the MIT license
- Language: Python3 ( > 3.9 )
- OS: Linux, macOS
- OS-level Dependencies: 
  - [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold)
  - [HMMER Suite](http://eddylab.org/software/hmmer)

## Installation

FlashFold can be installed on Linux and macOS. 

###### 🚨 *Important: If you are using macOS, please note that the structure prediction is 5-10 times slower compared to Linux with a GPU. <br> This is due to the absence of Nvidia GPU/CUDA drivers on macOS.*

If you are planning to use a GPU, it is recommended to check the following settings prior to installation:

<details> <summary> CUDA 12.1 or later (version 12.4 is recommended) and cudnn 9 are required. 
(<i>If you are planning to use a GPU)</i></summary>

- You can check the CUDA version using the following command: 

  ```sh
  nvcc --version
  ```

- DO N🚫T use `nvidia-smi` to check the version. ❌ <br> ✔️ See 
[NVIDIA CUDA Installation Guide for Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html) if 
you haven't installed it.

</details> 
<br>

<details> <summary> GNU compiler version is **9.0 or later** is required. </summary>

- You can check the GNU compiler version using the following command:
  ```sh
  gcc --version
  ```
  💡 If the version is 8.5.0 or older (e.g. CentOS 7, Rocky/Almalinux 8, etc.), install a new one and add `PATH` to it.

</details>
<br>

##### 📌 FlashFold can be installed using the following steps:

###### ✔ Step 1: Install Conda (*Skip this step if conda is already installed*)

Conda is a package manager that helps to install and manage dependencies. It can be downloaded and installed from:

- [Anaconda](https://www.anaconda.com/products/distribution) 
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
- [Conda-forge](https://conda-forge.org/) (recommended)

###### ✔ Step 2: Clone the git repository
  ```sh
  git clone https://github.com/chayan7/flashfold.git
  cd flashfold
  ```

###### ✔ Step 3: Install dependencies under conda environment

FlashFold internally uses `LocalColabFold` (local version of ColabFold) for structure prediction. 
The installation instructions for LocalColabFold can be found [here](https://github.com/YoshitakaMo/localcolabfold). 

To streamline the installation process for both Linux and macOS users, FlashFold provides a convenient installation 
script that sets up the required dependencies within a conda environment named `flashfold`.

  ```sh
  bash install.sh              # Install dependencies
  conda activate flashfold     # Activate the environment
  ```

###### ✔ Step 4: Install the package
  ```sh
  poetry install
  ```

###### ✔ Step 5: Run the tests
  ```sh
  poetry run pytest
  ```
  or, 
  ```sh
  pytest
  ```

## Workflow

FlashFold uses amino acid sequences to predict the structure of proteins and protein complexes. In order to achieve this,
FlashFold uses the following steps:

1. *Sequence Alignment*: FlashFold uses `jackhmmer` to generate a multiple sequence alignment (MSA) for the input 
sequence. FlashFold reduces the MSA generation time significantly by using a compact database.
2. *Structure Prediction*: The MSA is then formatted and used as an input for `colabfold_batch` to predict the structure.
3. *Model Refinement (optional)*: Based on user input, the predicted structure is refined using `OpenMM` and `OpenStructure`.
4. *Quality Metrics*: FlashFold provides a table of quality metrics for the predicted structures. For protein 
complexes, it uses the Predicted DockQ version 2 ([pDockQ2](https://doi.org/10.1093/bioinformatics/btad424)) script to 
calculate the quality of each interface.

## Application

- ###  Database
  In order to predict the structure of proteins and protein complexes, FlashFold requires a sequence database. The database
  is used for homology sequence detection as the input sequence to generate a multiple sequence alignment (MSA) . 
  FlashFold provides the following options:

  <br>
  <details><summary>Download in-built database</summary>

  FlashFold provides three in-built databases, that can be downloaded using the following command:
  ```sh
  flashfold download_db -i /path/to/database.json -o /path/to/downloaded_db/
  ```
  The `database.json` file can be found [here](https://github.com/chayan7/flashfold/blob/main/database.json). 
  User can avoid downloading a database by removing the database name and the download link in the json file.

  </details>
  <br>
  <details><summary>Create custom database</summary>

  FlashFold allows user to create custom database using the `create_db` subcommand. In this case, the input should be the 
  assembled genome data in [GenBank](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gbff/) format. 
  <br>In order to download the genome data from NCBI, FlashFold provides a convenient script `ncbi_data` that can be used as follows:<br><br>

    - For example, to download all the genbank files of <i>[Pseudomonas aeruginosa](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=287&annotated_only=true&refseq_annotation=true)</i> 
    form NCBI RefSeq, the following command can be used:
  
      ```sh
      flashfold ncbi_data  -n "Pseudomonas aeruginosa" -f gbff -s refseq -o /path/to/genbank_file_dir/ 
      ```

    - Or, user can download the genbank files of particular genome of interest from NCBI using accession numbers as input, see
    [example](https://github.com/chayan7/flashfold/blob/main/test/input/others/assembly_accessions.txt).

      ```sh
      flashfold ncbi_data  -i /path/to/assembly_accessions.txt -f gbff -o /path/to/genbank_file_dir/
      ```
  Once the genbank files are downloaded, the custom database can be created using the `create_db` subcommand as follows:
    ```sh
    flashfold create_db -p /path/to/genbank_file_dir/ -o /path/to/custom_db/
    ```
  </details>
  <br>
  <details><summary>Extend database</summary>
  
  FlashFold allows user to update or extend the current database with new information.
  - If user would like to extend or update database_1 with the information from database_2, it is possible by 
  using the `extend_db` subcommand.

    ```sh
    flashfold extend_db -m /path/to/database_1 -n /path/to/database_2
    ```
     Note that, only the database_1 will be updated with the new information from database_2.

  - It is also possible to extend the current database directly with the new collection of genbank files, using 
  the `extend_db` subcommand. 
  
    ```sh
    flashfold extend_db -m /path/to/database_to_be_extended -g /path/to/genbank_file_dir/
    ```
  </details>

- ###  Protein structure prediction
    FlashFold provides a subcommand `fold` to predict the structure of proteins and protein complexes. See details below:

  <br>
  <details><summary>Input file preparation</summary>

  - FlashFold takes amino acid sequence in FASTA format as input. Also, it can take multiple FASTA files as input 
  when `--batch` is set. The input file should follow the following guidelines:
    - It is recommended to keep the file name short and readable. Avoid using special characters in the file name.
    - It should be noted that, when `--batch` is set, the file name will be used as a directory to store results 
    under user provided output directory. If any special characters are found `except "." or "_"` in the file name, 
    it will be replaced with `"_"`.
    - File extension should be `.fasta`.

  - Additionally, FlashFold can take A3M file as input. The A3M file preferably should be generated using `FlashFold` 
  itself using the `--only_msa` option. User customised A3M file can be served as input as well. `--batch` option 
  is also applicable for A3M file input as FASTA.

  Few examples for FASTA sequence as input are shown below:

  **Monomer**
  ```
  >seq_1
  FHWDREGQADDSSSCWLRVASGWAGRNYGAIAIPRVGMEVLVTFLEGDPDQPLVTGCLFH
  REHPVPYELPGHKTRSVFKSLSSPGGGGYNELRIEDRKGQEQIFVHAQR
  ```
  
  **Protein complex**
  - Homo-dimer
    ```
    >seq_1
    FHWDREGQADDSSSCWLRVASGWAGRNYGAIAIPRVGMEVLVTFLEGDPDQPLVTGCLFH
    REHPVPYELPGHKTRSVFKSLSSPGGGGYNELRIEDRKGQEQIFVHAQR
    >seq_1
    FHWDREGQADDSSSCWLRVASGWAGRNYGAIAIPRVGMEVLVTFLEGDPDQPLVTGCLFH
    REHPVPYELPGHKTRSVFKSLSSPGGGGYNELRIEDRKGQEQIFVHAQR
    ```
  - Hetero-dimer
    ```
    >seq_1
    FHWDREGQADDSSSCWLRVASGWAGRNYGAIAIPRVGMEVLVTFLEGDPDQPLVTGCLFH
    REHPVPYELPGHKTRSVFKSLSSPGGGGYNELRIEDRKGQEQIFVHAQR
    >seq_2
    MTSWTLVTLVLLIILAAIRPEQLQVVAYKLVLVTLGAVAGYWIDRSLFPYVARPHECSAN
    LVVVGAWLRRGLIVLACILGLTLGL
    ```
  - Hetero-trimer
    ```
    >seq_1
    FHWDREGQADDSSSCWLRVASGWAGRNYGAIAIPRVGMEVLVTFLEGDPDQPLVTGCLFH
    REHPVPYELPGHKTRSVFKSLSSPGGGGYNELRIEDRKGQEQIFVHAQR
    >seq_2
    MTSWTLVTLVLLIILAAIRPEQLQVVAYKLVLVTLGAVAGYWIDRSLFPYVARPHECSAN
    LVVVGAWLRRGLIVLACILGLTLGL
    >seq_3
    MAFQADRFLWFNSSSGQTVAPVSIVGGQMFINTAMIQDGSITNAKIGNVIQSTALGANGE
    PLWKLDKAGSLTMNSATSGGFMRQTAEAVKVYDANLVLRVQIGNLDA
    ```
    </details>
    <br>
    <details><summary>Commands</summary>
    
    FlashFold offers subcommand `fold` to predict the structure of proteins and protein complexes. FlashFold uses 
    different algorithm and model for monomer and multimer prediction. However, the user does not need to specify it 
    because FlashFold can automatically detect based on the input sequence. 
    
    Few examples are shown below:
    
    __Beginner__
    
    ```shell
    flashfold fold -q /path/to/query.fasta -d /path/to/database/ -o /path/to/output/ -t number_of_threads
    ```
    __Moderate__
    
    ```shell
    flashfold fold -q /path/to/query.fasta -d /path/to/database/ -o /path/to/output/ -t number_of_threads --only_msa   
    ```
    
    __Advanced__
    
    ```shell
    flashfold fold -q /path/to/query.a3m -o /path/to/output/  
    ```
    
    __Expert__
    
    ```shell
    flashfold fold -q /path/to/dir/multiple_fasta_files --batch -d /path/to/database/ -o /path/to/output/ -t number_of_threads
    ```
    </details>
  

## Acknowledgements

FlashFold utilizes and/or references the following separate libraries
and packages:

*   [Abseil](https://github.com/abseil/abseil-py)
*   [Alphafold](https://github.com/google-deepmind/alphafold)
*   [Colabfold](https://github.com/sokrypton/ColabFold)
*   [HMMER Suite](http://eddylab.org/software/hmmer)
*   [Biopython](https://biopython.org)
*   [Chex](https://github.com/deepmind/chex)
*   [Haiku](https://github.com/deepmind/dm-haiku)
*   [Immutabledict](https://github.com/corenting/immutabledict)
*   [JAX](https://github.com/google/jax/)
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


