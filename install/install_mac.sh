#!/bin/bash -e

# check commands
type wget 2>/dev/null || { echo -e "Please install wget using Homebrew:\n\tbrew install wget" ; exit 1 ; }
type hhsearch 2>/dev/null || { echo -e "Please install hh-suite using Homebrew:\n\tbrew install brewsci/bio/hh-suite" ; exit 1 ; }
type kalign 2>/dev/null || { echo -e "Please install kalign using Homebrew:\n\tbrew install kalign" ; exit 1 ; }
type mmseqs 2>/dev/null || { echo -e "Please install mmseqs2 using Homebrew:\n\tbrew install mmseqs2" ; exit 1 ; }
type jackhmmer 2>/dev/null || { echo -e "Please install jackhmmer using Homebrew:\n\tbrew install hmmer" ; exit 1 ; }


CURRENTPATH="$(pwd)"
ENV_YML_FILE_PATH="${CURRENTPATH}/envs/mac-environment.yml"
FLASHFOLD_CONDA_ENV_DIR="${CURRENTPATH}/src/foldflash/.env"
FLASHFOLD_ENV_NAME="flashfold"

mkdir -p "${FLASHFOLD_CONDA_ENV_DIR}"
cd "${FLASHFOLD_CONDA_ENV_DIR}"
wget -q -P . https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash ./Miniforge3-MacOSX-arm64.sh -b -u -p "${FLASHFOLD_CONDA_ENV_DIR}/conda"
rm ${FLASHFOLD_CONDA_ENV_DIR}/Miniforge3-MacOSX-arm64.sh

source "${FLASHFOLD_CONDA_ENV_DIR}/conda/etc/profile.d/conda.sh"
export PATH="${FLASHFOLD_CONDA_ENV_DIR}/conda/condabin:${PATH}"
conda update -n base conda -y
conda env create -f "${ENV_YML_FILE_PATH}"
conda activate "${FLASHFOLD_ENV_NAME}"

# Download the updater
wget -qnc -O "$FLASHFOLD_CONDA_ENV_DIR/update_M1mac.sh" \
    https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/update_M1mac.sh
chmod +x "$FLASHFOLD_CONDA_ENV_DIR/update_M1mac.sh"

# Download weights
"$FLASHFOLD_CONDA_ENV_DIR/conda/envs/$FLASHFOLD_ENV_NAME/bin/python3" -m colabfold.download
echo "Download of alphafold2 weights finished."
echo "-----------------------------------------"
echo "Installation of FoldFlash finished."