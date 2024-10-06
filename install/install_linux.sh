#!/bin/bash -e

CURRENTPATH=`pwd`
CURRENT_CONDA_PATH=$(conda info --base)
ENV_YML_FILE_PATH="${CURRENTPATH}/envs/linux-environment.yml"
FLASHFOLD_ENV_NAME="flashfold"
FLASHFOLD_CONDA_ENV_DIR="${CURRENT_CONDA_PATH}/envs/${FLASHFOLD_ENV_NAME}"
FLASHFOLD_ENV_SITE_PACKAGES="${FLASHFOLD_CONDA_ENV_DIR}/lib/python3.10/site-packages/"


conda env update -f "${ENV_YML_FILE_PATH}"
conda activate "${FLASHFOLD_ENV_NAME}"

# Download the updater
wget -qnc -O "$FLASHFOLD_CONDA_ENV_DIR/update_linux.sh" \
    https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/update_linux.sh
chmod +x "$FLASHFOLD_CONDA_ENV_DIR/update_linux.sh"

# Download weights
"${FLASHFOLD_CONDA_ENV_DIR}/bin/python3" -m colabfold.download
echo "Download of alphafold2 weights finished."
echo "-----------------------------------------"
echo "Installation of FoldFlash finished."