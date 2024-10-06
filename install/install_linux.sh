#!/bin/bash -e

CURRENTPATH=`pwd`
CURRENT_CONDA_PATH=$(conda info --base)
ENV_YML_FILE_PATH="${CURRENTPATH}/envs/linux-environment.yml"
FLASHFOLD_ENV_NAME="flashfold"
FLASHFOLD_CONDA_ENV_DIR="${CURRENT_CONDA_PATH}/envs/${FLASHFOLD_ENV_NAME}"
FLASHFOLD_ENV_SITE_PACKAGES="${FLASHFOLD_CONDA_ENV_DIR}/lib/python3.10/site-packages/"


conda env create -f "${ENV_YML_FILE_PATH}"
conda activate "${FLASHFOLD_ENV_NAME}"

# Download the updater
# Download the updater
wget -qnc -O "$FLASHFOLD_CONDA_ENV_DIR/update_linux.sh" \
    https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/update_linux.sh
chmod +x "$FLASHFOLD_CONDA_ENV_DIR/update_linux.sh"

pushd "${FLASHFOLD_ENV_SITE_PACKAGES}/colabfold"
# Use 'Agg' for non-GUI backend
sed -i -e "s#from matplotlib import pyplot as plt#import matplotlib\nmatplotlib.use('Agg')\nimport matplotlib.pyplot as plt#g" plot.py
# modify the default params directory
sed -i -e "s#appdirs.user_cache_dir(__package__ or \"colabfold\")#\"${FLASHFOLD_ENV_SITE_PACKAGES}/colabfold\"#g" download.py
# suppress warnings related to tensorflow
sed -i -e "s#from io import StringIO#from io import StringIO\nfrom silence_tensorflow import silence_tensorflow\nsilence_tensorflow()#g" batch.py
# remove cache directory
rm -rf __pycache__
popd

# Download weights
"$FLASHFOLD_CONDA_ENV_DIR/bin/python3" -m colabfold.download
echo "Download of alphafold2 weights finished."
echo "-----------------------------------------"
echo "Installation of FoldFlash finished."