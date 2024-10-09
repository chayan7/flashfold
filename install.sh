#!/bin/bash -e

# Function to check and install dependencies
check_and_install() {
    if command -v "$1" &> /dev/null; then
        echo -e "\t- $1 is installed."
    else
        if command -v brew &> /dev/null; then
            # shellcheck disable=SC2162
            read -p "\t- Do you want $1 to be installed now using Homebrew (recommended)? (y/n): " choice
            if [ "$choice" == "y" ]; then
                brew install "$2"
            else
                echo -e "\t- Please install $1 using: 'brew install $2' (recommended) and re-run this script."
                exit 1
            fi
        else
            echo -e "\t- Warning: $1 is not installed (but required)."
            echo -e "\t\t- Please install $1 using: 'brew install $2' (recommended) and re-run this script."
            echo -e "\t\t- Homebrew (brew) is not installed either."
            exit 1
        fi
    fi
}

CURRENT_PATH="$(pwd)"

uname_out="$(uname -s)"
# Check machine OS type
case "${uname_out}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${uname_out}"
esac

arch_name="$(uname -m)"
if [ "$machine" == "Mac" ]; then
    if [ "${arch_name}" = "x86_64" ]; then
        if [ "$(sysctl -in sysctl.proc_translated)" = "1" ]; then
            echo "Running on Rosetta 2"
        else
            echo "Running on native Intel"
        fi
        OS_TYPE="mac-intel"
        ENV_YML_FILE_PATH="${CURRENT_PATH}/envs/mac-intel-environment.yml"
    elif [ "${arch_name}" = "arm64" ]; then
        echo "Running on Apple Silicon"
        OS_TYPE="mac-silicon"
        ulimit -n 99999
        ENV_YML_FILE_PATH="${CURRENT_PATH}/envs/mac-silicon-environment.yml"
    else
        echo "Unknown architecture: ${arch_name}"
        exit 1
    fi
elif [ "$machine" == "Linux" ]; then
    echo "Running on Linux"
    OS_TYPE="linux"
    ENV_YML_FILE_PATH="${CURRENT_PATH}/envs/linux-environment.yml"
fi


# macOS specific commands
if [ "$machine" == "Mac" ]; then
    check_and_install wget wget
    check_and_install gnu-sed gnu-sed
    check_and_install hhsearch brewsci/bio/hh-suite
    check_and_install kalign brewsci/bio/kalign
    check_and_install jackhmmer hmmer
    if [ "$OS_TYPE" == "mac-silicon" ]; then
        check_and_install cmake cmake
    fi
fi

# Conda installation & version check
if command -v conda &> /dev/null; then
    current_version=$(conda --version | awk '{print $2}')
    echo -e "\t- Continuing with the existing conda version: $current_version"
else
    echo -e "\t- Warning: conda is not installed (but required)."
    echo -e "\t- Please first install conda (as instructed in the README) and later, re-run this script."
    exit 1
fi

# Conda environment setup
CURRENT_CONDA_PATH=$(conda info --base)
FLASHFOLD_ENV_NAME="flashfold"
FLASHFOLD_CONDA_ENV_DIR="${CURRENT_CONDA_PATH}/envs/${FLASHFOLD_ENV_NAME}"

# Update conda environment
conda env update -f "${ENV_YML_FILE_PATH}"

# Download weights
"${FLASHFOLD_CONDA_ENV_DIR}/bin/python3" -m colabfold.download
echo "Download of alphafold2 weights finished."
echo "-----------------------------------------"
echo -e "\nInstallation of FlashFold dependencies is complete.\n"