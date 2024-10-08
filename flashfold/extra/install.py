import subprocess
import sys
from pathlib import Path


def check_and_install(command, brew_package):
    if subprocess.call(['which', command], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
        print(f"\t- {command} is installed.")
    else:
        if subprocess.call(['which', 'brew'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
            choice = input(f"\t- Do you want {command} to be installed now using Homebrew (recommended)? (y/n): ")
            if choice.lower() == 'y':
                subprocess.check_call(['brew', 'install', brew_package])
            else:
                print(f"\t- Please install {command} using: 'brew install {brew_package}' (recommended) and re-run this script.")
                sys.exit(1)
        else:
            print(f"\t- Warning: {command} is not installed (but required).")
            print(f"\t\t- Please install {command} using: 'brew install {brew_package}' (recommended) and re-run this script.")
            print(f"\t\t- Homebrew (brew) is not installed either.")
            sys.exit(1)


def main():
    current_path = Path.cwd()
    uname_out = subprocess.check_output(['uname', '-s']).strip().decode()
    arch_name = subprocess.check_output(['uname', '-m']).strip().decode()

    if uname_out == 'Darwin':
        if arch_name == 'x86_64':
            os_type = 'mac-intel'
            env_yml_file_path = current_path / 'envs' / 'mac-intel-environment.yml'
        elif arch_name == 'arm64':
            os_type = 'mac-silicon'
            subprocess.check_call(['ulimit', '-n', '99999'])
            env_yml_file_path = current_path / 'envs' / 'mac-silicon-environment.yml'
        else:
            print(f"Unknown architecture: {arch_name}")
            sys.exit(1)
    elif uname_out == 'Linux':
        os_type = 'linux'
        env_yml_file_path = current_path / 'envs' / 'linux-environment.yml'
    else:
        print(f"Unsupported OS: {uname_out}")
        sys.exit(1)

    if uname_out == 'Darwin':
        check_and_install('wget', 'wget')
        check_and_install('hhsearch', 'brewsci/bio/hh-suite')
        check_and_install('kalign', 'kalign')
        check_and_install('jackhmmer', 'hmmer')

    try:
        subprocess.check_call(['conda', '--version'])
    except subprocess.CalledProcessError:
        print("\t- Warning: conda is not installed (but required).")
        print("\t- Please first install conda (as instructed in the README) and later, re-run this script.")
        sys.exit(1)

    current_conda_path = subprocess.check_output(['conda', 'info', '--base']).strip().decode()
    flashfold_env_name = 'flashfold'
    flashfold_conda_env_dir = Path(current_conda_path) / 'envs' / flashfold_env_name

    subprocess.check_call(['conda', 'env', 'update', '-f', str(env_yml_file_path)])

    if os_type == 'mac-silicon':
        subprocess.check_call(['wget', '-qnc', '-O', str(flashfold_conda_env_dir / 'update_M1mac.sh'),
                               'https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/update_M1mac.sh'])
        (flashfold_conda_env_dir / 'update_M1mac.sh').chmod(0o755)

    if os_type == 'mac-intel':
        subprocess.check_call(['wget', '-qnc', '-O', str(flashfold_conda_env_dir / 'update_intelmac.sh'),
                               'https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/update_intelmac.sh'])
        (flashfold_conda_env_dir / 'update_intelmac.sh').chmod(0o755)

    if os_type == 'linux':
        subprocess.check_call(['wget', '-qnc', '-O', str(flashfold_conda_env_dir / 'update_linux.sh'),
                               'https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/update_linux.sh'])
        (flashfold_conda_env_dir / 'update_linux.sh').chmod(0o755)

    subprocess.check_call([str(flashfold_conda_env_dir / 'bin' / 'python3'), '-m', 'colabfold.download'])
    print("Download of alphafold2 weights finished.")
    print("-----------------------------------------")
    print("\nInstallation of FlashFold dependencies is complete.\n")


if __name__ == '__main__':
    main()