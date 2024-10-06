from setuptools import setup, find_packages, Command
from setuptools.command.install import install
import platform
import subprocess
import os


class CustomInstallCommand(install, Command):
    def run(self) -> None:
        os_name = platform.system()
        if os_name == "Darwin":  # macOS
            subprocess.check_call(["sh", "install/install_mac.sh"])
        elif os_name == "Linux":
            subprocess.check_call(["sh", "install/install_linux.sh"])
        else:
            raise OSError("Unsupported operating system: {}".format(os_name))

        install.run(self)


command_class = {"install": CustomInstallCommand}

setup(
    name='flashfold',
    description='A tool for protein and protein complex structure prediction.',
    version='1.0.0',
    author='Chayan Kumar Saha',
    author_email='chayan_saha@outlook.com',
    license='MIT',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'flashfold = flashfold.main:main'
        ]
    },
    zip_safe=False,
    cmdclass=command_class,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Linux',
    ],
    python_requires='>=3.6',
)
