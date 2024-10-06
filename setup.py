from setuptools import setup, find_packages, Command
from setuptools.command.install import install
import subprocess


class CustomInstallCommand(install, Command):
    def run(self) -> None:
        # Install dependencies
        subprocess.check_call(["sh", "install.sh"])
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
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
        "Operating System :: MacOS",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta",
    ],
    python_requires='>=3.6',
)
