from setuptools import setup, find_packages

setup(
    name='flashfold',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'setuptools~=69.0.3',
        'wget~=3.2',
        'biopython~=1.83',
        'ujson~=5.10.0',
        'numpy~=1.26.4',
        'pandas~=2.2.2',
        'requests~=2.32.3',
        'tqdm~=4.66.5',
        'colabfold[alphafold]',
        'jax[cuda12]==0.4.28',
        'tensorflow',
        'silence_tensorflow'
    ],
    entry_points={
        'console_scripts': [
            'flashfold = flashfold.main:main'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)