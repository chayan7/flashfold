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
        'tqdm~=4.66.5'
    ],
    dependency_links=[
        'git+https://github.com/sokrypton/ColabFold#egg=colabfold[alphafold]',
        'https://storage.googleapis.com/jax-releases/jax_cuda_releases.html'
    ],
    extras_require={
        'full': [
            'colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold',
            'jax[cuda]>=0.3.8,<0.4',
            'kalign2==2.04',
            'hhsuite==3.3.0',
            'hmmer',
            'openmm==8.0.0',
            'perl==5.32.1'
            'pdbfixer'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)