from setuptools import setup, find_packages

setup(
    name='flashfold',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        "flashfold"
    ],
    entry_points={
        'console_scripts': [
            'flashfold = flashfold.main:main'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Linux',
    ],
    python_requires='>=3.6',
)
