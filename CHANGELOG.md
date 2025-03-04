# Changelog

All notable changes to this project will be documented in this file.


## [1.1.2] - 2025-03-04

#### Added
- NA

#### Changed
- Version number is updated.
- Fastens MSA filtering process.

#### Fixed
- Processing of get_combined_a3m_records function that combines A3M record after deduplication.

#### Removed
- NA


## [1.1.1] - 2025-03-03

#### Added
- Tempdir for storing CD-HIT output files.

#### Changed
- Version number is updated.

#### Fixed
- Deletes temp file upon finishing the deduplication.

#### Removed
- NA


## [1.1.0] - 2025-02-28

#### Added
- Python script to reformat alignment files.
- CD-HIT based diversity aware deduplication method.
- Optional argument for actual interface pTM (actifpTM). 

#### Changed
- Version number is updated.
- Alignment formatting and deduplication process.
- Summary report is updated.
- Relaxed input checking criteria for MSA (a3m) file.

#### Fixed
- Reduction of hits when required, while maintaining diversity.

#### Removed
- Perl script for reformating the stockholm file.


## [1.0.9] - 2025-02-20

#### Added
- Shows relaxed structure in summary when available.
- Shows log message for batch query.

#### Changed
- Version number is updated.
- Summary report is updated.

#### Fixed
- Log message for batch query.

#### Removed
- NA


## [1.0.8] - 2025-02-15

#### Added
- Filtering when hits are too many.

#### Changed
- Version number is updated.
- MSA construction is updated for both monomer and multimer prediction.

#### Fixed
- A3M file size issue is fixed for complex prediction.
- Issue with downloading the database.
- Issue with the html table sorting in the table.
- Issue with validation of the input files (DNA/RNA).

#### Removed
- A function named get_to_be_downloaded_files.


## [1.0.7] - 2025-02-07

#### Added
- Functionality for generating an interactive html report.

#### Changed
- Version number is updated.
- Headers of the score.tsv file are updated.
- README.md file is updated with application of the `summary` subcommand.

#### Fixed
- Scoring table for the monomers.
- Time record log for batch processing.

#### Removed
- NA


## [1.0.6] - 2025-01-30

#### Added
- NA

#### Changed
- Version number is updated.

#### Fixed
- Unnecessary database loading issue is fixed.

#### Removed
- NA


## [1.0.5] - 2025-01-26

#### Added
- Doi from biorxiv to README.md file.

#### Changed
- README.md file is updated.
- Version number is updated.
- Git workflow file name is updated.

#### Fixed
- NA

#### Removed
- NA


## [1.0.4] - 2025-01-22

#### Added
- Added link of the databases to json file.

#### Changed
- README.md file is updated.

#### Fixed
- NA

#### Removed
- NA


## [1.0.3] - 2024-11-02

#### Added
- Version number is added to the help message.

#### Changed
- README.md file is updated with more examples.

#### Fixed
- Alignment file copy issue is fixed.
- A3M file input issue (database requirement) is fixed.

#### Removed
- NA


## [1.0.2] - 2024-10-28

#### Added
- `test` folder with test files.
- Used can run test scripts using `pytest` command.
- Logo added to the README.md file.

#### Changed
- Updated the README.md file with more examples.
- Updated the arguments for `ncbi_data` subcommand. 
- Collects `Hmmer` from conda biocore for macOS.
- Version number is updated.

#### Fixed
- Cutoff value is now only applicable for multimer prediction.
- All subcommands are fixed.

#### Removed
- No longer depends on `HH-suite`, `Kalign` and `mmseqs2`.


## [1.0.1] - 2024-10-09

#### Added
- CHANGELOG.md file to keep track of changes.

#### Changed
- Updated the installation instructions in the README.md file.
- updated the `install.sh` script.
- Version number is updated.

#### Fixed
- N/A

#### Removed
- N/A


## [1.0.0] - 2024-10-08

#### Added
- Initial release of FlashFold.

#### Changed
- N/A

#### Fixed
- N/A

#### Removed
- N/A