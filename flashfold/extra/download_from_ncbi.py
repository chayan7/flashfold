# Author: Chayan Kumar Saha

import argparse
import datetime
import os
import sys
import subprocess
from subprocess import Popen, PIPE
from typing import List, Dict
import shutil
import time


def summary_error_message(stderr: str, input_taxon: str) -> None:
    """
    Parses the stderr string to extract and print suggested taxa names when an error occurs.

    Args:
        stderr (str): The standard error output containing error messages and suggested taxa names.
        input_taxon (str): The taxon name input by the user that caused the error.

    Returns:
        None
    """
    suggested_taxa: List[str] = []
    for line in stderr.splitlines():
        if line.startswith("Error:") or line.startswith("Use datasets"):
            pass
        elif line == "":
            pass
        else:
            taxon_name = line.split("(")[0].strip()
            taxon = f"\t- {taxon_name}"
            suggested_taxa.append(taxon)
    print(f"\nError: '{input_taxon}' - is incorrect. Please try again with one of the suggested organism names:")
    print("\n".join(suggested_taxa))
    return None


def get_output_dir_path(output: str) -> str:
    """
    Ensures the specified output directory exists and is empty, creating it if necessary.

    Args:
        output (str): The path to the output directory.

    Returns:
        str: The absolute path to the output directory.

    Raises:
        SystemExit: If the specified output directory is not empty.
    """
    if os.path.exists(output):
        if not any(os.scandir(output)):
            return os.path.abspath(output)
        else:
            print(f"Provided output directory '{output}' is not empty, please use a different output directory.")
            sys.exit()
    else:
        print(f"\nCreating output directory: '{output}'\n")
        os.makedirs(output)
        return os.path.abspath(output)


def get_accessions_from_summary(summary_file: str) -> List[str]:
    """
    Extracts accession numbers from a summary file.

    This function reads a summary file and extracts accession numbers that start with "GCF_" or "GCA_".
    The accession numbers are expected to be the first element in a tab-separated line.

    Args:
        summary_file (str): The path to the summary file.

    Returns:
        List[str]: A list of accession numbers extracted from the summary file.
    """
    accessions: List[str] = []
    with open(summary_file, 'r') as sum_file:
        for line in sum_file:
            if line.startswith("GCF_") or line.startswith("GCA_"):
                accession = line.split("\t")[0]
                accessions.append(accession)
    return accessions


def create_batches(main_set_of_items: list, max_items_per_batch: int) -> list:
    """
    Splits a list of items into smaller batches of a specified maximum size.

    Args:
        main_set_of_items (list): The list of items to be divided into batches.
        max_items_per_batch (int): The maximum number of items allowed in each batch.

    Returns:
        list: A list of batches, where each batch is a list containing up to `max_items_per_batch` items.
    """
    subsets: list = []
    current_subset: list = []
    for item in main_set_of_items:
        current_subset.append(item)
        if len(current_subset) == max_items_per_batch:
            subsets.append(current_subset)
            current_subset = []
    if current_subset:
        subsets.append(current_subset)
    return subsets


def parse_formats(input_formats: str) -> List[str]:
    """
    Parses a comma-separated string of file formats and validates them against a list of valid formats.

    Args:
        input_formats (str): A comma-separated string of file formats.

    Returns:
        List[str]: A list of valid file formats.

    Raises:
        argparse.ArgumentTypeError: If any of the provided formats are invalid.

    Notes:
        If the input string contains 'all', the function returns the entire list of valid formats.
    """
    formats = input_formats.split(',')
    for file_format in formats:
        format_name = file_format.strip()
        if format_name not in valid_formats:
            if format_name == 'all':
                return valid_formats
            else:
                raise argparse.ArgumentTypeError(
                    f"Invalid format: '{file_format}'.\nValid choices are: {', '.join(valid_formats)}.")
    return [file_format.strip() for file_format in formats]


def process_command_batch(command_batch: List[str]) -> bool:
    """
    Executes a batch of shell commands and retries any failed commands.

    Args:
        command_batch (List[str]): A list of shell commands to be executed.

    Returns:
        bool: True if all commands succeed, False if any command fails after retries.
    """
    failed_commands = []
    process_list = [Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) for cmd in command_batch]
    for process in process_list:
        process.communicate()
        if process.returncode != 0:
            print(f"\t\tError in command: {process.args}")
            failed_commands.append(process.args)
    while failed_commands:
        current_failed_commands = failed_commands
        failed_commands = []
        for cmd in current_failed_commands:
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            while process.poll() is None:
                time.sleep(0.1)  # Wait a bit before checking again
            process.communicate()
            if process.returncode != 0:
                failed_commands.append(process.args)
            else:
                print(f"\t\tCommand succeeded on retry: {process.args}")
    return len(failed_commands) == 0


# Define valid choices
valid_formats = ['genome', 'protein', 'cds', 'gff3', 'gtf', 'gbff']

format_to_search_ext = {
    'genome': 'GC*_genomic.fna',
    'protein': '*.faa',
    'cds': 'cds_*.fna',
    'gff3': '*.gff',
    'gtf': '*.gtf',
    'gbff': '*.gbff',
}

format_to_file_ext = {
    'genome': '_genomic.fna',
    'protein': '.faa',
    'cds': '_cds.fna',
    'gff3': '.gff',
    'gtf': '.gtf',
    'gbff': '.gbff',
}


def main() -> None:
    usage = '''
        Description: 

            This script retrieves genome assembly data from NCBI. 
            Users can download the data in the following formats:
            (1) genome - Genomic sequences
            (2) protein - Amino acid sequences
            (3) cds - Nucleotide coding sequences
            (4) gff3 - General feature file (GFF3)
            (5) gtf - Gene transfer format (GTF)
            (6) gbff - GenBank flat file (GBFF)
        '''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-a", "--api-key", metavar="<String>", help="Specify an NCBI API key. "
                                                                    "To get this key kindly check: https://"
                                                                    "ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/"
                                                                    "new-api-keys-for-the-e-utilities/")

    # # Create a mutually exclusive group for -i and -n
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", metavar="<FILE_In>",
                       help="a text file that contains one NCBI assembly accession per line.")
    group.add_argument("-n", "--name", metavar="<String>",
                       help="name of organism eg., 'Pseudomonas aeruginosa' etc.")
    parser.add_argument("-f", "--format", metavar="<String>", type=parse_formats, required=True,
                        help="file format to be downloaded (comma-separated for multiple). "
                             "Valid options: 'genome', 'protein', 'cds', 'gff3', 'gtf', 'gbff', 'all'.")
    parser.add_argument("-o", "--output", metavar="<Output_Dir>", required=True,
                        help="path that will contain output.")
    parser.add_argument("-s", "--source", metavar="<String>", choices=['RefSeq', 'GenBank', 'all'],
                        default='all', help="assembly source selection. Valid options: 'RefSeq', 'GenBank' or 'all' "
                                            "(default).")
    args = parser.parse_args()

    today = datetime.date.today()
    output_directory = get_output_dir_path(args.output)

    api_key_parameter = f"--api-key {args.api_key}" if args.api_key else ""
    batch_length = 8 if args.api_key else 2
    include_parameter = f"--include {','.join(args.format)}" if args.format else ""

    assembly_accessions: List[str] = []
    if args.input:
        with open(args.input, 'r') as acc_in:
            for line in acc_in:
                if line.startswith("GCA_") or line.startswith("GCF_"):
                    assembly_accession = line.split("\t")[0].strip()
                    if assembly_accession not in assembly_accessions:
                        assembly_accessions.append(assembly_accession)
    elif args.name:
        summary_file_path = os.path.join(output_directory, "summary.tsv")
        summary_fields = (f"accession,organism-name,organism-tax-id,assminfo-level,checkm-completeness,"
                          f"assminfo-release-date,checkm-completeness,source_database")
        summary_command = (f"datasets summary genome taxon '{args.name}' --assembly-source {args.source} "
                           f"--as-json-lines {api_key_parameter} | dataformat tsv genome --fields {summary_fields} "
                           f"> {summary_file_path}")
        try:
            summary = subprocess.run(summary_command, shell=True, capture_output=True, text=True)
            if summary.stderr:
                summary_error_message(summary.stderr, args.name)
            else:
                ncbi_accessions = get_accessions_from_summary(summary_file_path)
                assembly_accessions.extend(ncbi_accessions)
        except subprocess.CalledProcessError as e:
            print(f"Error running datasets: {e}")

    if len(assembly_accessions) == 0:
        print("No assembly accessions found. Please try again with a valid input.")
        sys.exit()

    accession_counter = 0

    log_file_name = f"download_{today.strftime('%Y%m%d')}.log"
    log_file_path = os.path.join(output_directory, log_file_name)

    for file_format in args.format:
        out_dir_per_file_format = os.path.join(output_directory, f"{file_format}")
        os.makedirs(out_dir_per_file_format)

    temp_dir = os.path.join(output_directory, "temp")
    os.makedirs(temp_dir)

    command_to_zip_file: Dict[str, str] = {}
    for accession in assembly_accessions:
        zip_file_path = os.path.join(temp_dir, f"{accession}.zip")
        download_command = (f"datasets download genome accession {accession} {include_parameter} "
                            f"{api_key_parameter} --filename {zip_file_path}")
        command_to_zip_file[download_command] = zip_file_path

    total_command_batches = create_batches(list(command_to_zip_file.keys()), batch_length)

    with open(log_file_path, 'w') as log_file:
        for command_batch in total_command_batches:
            is_downloaded = process_command_batch(command_batch)
            if is_downloaded:
                for command in command_batch:
                    zip_file_path = command_to_zip_file[command]
                    accession = os.path.basename(zip_file_path).split('.')[0]
                    if os.path.exists(zip_file_path):
                        accession_counter += 1
                        print(f"\t-- Downloaded {accession} ({accession_counter}/{len(assembly_accessions)})")
                        for file_format in args.format:
                            search_extension = format_to_search_ext[file_format]
                            file_ext = format_to_file_ext[file_format]
                            unzip_command = (
                                f"unzip -p {zip_file_path} ncbi_dataset/data/{accession}/{search_extension} "
                                f"> {output_directory}/{file_format}/{accession}{file_ext}")
                            try:
                                subprocess.run(unzip_command, shell=True, capture_output=True, text=True)
                            except subprocess.CalledProcessError as e:
                                print(f"Error while preparing {accession} {file_format}: {e}")
                        # noinspection PyTypeChecker
                        print(accession, file=log_file)
                    else:
                        print(f"Error while downloading {accession}")

    shutil.rmtree(temp_dir)
    print("Done!" if accession_counter == len(assembly_accessions) else "Completed with some error!")
    return None


if __name__ == "__main__":
    main()
