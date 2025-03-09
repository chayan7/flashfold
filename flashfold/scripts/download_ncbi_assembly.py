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
from flashfold.utils import wget_file_from_url, create_batches


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
    path_to_output = os.path.realpath(output)
    if os.path.exists(path_to_output):
        if not any(os.scandir(path_to_output)):
            return path_to_output
        else:
            print(f"Provided output directory '{path_to_output}' is not empty, please use a different directory.")
            sys.exit()
    else:

        print(f"\nCreating output directory: '{path_to_output}'\n")
        os.makedirs(path_to_output)
        return path_to_output


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


def parse_source(source: str) -> str:
    """
    Parses the source parameter and returns the corresponding source database name.

    Args:
        source (str): The source parameter provided by the user.

    Returns:
        str: The name of the source database.

    Raises:
        argparse.ArgumentTypeError: If the source parameter is invalid.
    """
    if source == 'refseq':
        return 'RefSeq'
    elif source == 'genbank':
        return 'GenBank'
    else:
        return source


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


def process_command_batch(command_batch: List[str]) -> Dict[str, List[str]]:
    """
    Executes a batch of shell commands.

    Args:
        command_batch (List[str]): A list of shell commands to be executed.

    Returns:
        Dict[str, List[str]]: A dictionary with keys 'processed' and 'failed' containing lists of commands.
    """
    failed_commands: List[str] = []
    processed_commands: List[str] = []
    process_list = [Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) for cmd in command_batch]
    for process in process_list:
        process.communicate()
        if process.returncode != 0:
            failed_process = str(process.args)
            print(f"\t\tError while download: {failed_process}")
            failed_commands.append(failed_process)
        else:
            successful_process = str(process.args)
            processed_commands.append(successful_process)

    # retry failed commands
    failed_commands_after_retry: List[str] = []
    for failed_command in failed_commands:
        print(f"\t\tRetrying failed attempt: {failed_command}")
        time.sleep(10)
        process = Popen(failed_command, shell=True, stdout=PIPE, stderr=PIPE)
        process.communicate()
        if process.returncode != 0:
            failed_process = str(process.args)
            print(f"\t\tError while retrying download: {failed_process}")
            failed_commands_after_retry.append(failed_process)
        else:
            successful_process = str(process.args)
            processed_commands.append(successful_process)

    return {'processed': processed_commands, 'failed': failed_commands_after_retry}


def is_reference_or_refseq_annotated(info_line: str, sub_dir: str) -> bool:
    """
    Checks if the info line contains the string 'reference genome' or 'NCBI RefSeq'.
    Args:
        info_line: string containing the information line from the assembly summary file.
        sub_dir:    string containing the subdirectory name of the assembly summary file.

    Returns:
        bool: True if the info line contains 'reference genome' or 'NCBI RefSeq', False otherwise.
    """
    search_string = "reference genome" if sub_dir != "viral" else "NCBI RefSeq"
    split_line = info_line.rstrip().split("\t")
    for item in split_line:
        if search_string == item:
            return True
    return False


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


def download_ncbi_data(args) -> None:
    if args.input and args.source:
        print("The --source (-s) can not be applied with --input (-i) , only applicable with --name (-n).")
        sys.exit()

    today = datetime.date.today()
    output_directory = get_output_dir_path(args.output)

    temp_dir = os.path.join(output_directory, "temp")
    os.makedirs(temp_dir)

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
        source = parse_source(args.source) if args.source else "RefSeq"
        summary_command = (f"datasets summary genome taxon '{args.name}' --assembly-source {source} "
                           f"--as-json-lines {api_key_parameter} | dataformat tsv genome --fields {summary_fields} "
                           f"> {summary_file_path}")
        try:
            summary = subprocess.run(summary_command, shell=True, capture_output=True, text=True)
            accessions = get_accessions_from_summary(summary_file_path)
            if summary.stderr.startswith("Error:") and len(accessions) == 0:
                summary_error_message(summary.stderr, args.name)
            else:
                assembly_accessions.extend(accessions)
        except subprocess.CalledProcessError as e:
            print(f"Error running datasets: {e}")
    elif args.reference:
        query = args.reference
        # RefSeq file to be downloaded
        refseq_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{query}/assembly_summary.txt"
        print(f"\nDownloading data from RefSeq database ... \n")
        refseq_file = wget_file_from_url(refseq_url, temp_dir)
        if os.path.isfile(refseq_file):
            with open(refseq_file, 'r') as refseq_in:
                for line in refseq_in:
                    if is_reference_or_refseq_annotated(line, query):
                        accession = line.split("\t")[0]
                        if accession.startswith("GCF_"):
                            assembly_accessions.append(accession)
    else:
        sys.exit()

    if len(assembly_accessions) == 0:
        print("No assembly accessions found. Please try again with a valid input.")
        sys.exit()

    accession_counter = 0

    log_file_name = f"download_{today.strftime('%Y%m%d')}.log"
    log_file_path = os.path.join(output_directory, log_file_name)

    for file_format in args.format:
        out_dir_per_file_format = os.path.join(output_directory, f"{file_format}")
        os.makedirs(out_dir_per_file_format)

    command_to_zip_file: Dict[str, str] = {}
    for accession in assembly_accessions:
        zip_file_path = os.path.join(temp_dir, f"{accession}.zip")
        download_command = (f"datasets download genome accession {accession} {include_parameter} "
                            f"{api_key_parameter} --filename {zip_file_path}")
        command_to_zip_file[download_command] = zip_file_path

    total_command_batches = create_batches(list(command_to_zip_file.keys()), batch_length)
    with open(log_file_path, 'w') as log_file:
        for command_batch in total_command_batches:
            results_after_processing = process_command_batch(command_batch)
            processed_commands = results_after_processing['processed']  # List of processed commands
            failed_commands = results_after_processing['failed']  # List of failed commands
            for processed_command in processed_commands:
                zip_file_path = command_to_zip_file[processed_command]
                accession = os.path.basename(zip_file_path)[:-4]
                if os.path.exists(zip_file_path):
                    accession_counter += 1
                    print(f"\t-- Downloaded {accession} ({accession_counter}/{len(assembly_accessions)})")
                    for file_format in args.format:
                        search_extension = format_to_search_ext[file_format]
                        file_ext = format_to_file_ext[file_format]
                        unzip_command = (
                            f"unzip -p {zip_file_path} ncbi_dataset/data/{accession}/{search_extension} "
                            f"> '{output_directory}/{file_format}/{accession}{file_ext}'")
                        try:
                            subprocess.run(unzip_command, shell=True, capture_output=True, text=True)
                        except subprocess.CalledProcessError as e:
                            print(f"Error while preparing {accession} {file_format}: {e}")
                    # noinspection PyTypeChecker
                    print(accession, "passed", sep="\t", file=log_file)
                else:
                    print(f"Error while downloading {accession}")
            for failed_command in failed_commands:
                print(f"\t\tFailed to download: {failed_command}")
                zip_file_path_failed = command_to_zip_file[failed_command]
                accession_failed = os.path.basename(zip_file_path_failed)[:-4]
                # noinspection PyTypeChecker
                print(accession_failed, "failed to download", sep="\t", file=log_file)

    shutil.rmtree(temp_dir)
    print("Done!" if accession_counter == len(assembly_accessions) else "Completed with some error!")
    return None
