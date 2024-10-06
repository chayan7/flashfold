# Author: Chayan Kumar Saha

import argparse
import wget
import os.path
import sys
import datetime
from datetime import date
from collections import defaultdict
from typing import TextIO

today = datetime.date.today()

# Get the current directory
current_directory = os.getcwd()


def parse_strain(value: str) -> int:
    try:
        n = int(value)
        if n < 1:
            raise argparse.ArgumentTypeError("Value must be a positive integer greater than or equal to 1")
        else:
            return n
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid integer value")


def parse_sub_directory(value: str) -> str:
    input_sub_directory = value.lower()
    desired_sub_directory = ["bacteria", "archaea", "fungi", "invertebrate", "plant", "protozoa",
                             "vertebrate_mammalian", "vertebrate_other"]
    if input_sub_directory not in desired_sub_directory:
        raise argparse.ArgumentTypeError("Input should be any of them: bacteria/archaea/fungi/invertebrate/plant/"
                                         "protozoa/vertebrate_mammalian/vertebrate_other")
    else:
        return input_sub_directory


def is_substring_present(target_string: str, query_string: str):
    if query_string in target_string:
        return True
    else:
        return False


def get_include_assemblies(input_file: TextIO) -> set:
    assembly_set = set()
    for each_line in input_file:
        processed_line = each_line.rstrip().split("\t")
        assembly_set.add(processed_line[0])
    input_file.close()
    return assembly_set


def wget_file_from_url(url: str) -> None:
    # Download the file to the current directory
    try:
        filename = wget.download(url, out=current_directory)
        print(f"\nFile '{filename}' downloaded successfully to {current_directory}.")
    except Exception as e:
        print(f"\nError downloading file: {e}")


# Function to get the updated seqs that have been released
def released_day_from_today(end: str) -> int:
    [year, month, day] = map(int, end.split('/'))
    end_date = date(year, month, day)
    return (today - end_date).days


def get_must_include_assembly(assembly_list: list, assembly_set_for_inclusion: set) -> list:
    include_list = []
    for assembly in assembly_list:
        if assembly in assembly_set_for_inclusion:
            include_list.append(assembly)
    return include_list


def has_sp_in_taxa(taxa: str) -> bool:
    sp_formats = ["Sp.", "sp.", "sp.,", "species"]
    taxa_name_parts = taxa.split(' ')[1:]
    return any(part in sp_formats for part in taxa_name_parts)


def get_best_updated_assemblies(gcf_accession_set: set) -> list:
    if len(gcf_accession_set) <= 1:
        return list(gcf_accession_set)
    else:
        gcf_to_assembly_level = {}
        gcf_to_released_date = {}
        for gcf_accession in gcf_accession_set:
            assembly_level = gcf_to_assembly_info[gcf_accession][11]
            released_date = gcf_to_assembly_info[gcf_accession][14]
            gcf_to_released_date[gcf_accession] = released_day_from_today(released_date)
            if assembly_level == "Complete Genome":
                gcf_to_assembly_level[gcf_accession] = 1
            elif assembly_level == "Chromosome":
                gcf_to_assembly_level[gcf_accession] = 2
            elif assembly_level == "Scaffold":
                gcf_to_assembly_level[gcf_accession] = 3
            elif assembly_level == "Contig":
                gcf_to_assembly_level[gcf_accession] = 4
        gcf_sorted_by_released_date: list = sorted(gcf_to_released_date, key=gcf_to_released_date.get)
        gcf_to_score = {}  # scoring based on release date and assembly level
        for gcf_accession in gcf_sorted_by_released_date:
            gcf_to_score[gcf_accession] = (0.5 * gcf_sorted_by_released_date.index(gcf_accession)
                                           + int(gcf_to_assembly_level[gcf_accession]))
        gcf_ranked_by_score: list[str] = sorted(gcf_to_score, key=gcf_to_score.get)
        return gcf_ranked_by_score


def get_short_listed_assembly(best_assembly_list: list, to_include_assembly_set: set, assembly_info: dict,
                              limit: int) -> list:
    shortlisted_assembly = get_must_include_assembly(best_updated_assemblies, to_include_assembly_set)
    for assembly in best_assembly_list:
        assembly_taxa = assembly_info[assembly][7]
        if not has_sp_in_taxa(assembly_taxa):
            if assembly not in shortlisted_assembly:
                shortlisted_assembly.append(assembly)

    if len(shortlisted_assembly) == 0:
        shortlisted_assembly.append(best_assembly_list[0])

    return shortlisted_assembly[:limit]


usage = ''' Description:  This program uses NCBI FTP and will output a tab delimited text file with the Accession 
Number of the latest and best possible sequence assembly; including some useful information. '''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--include", type=argparse.FileType('r'),
                    help="List of assembly number eg. from COG database or reference genomes that must be included in "
                         "the final output")
# FIXME: use options rather than type function
parser.add_argument("-g", "--get", type=parse_sub_directory, required=True,
                    help="Input should be any of them : bacteria/archaea/fungi/invertebrate/plant/protozoa/"
                         "vertebrate_mammalian/vertebrate_other")
parser.add_argument("-s", "--strain", type=parse_strain, default=1, help="Number of strain 1,2,3.. default 1")
parser.add_argument("-o", "--output", type=str, required=True, help="Output prefix [eg. output]")
parser.add_argument("-k", "--keep", action="store_true",
                    help=" If you want to keep the intermediate files eg. assembly_summary.txt use [-k]. By default it "
                         "will remove. ")
parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0.0')
args = parser.parse_args()
parser.parse_args()

number_of_strain_per_taxa = args.strain
query = args.get.lower()
include_assembly_set = get_include_assemblies(args.include)

gcf_to_assembly_info = {}  # create a dictionary where accessionNumber is a key and other all information as value
species_taxid_set = set()  # Species TaxId set created to get unique Species TaxId

# RefSeq file to be downloaded
refseq_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{query}/assembly_summary.txt"
refseq_file = os.path.join(current_directory, "assembly_summary.txt")

if os.path.isfile(refseq_file):
    pass
else:
    print(f"\nDownloading data from RefSeq database ... \n")
    wget_file_from_url(refseq_url)

with open('assembly_summary.txt', 'r') as fileIn:
    for line in fileIn:
        if line[0] != '#':
            if is_substring_present(line, "ftp.ncbi.nlm.nih.gov"):
                Line = line.rstrip().split('\t')
                gcf_to_assembly_info[Line[0]] = Line
                species_taxid_set.add(Line[6])
print(f"\n- Total unique species TaxId(s) = {len(species_taxid_set)}")

# Collect all assemblies per taxId
taxid_to_gcf_set = defaultdict(set)
for gcf in gcf_to_assembly_info:
    tax_id = gcf_to_assembly_info[gcf][6]
    taxid_to_gcf_set[tax_id].add(gcf)

# Shortlist the candidate assemblies per taxId based on criteria
taxid_to_shortlisted_gcf_set = defaultdict(list)
for tax_id, gcf_set in taxid_to_gcf_set.items():
    best_updated_assemblies = get_best_updated_assemblies(gcf_set)
    if len(best_updated_assemblies) == 1:
        if tax_id not in taxid_to_shortlisted_gcf_set:
            taxid_to_shortlisted_gcf_set[tax_id].extend(best_updated_assemblies)
    elif len(best_updated_assemblies) > 1:
        shortlisted_assemblies = get_short_listed_assembly(best_updated_assemblies, include_assembly_set,
                                                           gcf_to_assembly_info, number_of_strain_per_taxa)
        if tax_id not in taxid_to_shortlisted_gcf_set:
            taxid_to_shortlisted_gcf_set[tax_id].extend(shortlisted_assemblies)

# Write out the results
i = 0
out_file_name = f"{args.output}_{query}_{today.strftime('%Y%m%d')}.tsv"
with open(out_file_name, 'w') as file_out:
    # noinspection PyTypeChecker
    print('#Accession_Number', 'Taxonomy', 'TaxId', 'Species_TaxId', 'Organism_name', 'Infraspecific_name',
          'Assembly_level', 'Genome_rep', 'Seq_rel_date', 'ftp_Path', sep='\t', file=file_out)
    for tax_id in taxid_to_shortlisted_gcf_set:
        for gcf in taxid_to_shortlisted_gcf_set[tax_id]:
            i += 1
            # noinspection PyTypeChecker
            print(gcf_to_assembly_info[gcf][0], query, gcf_to_assembly_info[gcf][5], gcf_to_assembly_info[gcf][6],
                  gcf_to_assembly_info[gcf][7], gcf_to_assembly_info[gcf][8], gcf_to_assembly_info[gcf][11],
                  gcf_to_assembly_info[gcf][13], gcf_to_assembly_info[gcf][14], gcf_to_assembly_info[gcf][19],
                  sep='\t', file=file_out)

j = 0
out_file_name_no_sp = f"{args.output}_{query}_no_sp_{today.strftime('%Y%m%d')}.tsv"
with open(out_file_name_no_sp, 'w') as no_sp_file_out:
    # noinspection PyTypeChecker
    print('#Accession_Number', 'Taxonomy', 'TaxId', 'Species_TaxId', 'Organism_name', 'Infraspecific_name',
          'Assembly_level', 'Genome_rep', 'Seq_rel_date', 'ftp_Path', sep='\t', file=no_sp_file_out)
    for tax_id in taxid_to_shortlisted_gcf_set:
        for gcf in taxid_to_shortlisted_gcf_set[tax_id]:
            if not has_sp_in_taxa(gcf_to_assembly_info[gcf][7]):
                j += 1
                # noinspection PyTypeChecker
                print(gcf_to_assembly_info[gcf][0], query, gcf_to_assembly_info[gcf][5], gcf_to_assembly_info[gcf][6],
                      gcf_to_assembly_info[gcf][7], gcf_to_assembly_info[gcf][8], gcf_to_assembly_info[gcf][11],
                      gcf_to_assembly_info[gcf][13], gcf_to_assembly_info[gcf][14], gcf_to_assembly_info[gcf][19],
                      sep='\t', file=no_sp_file_out)

assembly_file = "./assembly_summary.txt"  # for deletion of the downloaded file from ftp

# if file exists, delete it
if os.path.isfile(assembly_file):
    if args.keep:
        pass
    else:
        os.remove(assembly_file)
else:
    print("Error: %s file not found" % assembly_file)

print(f"- Total filtered assemblies with minimum sp. labelled taxa = {i}")
print(f"- Total filtered assemblies without sp. labelled taxa = {j}")
print(f"\nThanks for waiting, your output file is ready to use!\n\tCheck output here: "
      f"{os.path.abspath(os.path.curdir)}/\n")
sys.exit()
