# Author: Chayan Kumar Saha

import argparse
import ftplib
from ftplib import FTP
from urllib.parse import urlparse
import os
import os.path
import sys
import datetime
import wget

# Get the current directory
current_directory = os.getcwd()

# RefSeq and GenBank files to be downloaded
refseq_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
genbank_url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"

# Get today 
today = datetime.date.today()


def validate_ftp_link(url: str) -> dict:
    parsed_url = urlparse(url)
    netloc = parsed_url.netloc
    if netloc.startswith("ftp"):
        try:
            ftp_from_url = FTP(netloc)
            ftp_from_url.connect()
            ftp_from_url.quit()
            path = parsed_url.path
            user = "anonymous"
            acct = f"{user}@{netloc}"
            return {"netloc": netloc, "user": user, "acct": acct, "path": path}
        except Exception as e:
            raise argparse.ArgumentTypeError("FTP server validation failed:", e)
    else:
        raise argparse.ArgumentTypeError("FTP server validation failed.")
    

def get_cog_db_assemblies(cog_db_file_list: list[str]) -> set:
    cog_assemblies = set()
    if "cog-20.org.csv" in cog_db_file_list:
        ftp.retrbinary('RETR ' + 'cog-20.org.csv', open('cog-20.org.txt', 'wb').write)
        with open('cog-20.org.txt', 'r') as fileIn:
            for line in fileIn:
                if line[0] != '#':
                    processed_line = line.rstrip().split(',')
                    cog_assemblies.add(processed_line[0])
    else:
        pass
    
    os.remove("cog-20.org.txt")
    return cog_assemblies


def wget_file_from_url(url: str) -> None:
    # Download the file to the current directory
    try:
        filename = wget.download(url, out=current_directory)
        print(f"\nFile '{filename}' downloaded successfully to {current_directory}.")
    except Exception as e:
        print(f"\nError downloading file: {e}")


def get_refseq_acc_set_by_cog_acc(total_cog_assemblies: set) -> dict:
    refseq_file = os.path.join(current_directory, "assembly_summary_refseq.txt")
    genbank_file = os.path.join(current_directory, "assembly_summary_genbank.txt")
    
    cog_acc_to_refseq_acc = {}
    if not os.path.isfile(refseq_file):
        print(f"\nDownloading data from RefSeq database ... \n")
        wget_file_from_url(refseq_url)
    else:
        pass
    
    refseq_gcf_set = set()
    with open(refseq_file, 'r') as fileIn:
        for line in fileIn:
            if line[0] != '#':
                processed_line = line.rstrip().split('\t')
                refseq_gcf_set.add(processed_line[0])
                if processed_line[0] in total_cog_assemblies:
                    cog_acc_to_refseq_acc[processed_line[0]] = processed_line[0]
                elif processed_line[4] == "reference genome":
                    cog_acc_to_refseq_acc[processed_line[0]] = processed_line[0]
                    
    assemblies_missing_in_refseq = total_cog_assemblies - set(cog_acc_to_refseq_acc.keys())
    
    if len(assemblies_missing_in_refseq) > 0:
        if not os.path.isfile(genbank_file):
            print(f"\nDownloading data from GenBank database ... \n")
            wget_file_from_url(genbank_url)
        else:
            pass
        
        with open(genbank_file, 'r') as fileIn:
            for line in fileIn:
                if line[0] != '#':
                    processed_line = line.rstrip().split('\t')
                    if processed_line[0] in assemblies_missing_in_refseq and processed_line[18] == "identical":
                        if processed_line[17] in refseq_gcf_set:
                            cog_acc_to_refseq_acc[processed_line[0]] = processed_line[17]

    return cog_acc_to_refseq_acc


usage = '''Description: Uses information from the latest COG Database \
(https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/) and finds the corresponding RefSeq assembly accessions.  '''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-f", "--ftp", type=validate_ftp_link,
                    default="https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data/",
                    help="FTP link of updated COG database")
parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0.0')
args = parser.parse_args()
parser.parse_args()


ftp_link: dict = args.ftp


# FTP login to COG database
ftp = ftplib.FTP(ftp_link["netloc"], ftp_link["user"], ftp_link["acct"])
ftp.cwd(ftp_link["path"])  # move to data directory
name_list = ftp.nlst()  # get file/directory names within the directory
cog_assemblies_set = get_cog_db_assemblies(name_list)
cog_to_refseq_assemblies = get_refseq_acc_set_by_cog_acc(cog_assemblies_set)

cog_file = f"cog_in_refseq_{today.strftime('%Y%m%d')}.txt"
missed_cog_file = f"cog_missed_{today.strftime('%Y%m%d')}.txt"


with open(cog_file, "w") as file_out:
    for acc in cog_to_refseq_assemblies:
        if acc in cog_assemblies_set:
            # noinspection PyTypeChecker
            print(cog_to_refseq_assemblies[acc], "COG", sep="\t", file=file_out)
        else:
            # noinspection PyTypeChecker
            print(cog_to_refseq_assemblies[acc], "Reference", sep="\t", file=file_out)

i = 0
with open(missed_cog_file, "w") as missed_out:
    for acc in cog_assemblies_set:
        if acc not in cog_to_refseq_assemblies:
            i += 1
            # noinspection PyTypeChecker
            print(acc, file=missed_out)

print(f"Out of total {len(cog_assemblies_set)} COG assemblies, {i} could not be found in RefSeq and GenBank database.")

sys.exit()
