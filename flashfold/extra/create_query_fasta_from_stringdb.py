# Author: Chayan Kumar Saha

import argparse
import os
import sys
import shutil
# noinspection PyPackageRequirements
from Bio import SeqIO
from typing import List, Set, Dict


def manage_output_dir(user_provided_dir: str, clean: bool) -> str:
    if not os.path.exists(user_provided_dir):
        os.makedirs(user_provided_dir)
        print(f"Created new output directory: {os.path.abspath(user_provided_dir)}")
        return os.path.abspath(user_provided_dir)
    else:
        if not os.path.isdir(user_provided_dir):
            raise ValueError(
                f"Error: The provided output path '{user_provided_dir}' is not a directory.")
        else:
            # List the directory contents; True if the directory is empty, False otherwise.
            if not any(os.scandir(user_provided_dir)):
                return os.path.abspath(user_provided_dir)
            else:
                if not clean:
                    print(f"Error: The provided output path: {os.path.abspath(user_provided_dir)} is not empty. ")
                    sys.exit()
                else:
                    for filename in os.listdir(user_provided_dir):
                        file_path = os.path.join(user_provided_dir, filename)
                        shutil.rmtree(file_path)
                    return os.path.abspath(user_provided_dir)


def make_fasta_from_pairs(out_path: str, list_of_pairs: List[List[str]], seq_dict: Dict) -> None:
    if len(list_of_pairs) > 0:
        os.makedirs(out_path)
        for pair in list_of_pairs:
            prot_1 = pair[0]
            prot_2 = pair[1]
            prot_1_id = prot_1.split(".")[1]
            prot_2_id = prot_2.split(".")[1]
            fasta_file_name = f"{prot_1_id}_{prot_2_id}.fasta"
            fasta_path = os.path.join(out_path, fasta_file_name)
            with open(fasta_path, "w") as exp_out:
                # noinspection PyTypeChecker
                print(f">{prot_1_id}\n{seq_dict[prot_1]}\n>{prot_2_id}\n{seq_dict[prot_2]}",
                      file=exp_out)
    else:
        pass


def create_id_pairs(a: str, b: str) -> List[List[str]]:
    return [[a, b], [b, a]]


def create_new_pair_fasta(out_dir: str, query: str, relevant_seq_id_list: List[str], irrelevant_seq_id_list: List[str],
                          sequence_dict: Dict, limit: int) -> None:

    potential_partner_ids = []
    for seq_id in sorted(sequence_dict):
        if seq_id not in relevant_seq_id_list:
            potential_partner_ids.append(seq_id)
    
    for irrelevant_id in irrelevant_seq_id_list:
        if irrelevant_id not in potential_partner_ids:
            potential_partner_ids.append(irrelevant_id)
            
    candidate_limit = len(potential_partner_ids) if limit == 0 else limit
    potential_partner_candidates = potential_partner_ids[:candidate_limit]

    candidate_pair_list: List[List[str]] = []
    for partner_id in potential_partner_candidates:
        if "." not in partner_id:
            print("Error: protein id is not correct.")
            sys.exit()
        prot_1 = f"{partner_id.split('.')[0]}.{query}"
        prot_2 = partner_id
        for pair in create_id_pairs(prot_1, prot_2):
            candidate_pair_list.append(pair)

    make_fasta_from_pairs(out_dir, candidate_pair_list, sequence_dict)


def main():
    usage = '''Description: Get Assembly with column 1 GCF number and final column FTP'''

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-f", "--fasta", metavar="<FILE_In>", required=True,
                        help="Fasta file from String DB that contains all protein sequences of a species. "
                             "Hints: 208964.protein.sequences.v12.0.fa")
    parser.add_argument("-l", "--link", metavar="<FILE_In>", required=True,
                        help="Detailed link file from String DB that contains all protein links of a species. "
                             "Hints: 208964.protein.links.full.v12.0.txt")
    parser.add_argument("-n", "--number", metavar="<Integer>", type=int,
                        help="Number of query fasta files to be generated")
    parser.add_argument("-q", "--query", metavar="<String>", type=str, required=True,
                        help="Name of query protein eg. PA001 or PA1718")
    parser.add_argument("-o", "--output", metavar="<Output_Dir>", required=True,
                        help="Path to output directory")
    parser.add_argument("-d", "--delete", metavar="<Boolean>", type=bool, default=False,
                        help="Delete existing files in the provided output directory")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s 1.0.1')
    args = parser.parse_args()

    id_sequence_dict: Dict[str, str] = {}
    parsed_sequence = SeqIO.parse(args.fasta, "fasta")
    for record in parsed_sequence:
        id_sequence_dict[str(record.id)] = str(record.seq)

    query = args.query
    limit = args.number if args.number else 0

    output_dir = manage_output_dir(args.output, args.delete)
    exp_out_dir = os.path.join(output_dir, "validated_query_fasta")
    predicted_out_dir = os.path.join(output_dir, "predicted_query_fasta")
    new_out_dir = os.path.join(output_dir, "new_query_fasta")

    exp_validated_links: List[List[str]] = []
    predicted_links: List[List[str]] = []
    relevant_seq_ids: Set[str] = set()
    irrelevant_seq_ids: Set[str] = set()
    relevant_seq_id_list: List[str] = []
    irrelevant_seq_id_list: List[str] = []

    with open(args.link, "r") as link_file:
        """
        Header of the link file:
        
        0    protein1
        1    protein2
        2    neighborhood
        3    neighborhood_transferred
        4    fusion cooccurence
        5    homology
        6    coexpression
        7    coexpression_transferred
        8    experiments
        9    experiments_transferred
        10   database database_transferred
        11   textmining
        12   textmining_transferred
        13   combined_score
        
        """
        for line in link_file:
            if line[0] != "#":
                if query in line:
                    split_line = line.rstrip().split()
                    exp_score = int(split_line[8])
                    experiments_transferred_score = int(split_line[9])
                    score = exp_score + experiments_transferred_score
                    if score != 0:
                        exp_validated_links.append(split_line[:2])
                        if split_line[0] not in relevant_seq_ids:
                            relevant_seq_id_list.append(split_line[0])
                            relevant_seq_ids.add(split_line[0])
                        if split_line[1] not in relevant_seq_ids:
                            relevant_seq_id_list.append(split_line[1])
                            relevant_seq_ids.add(split_line[1])
                    else:
                        predicted_links.append(split_line[:2])
                        if split_line[0] not in relevant_seq_ids:
                            relevant_seq_id_list.append(split_line[0])
                            relevant_seq_ids.add(split_line[0])
                        if split_line[1] not in relevant_seq_ids:
                            relevant_seq_id_list.append(split_line[1])
                            relevant_seq_ids.add(split_line[1])
                else:
                    if not line.startswith("protein"):
                        irrelevant_split_line = line.rstrip().split()
                        if irrelevant_split_line[0] not in irrelevant_seq_ids:
                            irrelevant_seq_id_list.append(irrelevant_split_line[0])
                            irrelevant_seq_ids.add(irrelevant_split_line[0])
                        if irrelevant_split_line[1] not in irrelevant_seq_ids:
                            irrelevant_seq_id_list.append(irrelevant_split_line[1])
                            irrelevant_seq_ids.add(irrelevant_split_line[1])

    make_fasta_from_pairs(exp_out_dir, exp_validated_links, id_sequence_dict)
    make_fasta_from_pairs(predicted_out_dir, predicted_links, id_sequence_dict)
    create_new_pair_fasta(new_out_dir, query, relevant_seq_id_list, irrelevant_seq_id_list, id_sequence_dict, limit)


if __name__ == "__main__":
    main()
