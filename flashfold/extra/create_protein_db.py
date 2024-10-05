# Author: Chayan Kumar Saha

import argparse
import os
import sys
from collections import defaultdict
from utils import extract_protein_sequences, is_valid_path, get_hash_to_files_with_extensions_from_dir, \
    load_json_file, write_dict_to_json_as_file, current_time, seq_to_fasta, write_dict_of_set_to_json_as_file


gbk_to_db_history = "history.json"
genbank_file_extensions = [".gbff", ".gbk"]


def update_gbk_file_dict(database_dir: str, gbk_file_dir: str) -> dict:
    if not is_valid_path(database_dir):
        os.makedirs(database_dir)

    history_file_path = os.path.join(os.path.abspath(database_dir), gbk_to_db_history)
    prev_gbk_files = load_json_file(history_file_path)
    new_gbk_files = get_hash_to_files_with_extensions_from_dir(gbk_file_dir, genbank_file_extensions)
    
    # Find new files
    to_be_updated_files = {file_hash: file_path for file_hash, file_path in new_gbk_files.items()
                           if file_hash not in prev_gbk_files}
    
    # Update stored directory info with current info
    if len(to_be_updated_files) > 0:
        new_gbk_files.update(prev_gbk_files)
        write_dict_to_json_as_file(new_gbk_files, history_file_path)
    return to_be_updated_files


def main():
    parser = argparse.ArgumentParser(description="Extract protein sequences from '.gbff' or '.gbk' "
                                                 "files and create sequence database.")
    parser.add_argument("-d", "--directory", required=True, help="Path to the directory that contains "
                                                                 "'.gbff' or '.gbk' files.")
    parser.add_argument("-db", "--database_dir", required=True, help="Path to the directory that will "
                                                                     "contain sequence database.")
    args = parser.parse_args()

    database_output_directory = args.database_dir
    print(f"\nStarted at {current_time()}")
    new_gbk_files = update_gbk_file_dict(database_output_directory, args.directory)
    if len(new_gbk_files) == 0:
        print(f"Database already contains sequences from: {os.path.abspath(args.directory)}\n")
        sys.exit()
    else:
        print(f"\n\tCreating database ... \n")
        history_file_path = os.path.join(os.path.abspath(database_output_directory), gbk_to_db_history)
        protein_hash_keys_to_gbk_hashes_file_path = os.path.join(os.path.abspath(database_output_directory),
                                                                 "protein_to_gbks.json")
        protein_hash_keys_to_accessions_file_path = os.path.join(os.path.abspath(database_output_directory),
                                                                 "prot_hash_to_accession.json")
        sequence_file = os.path.join(os.path.abspath(database_output_directory), "sequence_db.fasta")
        gbk_files = load_json_file(history_file_path)
        protein_hash_keys_to_gbk_hashes = defaultdict(set)
        protein_hash_keys_to_accessions = defaultdict(set)
        prot_keys = set()
        gbk = 0
        total = len(gbk_files)
        with open(sequence_file, 'w') as fasta_out:
            for gbk_hash, gbk_path in gbk_files.items():
                gbk += 1
                if gbk % 1 == 0:
                    print(f"\t > {current_time()} {gbk} of {total} processed ... ")
                if gbk % len(gbk_files) == 0:
                    pass
                protein_sequences = extract_protein_sequences(gbk_path)
                for accession, gene, prot_product, sequence, protein_hash in protein_sequences:
                    protein_hash_keys_to_accessions[protein_hash].add(accession)
                    protein_hash_keys_to_gbk_hashes[protein_hash].add(gbk_hash)
                    if protein_hash not in prot_keys:
                        # noinspection PyTypeChecker
                        print(seq_to_fasta(accession, gene, prot_product, protein_hash, sequence), file=fasta_out)
                        prot_keys.add(protein_hash)

        write_dict_of_set_to_json_as_file(protein_hash_keys_to_gbk_hashes, protein_hash_keys_to_gbk_hashes_file_path)
        write_dict_of_set_to_json_as_file(protein_hash_keys_to_accessions, protein_hash_keys_to_accessions_file_path)
        
        print(f"\nCompleted at {current_time()}\n")


if __name__ == "__main__":
    main()
