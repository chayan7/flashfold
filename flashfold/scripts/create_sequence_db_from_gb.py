# Author: Chayan Kumar Saha

import os
import sys
from collections import defaultdict
from flashfold.utils import extract_protein_sequences, is_valid_path, get_hash_to_files_with_extensions_from_dir, \
    current_time, seq_to_fasta, write_dict_of_set_to_json_as_file


genbank_file_extensions = [".gbff", ".gbk"]


def create_protein_db_from_gbk(database_output_directory: str, gbk_dir: str) -> None:

    gbk_files = get_hash_to_files_with_extensions_from_dir(gbk_dir, genbank_file_extensions)

    total = len(gbk_files)
    if total == 0:
        print(f"\n Error: Can not find genbank file with valid extensions eg. {genbank_file_extensions} in the "
              f"input directory: {gbk_dir} ")
        sys.exit()

    print(f"\n\tCreating database ... \n")

    if not is_valid_path(database_output_directory):
        os.makedirs(database_output_directory)

    print(f"\nStarted at {current_time()}")

    protein_hash_keys_to_gbk_hashes_file_path = os.path.join(os.path.abspath(database_output_directory),
                                                             "protein_to_gbks.json")
    protein_hash_keys_to_accessions_file_path = os.path.join(os.path.abspath(database_output_directory),
                                                             "prot_hash_to_accession.json")
    sequence_file = os.path.join(os.path.abspath(database_output_directory), "sequence_db.fasta")
    protein_hash_keys_to_gbk_hashes = defaultdict(set)
    protein_hash_keys_to_accessions = defaultdict(set)

    prot_keys = set()
    gbk = 0
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
