# Author: Chayan Kumar Saha

import os
import sys
from flashfold.utils import extract_protein_sequences, get_hash_to_files_with_extensions_from_dir, \
    current_time, seq_to_fasta, write_dict_of_set_to_json_as_file, is_valid_database_dir, \
    get_filename_to_path_set_by_directory, load_json_file, get_new_elements_from_second_list, \
    get_directory_from_file_path, SequenceDbFasta

genbank_file_extensions = [".gbff", ".gbk"]


def extend_main_sequence_db(main_db_path: str, new_genbank_path: str, new_db_path: str) -> None:
    if not is_valid_database_dir(main_db_path):
        sys.exit()
    protein_hash_key_to_gbk_hashes_file = os.path.join(os.path.abspath(main_db_path), "protein_to_gbks.json")
    protein_hash_key_to_accessions_file = os.path.join(os.path.abspath(main_db_path), "prot_hash_to_accession.json")
    main_sequence_file = os.path.join(os.path.abspath(main_db_path), "sequence_db.fasta")
    main_protein_hash_to_gbks = load_json_file(protein_hash_key_to_gbk_hashes_file)
    main_protein_hash_to_accessions = load_json_file(protein_hash_key_to_accessions_file)
    if new_genbank_path:
        gbk_files = get_hash_to_files_with_extensions_from_dir(new_genbank_path, genbank_file_extensions)
        total = len(gbk_files)
        if total == 0:
            print(f"\n Error: Can not find genbank file with valid extensions eg. {genbank_file_extensions} in the "
                  f"input directory: {new_genbank_path} ")
            sys.exit()
        print(f"\n-- Extending main database ... \n")
        print(f"\n\t-- Started at {current_time()}")
        gbk = 0
        for gbk_hash, gbk_path in gbk_files.items():
            protein_sequences = extract_protein_sequences(gbk_path)
            gbk += 1
            if gbk % 1 == 0:
                print(f"\t > {current_time()} updating {gbk} of {total} ... ")
            if gbk % len(gbk_files) == 0:
                print(f"\t > {current_time()} updating {gbk} of {total} ... ")
            for accession, gene, prot_product, sequence, protein_hash in protein_sequences:
                if protein_hash in main_protein_hash_to_gbks:
                    if gbk_hash not in main_protein_hash_to_gbks[protein_hash]:
                        main_protein_hash_to_gbks[protein_hash].append(gbk_hash)
                else:
                    main_protein_hash_to_gbks[protein_hash] = [gbk_hash]
                    main_protein_hash_to_accessions[protein_hash] = [accession]
                    with open(main_sequence_file, "a") as seq_file:
                        new_fasta = seq_to_fasta(accession, gene, prot_product, protein_hash, sequence)
                        seq_file.write(f"{new_fasta}\n")
    elif new_db_path:
        if not is_valid_database_dir(new_db_path):
            sys.exit()
        print(f"\n-- Extending main database ... \n")
        print(f"\n\t-- Started at {current_time()}")
        database_files = get_filename_to_path_set_by_directory(new_db_path, [".json"])
        protein_to_gbks_json_set = database_files["protein_to_gbks.json"]
        for protein_to_gbks_json_path in protein_to_gbks_json_set:
            protein_to_gbks = load_json_file(protein_to_gbks_json_path)
            sequence_file = os.path.join(get_directory_from_file_path(protein_to_gbks_json_path), "sequence_db.fasta")
            parsed_sequence_file = SequenceDbFasta(sequence_file)
            protein_to_accessions_file_path = os.path.join(get_directory_from_file_path(protein_to_gbks_json_path),
                                                           "prot_hash_to_accession.json")
            protein_to_accessions = load_json_file(protein_to_accessions_file_path)
            for protein_hash, gbk_hashes in protein_to_gbks.items():
                if protein_hash in main_protein_hash_to_gbks:
                    for gbk_hash in get_new_elements_from_second_list(main_protein_hash_to_gbks[protein_hash],
                                                                      protein_to_gbks[protein_hash]):

                        main_protein_hash_to_gbks[protein_hash].append(gbk_hash)
                else:
                    main_protein_hash_to_gbks[protein_hash] = protein_to_gbks[protein_hash]
                    main_protein_hash_to_accessions[protein_hash] = protein_to_accessions[protein_hash]
                    with open(main_sequence_file, "a") as seq_file:
                        new_fasta = parsed_sequence_file.get_fasta_by_protein_hash(protein_hash)
                        seq_file.write(f"{new_fasta}\n")

    write_dict_of_set_to_json_as_file(main_protein_hash_to_gbks, protein_hash_key_to_gbk_hashes_file)
    write_dict_of_set_to_json_as_file(main_protein_hash_to_accessions, protein_hash_key_to_accessions_file)

    print(f"\t-- Completed at {current_time()}\n")

    return None
