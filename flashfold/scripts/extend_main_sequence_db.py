import os
import sys
import concurrent.futures
from flashfold.utils import extract_protein_sequences, get_hash_to_files_with_extensions_from_dir, load_json_file, \
    current_time, create_fasta_for_db, SequenceDbFasta, is_valid_database_dir, write_dict_of_set_to_json_as_file

genbank_file_extensions = [".gbff", ".gbk"]


def process_gbk_file(gbk_hash, gbk_path, main_protein_hash_to_gbks, main_protein_hash_to_accessions):
    new_content = 0
    count_new_prot = 0
    new_fasta_entries = []
    new_gbk_set = set()
    protein_sequences = extract_protein_sequences(gbk_path)
    for accession, gene, prot_product, sequence, protein_hash in protein_sequences:
        if protein_hash in main_protein_hash_to_gbks:
            if gbk_hash not in main_protein_hash_to_gbks[protein_hash]:
                new_content += 1
                new_gbk_set.add(gbk_hash)
                main_protein_hash_to_gbks[protein_hash].append(gbk_hash)
        else:
            new_content += 1
            count_new_prot += 1
            new_gbk_set.add(gbk_hash)
            main_protein_hash_to_gbks[protein_hash] = [gbk_hash]
            main_protein_hash_to_accessions[protein_hash] = [accession]
            new_fasta = create_fasta_for_db(accession, gene, prot_product, protein_hash, sequence)
            new_fasta_entries.append(new_fasta)
    return new_content, count_new_prot, new_gbk_set, new_fasta_entries


def extend_main_sequence_db(args) -> None:
    if not is_valid_database_dir(args.main_db):
        sys.exit()

    main_db_path = os.path.realpath(args.main_db)

    protein_hash_key_to_gbk_hashes_file = os.path.join(main_db_path, "protein_to_gbks.json")
    protein_hash_key_to_accessions_file = os.path.join(main_db_path, "prot_hash_to_accession.json")
    main_sequence_file = os.path.join(os.path.abspath(main_db_path), "sequence_db.fasta")
    main_protein_hash_to_gbks = load_json_file(protein_hash_key_to_gbk_hashes_file)
    main_protein_hash_to_accessions = load_json_file(protein_hash_key_to_accessions_file)

    new_content = 0
    count_new_prot = 0
    new_gbk_set = set()
    new_fasta_entries = []

    if args.genbank_path:
        new_genbank_path = os.path.realpath(args.genbank_path)
        gbk_files = get_hash_to_files_with_extensions_from_dir(new_genbank_path, genbank_file_extensions)
        total = len(gbk_files)
        if total == 0:
            print(f"\n Error: Can not find genbank file with valid extensions eg. {genbank_file_extensions} in the "
                  f"input directory: {new_genbank_path} ")
            sys.exit()
        print(f"\n-- {current_time()} - Cross-checking new information with main database: {main_db_path} ...")

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(process_gbk_file, gbk_hash, gbk_path, main_protein_hash_to_gbks,
                                       main_protein_hash_to_accessions) for gbk_hash, gbk_path in gbk_files.items()]
            total_files = len(futures)
            for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
                result = future.result()
                new_content += result[0]
                count_new_prot += result[1]
                new_gbk_set.update(result[2])
                new_fasta_entries.extend(result[3])
                if i % 10 == 0:
                    print(f"\t-- Processed {i}/{total_files} new GenBank files ...")

    elif args.new_db:
        new_db_path = os.path.realpath(args.new_db)

        if new_db_path == main_db_path:
            print(f"\nError: The new database path is same as the main database path. \n")
            sys.exit()

        new_protein_hash_key_to_gbks_file_path = os.path.join(new_db_path, "protein_to_gbks.json")
        new_protein_hash_key_to_accessions_file_path = os.path.join(new_db_path, "prot_hash_to_accession.json")
        new_db_sequence_file = os.path.join(os.path.abspath(new_db_path), "sequence_db.fasta")

        new_db_protein_hash_to_gbks = load_json_file(new_protein_hash_key_to_gbks_file_path)
        new_db_protein_hash_to_accessions = load_json_file(new_protein_hash_key_to_accessions_file_path)
        new_db_parsed_sequence_file = SequenceDbFasta(new_db_sequence_file)

        print(f"\n-- {current_time()} - Cross-checking new information with main database: {main_db_path} ...")

        i = 0
        for protein_hash, gbk_hashes in new_db_protein_hash_to_gbks.items():
            i += 1
            if i % 10 == 0:
                print(f"\t-- Processed {i}/{len(new_db_protein_hash_to_gbks)} entries from new database ...")
            if protein_hash in main_protein_hash_to_gbks:
                for gbk_hash in list(set(new_db_protein_hash_to_gbks[protein_hash]) -
                                     set(main_protein_hash_to_gbks[protein_hash])):
                    new_content += 1
                    new_gbk_set.update(gbk_hash)
                    main_protein_hash_to_gbks[protein_hash].append(gbk_hash)
            else:
                new_gbk_contents = new_db_protein_hash_to_gbks[protein_hash]
                new_content += 1
                new_gbk_set.update(new_gbk_contents)
                count_new_prot += 1
                main_protein_hash_to_gbks[protein_hash] = new_gbk_contents
                main_protein_hash_to_accessions[protein_hash] = new_db_protein_hash_to_accessions[protein_hash]
                new_fasta_record = new_db_parsed_sequence_file.get_record_by_protein_hash(protein_hash)
                new_fasta_entries.append(new_fasta_record.fasta)

    if new_content == 0:
        print(f"\t-- No new content found. \n")
        print(f"\n-- Completed at {current_time()}\n")
        sys.exit()

    print(f"\n ∞ {current_time()} ∞ Found total {new_content} new content(s) [protein: {count_new_prot}, "
          f"gbk: {len(new_gbk_set)}]\n")

    if not args.yes:
        print(f"\n-- You are about to extend the main database {main_db_path} with new information. \n"
              f"Are you sure? [y/n] ", end="")
        user_input = input()
        if user_input.lower() == "n":
            print(f"\n-- Exiting the extension of the main database {main_db_path} \n")
            sys.exit()

    with open(main_sequence_file, "a") as seq_file:
        for entry in new_fasta_entries:
            seq_file.write(f"{entry}\n")

    write_dict_of_set_to_json_as_file(main_protein_hash_to_gbks, protein_hash_key_to_gbk_hashes_file)
    write_dict_of_set_to_json_as_file(main_protein_hash_to_accessions, protein_hash_key_to_accessions_file)

    print(f"\n-- Completed at {current_time()}\n")
    return
