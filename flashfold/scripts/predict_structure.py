# Author: Chayan Kumar Saha

import os.path
import sys
from flashfold.utils import is_valid_protein_fasta, is_valid_database_dir, manage_output_path, \
    join_list_elements_by_character, create_new_directory, is_installed, run_jackhmmer, \
    get_files_from_path_by_extension, run_colabfold, current_time, current_time_raw, update_time_log, \
    get_valid_sequence_records_from_fasta, get_input_fasta_features, Database, create_a3m_for_folding, \
    get_query_to_a3m_records, generate_score_matrix


def predict_3d_structure(fasta_file: str, database_path: str, out_dir: str, cpu: int, num_models: int, num_recycle: int,
                         stop_at_score: int, num_top: int, relax_max_iterations: int,
                         overwrite_existing_results: bool, cutoff: float) -> None:

    prediction_start_time = current_time_raw()

    if not is_valid_protein_fasta(fasta_file) or not is_valid_database_dir(database_path):
        sys.exit()

    package_root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    input_fasta_records = get_valid_sequence_records_from_fasta(fasta_file)
    query_fasta_features = get_input_fasta_features(input_fasta_records)
    is_monomer = len(query_fasta_features.seqs) == 1
    parent_result_path = manage_output_path(out_dir, overwrite_existing_results)
    create_new_directory(parent_result_path)
    
    # Time_log
    time_log_file = os.path.join(parent_result_path, "timings.txt")
    initial_msg = f"-- FlashFold analysis began at: {current_time()}"
    update_time_log(time_log_file, initial_msg, False)

    # Load sequence database
    sequence_database = Database(database_path)

    # Make fasta files for homology search
    child_fasta_path = os.path.join(parent_result_path, "flashfold_fasta")
    create_new_directory(child_fasta_path)
    hash_to_fasta = query_fasta_features.hash_to_fasta
    unique_fasta_file_path = []
    for hash_id in hash_to_fasta:
        write_fasta = os.path.join(child_fasta_path, f"{hash_id}.fasta")
        with open(write_fasta, "w", encoding="utf-8") as fasta_file_out:
            # noinspection PyTypeChecker
            print(hash_to_fasta[hash_id], file=fasta_file_out)
        unique_fasta_file_path.append(write_fasta)

    # Run homology searching against flashfold databases
    jackhmmer_is_installed = is_installed("jackhmmer")
    if not jackhmmer_is_installed:
        print("Missing Jackhmmer in system path: it is recommended to reinstall and run flashfold afterwards.")
        sys.exit()

    perl_is_installed = is_installed("perl")
    if not perl_is_installed:
        print("Missing Perl in system path: it is recommended to reinstall and run flashfold afterwards.")
        sys.exit()
    path_reformat_script = os.path.join(package_root_dir, "extra", "reformat.pl")
    child_alignment_path = os.path.join(parent_result_path, "flashfold_msa")
    create_new_directory(child_alignment_path)

    run_jackhmmer(path_reformat_script, unique_fasta_file_path, sequence_database.fasta_db, cpu, child_alignment_path)

    # Process homology search output
    process_time = current_time()
    print(f"\n-- {process_time} > Selecting the best homologue sequences for structure prediction")
    alignment_path = os.path.join(parent_result_path, "flashfold_msa")
    homology_summary_json_file = os.path.join(alignment_path, "homology_summary.json")
    sequence_database.process_homology_search_output(alignment_path, query_fasta_features.chain_seq_hashes,
                                                     homology_summary_json_file)

    # Make combined alignment file
    a3m_files = get_files_from_path_by_extension(alignment_path, ".a3m")
    a3m_records = get_query_to_a3m_records(a3m_files, query_fasta_features.chain_seq_hashes)
    filtered_a3m_path = os.path.join(parent_result_path, "flashfold_filtered")
    create_new_directory(filtered_a3m_path)
    create_a3m_for_folding(homology_summary_json_file, a3m_records, query_fasta_features, filtered_a3m_path)
    update_time_log(time_log_file, "Completed Step 1: Alignment", True)

    # Run_colabfold
    colabfold_is_installed = is_installed("colabfold_batch")
    if not colabfold_is_installed:
        print("Missing colabfold_batch: it is recommended to reinstall and run flashfold afterwards.")
        sys.exit()

    a3m_file_name = join_list_elements_by_character(query_fasta_features.accnrs, "-")
    combined_homology_search_output = os.path.join(filtered_a3m_path, f"{a3m_file_name}.a3m")
    structure_prediction_out_path = os.path.join(parent_result_path, "flashfold_predicted_structure")
    create_new_directory(structure_prediction_out_path)
    run_colabfold(is_monomer, combined_homology_search_output, structure_prediction_out_path, num_models,
                  num_recycle, stop_at_score, num_top, relax_max_iterations)
    update_time_log(time_log_file, "Completed Step 2: Structure prediction", True)
    generate_score_matrix(structure_prediction_out_path, cutoff, is_monomer)
    update_time_log(time_log_file, "Completed Step 3: Scoring", True)
    prediction_end_time = current_time_raw()
    prediction_duration = prediction_end_time - prediction_start_time
    prediction_duration_seconds = prediction_duration.total_seconds()
    program_timing_msg = f"\n-- Total time in seconds: {prediction_duration_seconds}"
    update_time_log(time_log_file, program_timing_msg, False)
