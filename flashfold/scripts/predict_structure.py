# Author: Chayan Kumar Saha

import os.path
import shutil
import sys
from typing import List, Dict, Tuple, Literal
from collections import defaultdict
from collections import namedtuple


from flashfold.utils import is_valid_protein_fasta, is_valid_database_dir, manage_output_path, \
    join_list_elements_by_character, create_new_directory, run_jackhmmer, \
    get_files_from_path_by_extension, run_colabfold, current_time, current_time_raw, update_time_log, \
    get_valid_sequence_records_from_fasta, get_input_fasta_features, Database, create_a3m_for_folding, \
    get_query_to_a3m_records, generate_score_matrix, get_filename_without_extension, replace_char_from_string, \
    JsonStructure, run_jobs_in_parallel, is_valid_protein_a3m, is_a3m_monomer

# Define the named tuple
Valid_Input = namedtuple('Valid_Input', ['file_ext', 'file_paths'])


def get_valid_input_files_with_type(dir_or_file: str, is_batch: bool) -> Valid_Input:
    """
    Retrieves a list of valid input files with extension, from the specified directory or file.

    Parameters:
    dir_or_file (str): Path to a directory or a single file.
    is_batch (bool): Flag indicating whether the input is a directory (batch mode) or a single file.

    Returns:
    List[str]: A list of paths to valid files and the file extension.

    Raises:
    SystemExit: If the input path is invalid.
    """
    file_paths: List[str] = []
    file_ext: Literal["fasta", "a3m"] = "fasta"
    if is_batch:
        if not os.path.isdir(dir_or_file):
            print(f"\n-- Error: For --batch option, the input query should be a directory with valid FASTA/A3M files. "
                  f"The input query is not a directory. Please try again with a valid path.\n")
            sys.exit()
        else:
            files_with_fasta_ext = get_files_from_path_by_extension(dir_or_file, ".fasta")
            files_with_a3m_ext = get_files_from_path_by_extension(dir_or_file, ".a3m")
            if len(files_with_fasta_ext) >= 1 and len(files_with_a3m_ext) >= 1:
                print(f"\n-- Error: Provided path contains both .fasta and .a3m files. Please use either FASTA or A3M.")
                sys.exit()
            if len(files_with_fasta_ext) >= 1:
                file_paths.extend(files_with_fasta_ext)
            elif len(files_with_a3m_ext) >= 1:
                file_ext = "a3m"
                file_paths.extend(files_with_a3m_ext)
            else:
                print(f"\n-- Error: No files detected with .fasta or .a3m extension.")
                sys.exit()
    else:
        if os.path.isfile(dir_or_file):
            if dir_or_file.endswith(".fasta"):
                file_paths.append(dir_or_file)
            elif dir_or_file.endswith(".a3m"):
                file_ext = "a3m"
                file_paths.append(dir_or_file)
            else:
                print(f"\n-- Error: No files detected with .fasta or .a3m extension.")
                sys.exit()
        else:
            print(f"\n-- Error: The input query is not a file. "
                  f"Please provide a valid FASTA/A3M file path.\n")
            sys.exit()

    # checking if file is valid
    invalid_files: List[str] = []
    if file_ext == "fasta":
        for fasta_file in file_paths:
            if not is_valid_protein_fasta(fasta_file):
                invalid_files.append(fasta_file)
    elif file_ext == "a3m":
        for a3m_file in file_paths:
            if not is_valid_protein_a3m(a3m_file):
                invalid_files.append(a3m_file)

    if len(invalid_files) != 0:
        format_invalids = join_list_elements_by_character(invalid_files, "\n")
        print(f"\n-- Error: The following files are not valid protein {file_ext.upper()} files: \n{format_invalids}\n"
              f"Please provide valid protein FASTA files and try again.\n")
        sys.exit()

    valid_files = Valid_Input(file_ext, file_paths)
    return valid_files


def get_results_and_temp_dirs(out_dir_path: str, rewrite: bool, is_temp_needed: bool) -> Tuple[str, str]:
    """
    Creates the output and temporary directories for the results.

    Parameters:
    out_dir_path (str): Path to the output directory.
    is_batch (bool): Flag indicating whether the input is a directory (batch mode) or a single file.
    rewrite (bool): Flag indicating whether to overwrite existing results.

    Returns:
    str: Path to the parent directory for the results.

    Raises:
    SystemExit: If the output directory path is invalid.
    """
    parent_result_path = manage_output_path(out_dir_path, rewrite)
    create_new_directory(parent_result_path)

    if not is_temp_needed:
        return parent_result_path, "NA"

    # Create a temporary directory
    temp_dir_path = os.path.join(parent_result_path, "temp")
    create_new_directory(temp_dir_path)

    return parent_result_path, temp_dir_path


def predict_3d_structure(args) -> None:
    # Time_log for the entire program
    prediction_start_time = current_time_raw()

    query = args.query
    is_batch = args.batch

    # Get input FASTA files
    valid_input_files = get_valid_input_files_with_type(query, is_batch)
    is_fasta = valid_input_files.file_ext == "fasta"

    if args.database and not is_fasta:
        print(f"\n-- Warning: Database path is only required and used when the query file is in FASTA format.")
        print(f"-- Tip: Skip using the database when query file is an MSA (.a3m).\n")

    if args.only_msa and not is_fasta:
        print(f"\n-- Error: Provided input file type is already MSA and --only_msa is set true")
        print(f"-- Tip: Either use FASTA file(s) as input or skip --only_msa for structure prediction :) ...\n")
        sys.exit()

    # Create output directory
    out_dir_path = os.path.realpath(args.output)
    overwrite = args.overwrite_existing_results
    parent_result_path, temp_dir_path = get_results_and_temp_dirs(out_dir_path, overwrite, is_fasta)
    
    time_log_file = os.path.join(out_dir_path, "batch_timings.txt") if is_batch \
        else os.path.join(parent_result_path, "timings.txt")
    log_text = "-batch" if is_batch else ""
    
    initial_msg = f"-- FlashFold{log_text} analysis began at: {prediction_start_time}"
    update_time_log(time_log_file, initial_msg, False)

    # Load the database if needed
    sequence_database = None
    if is_fasta:
        if not args.database:
            print(f"\n-- Error: Database path is required when the query file is in FASTA format.\n")
            shutil.rmtree(out_dir_path)
            sys.exit()

        database_path = os.path.realpath(args.database)

        # Check if the database directory is valid
        if not is_valid_database_dir(database_path):
            shutil.rmtree(out_dir_path)
            sys.exit()

        # Load database
        sequence_database = Database(database_path)

    # Create a JsonStructure object to store the features of the input FASTA files
    infile_features_json = JsonStructure()

    # Create a JsonStructure object to store the features required for structure prediction
    fold_features = JsonStructure()

    # Homology search Initialization
    query_hash_to_single_fasta_path = {}
    query_hash_to_msa_paths: Dict[str, List[str]] = defaultdict(list)

    for valid_file in valid_input_files.file_paths:
        valid_file_name = get_filename_without_extension(valid_file)
        sub_out_dir_path = replace_char_from_string(valid_file_name, "_")
        batch_out_dir_path = os.path.join(out_dir_path, sub_out_dir_path)

        if is_batch:
            parent_result_path = manage_output_path(batch_out_dir_path, overwrite)

        infile_features_json.add_entry(valid_file, "temp", temp_dir_path)

        #structure prediction out path
        structure_prediction_out_path = os.path.join(parent_result_path, "flashfold_structure")
        fold_features.add_entry(valid_file, "structure_prediction_path", structure_prediction_out_path)

        if not is_fasta:
            valid_a3m_file = valid_file
            is_monomer = is_a3m_monomer(valid_a3m_file)
            fold_features.add_entry(valid_a3m_file, "is_monomer", is_monomer)
            fold_features.add_entry(valid_a3m_file, "time_log", time_log_file)
            fold_features.add_entry(valid_a3m_file, "input_msa", valid_a3m_file)

        if is_fasta:
            valid_fasta_file = valid_file
            input_fasta_records = get_valid_sequence_records_from_fasta(valid_fasta_file)
            query_fasta_features = get_input_fasta_features(input_fasta_records)
            is_monomer = len(query_fasta_features.seqs) == 1

            infile_features_json.add_entry(valid_fasta_file, "is_monomer", is_monomer)
            infile_features_json.add_entry(valid_fasta_file, "features", query_fasta_features)
            infile_features_json.add_entry(valid_fasta_file, "parent_result_path", parent_result_path)

            infile_features_json.add_entry(valid_fasta_file, "time_log", time_log_file)

            # Make fasta files for homology search
            child_alignment_path = os.path.join(parent_result_path, "flashfold_msa")
            create_new_directory(child_alignment_path)
            infile_features_json.add_entry(valid_fasta_file, "msa", child_alignment_path)

            hash_to_fasta = query_fasta_features.hash_to_fasta
            for hash_id in hash_to_fasta:
                write_fasta = os.path.join(temp_dir_path, f"{hash_id}.fasta")
                query_hash_to_single_fasta_path[hash_id] = write_fasta
                query_hash_to_msa_paths[hash_id].append(child_alignment_path)
                with open(write_fasta, "w", encoding="utf-8") as fasta_file_out:
                    # noinspection PyTypeChecker
                    print(hash_to_fasta[hash_id], file=fasta_file_out)

    if is_fasta:
        # Run Jackhmmer
        unique_fasta_file_paths = list(query_hash_to_single_fasta_path.values())

        msa_start_log_text = f"Started Step 1: MSA construction for {len(valid_input_files.file_paths)} sequences" \
            if is_batch else f"Started Step 1: MSA construction"

        update_time_log(time_log_file, msa_start_log_text, True)
        run_jackhmmer(unique_fasta_file_paths, sequence_database.fasta_db, args.threads, temp_dir_path)

        # copy homology search output
        copy_alignment_files_commands = []
        for query_hash in query_hash_to_msa_paths:
            for msa_path in query_hash_to_msa_paths[query_hash]:
                copy_a3m_command = f"cp {temp_dir_path}/{query_hash}.a3m {msa_path}"
                copy_sto_command = f"cp {temp_dir_path}/{query_hash}.sto {msa_path}"
                copy_alignment_files_commands.append(copy_a3m_command)
                copy_alignment_files_commands.append(copy_sto_command)

        log_text = "files" if len(copy_alignment_files_commands) > 1 else "file"
        run_jobs_in_parallel(args.threads, 1, copy_alignment_files_commands,
                             f"Copying alignment {log_text}")

        # Process homology search output
        for each_fasta in infile_features_json.get_data():
            print(f"\n-- {current_time()} > Filtering MSA")
            each_fasta_features = infile_features_json.get_data()[each_fasta]
            is_monomer = each_fasta_features["is_monomer"]

            alignment_path = each_fasta_features["msa"]
            result_subdirectory = each_fasta_features["parent_result_path"]
            query_fasta_features = each_fasta_features["features"]
            time_log_file = each_fasta_features["time_log"]
            homology_summary_json_file = os.path.join(alignment_path, "homology_summary.json")
            sequence_database.process_homology_search_output(alignment_path, query_fasta_features.chain_seq_hashes,
                                                             homology_summary_json_file)

            # Make combined alignment file
            a3m_files = get_files_from_path_by_extension(alignment_path, ".a3m")
            a3m_records = get_query_to_a3m_records(a3m_files, query_fasta_features.chain_seq_hashes)
            filtered_a3m_path = os.path.join(result_subdirectory, "flashfold_filtered_msa")
            create_new_directory(filtered_a3m_path)
            create_a3m_for_folding(homology_summary_json_file, a3m_records, query_fasta_features, filtered_a3m_path)
            a3m_file_name = join_list_elements_by_character(query_fasta_features.accnrs, "-")
            filtered_homology_search_output = os.path.join(filtered_a3m_path, f"{a3m_file_name}.a3m")

            msa_end_log_text = f"Completed Step 1: MSA construction for {os.path.basename(result_subdirectory)} " \
                if is_batch else f"Completed Step 1: MSA construction"
            update_time_log(time_log_file, msa_end_log_text, True)

            structure_prediction_out_path = os.path.join(result_subdirectory, "flashfold_structure")
            fold_features.add_entry(each_fasta, "filtered_msa", filtered_homology_search_output)
            fold_features.add_entry(each_fasta, "structure_prediction_path", structure_prediction_out_path)
            fold_features.add_entry(each_fasta, "is_monomer", is_monomer)
            fold_features.add_entry(each_fasta, "time_log", time_log_file)
            #remove child msa path that contains msa per chain
            shutil.rmtree(alignment_path)

        # remove temp directory
        shutil.rmtree(temp_dir_path)

    if not args.only_msa:
        # Run structure prediction
        for each_file in fold_features.get_data():
            each_file_fold_features = fold_features.get_data()[each_file]
            is_monomer = each_file_fold_features["is_monomer"]
            filtered_a3m_file = each_file_fold_features["filtered_msa"] if is_fasta \
                else each_file_fold_features["input_msa"]
            prediction_out_path = each_file_fold_features["structure_prediction_path"]
            time_log_file = each_file_fold_features["time_log"]
            create_new_directory(prediction_out_path)

            fold_log_extra = f" for {os.path.basename(os.path.dirname(prediction_out_path))}" if is_batch else ""

            fold_start_log = f"Started Step 2: Structure prediction{fold_log_extra}"
            update_time_log(time_log_file, fold_start_log, True)
            run_colabfold(is_monomer, filtered_a3m_file, prediction_out_path, args.num_models,
                          args.num_recycle, args.stop_at_score, args.num_structure_relax, args.relax_max_iterations)

            fold_end_log = f"Completed Step 2: Structure prediction{fold_log_extra}"
            update_time_log(time_log_file, fold_end_log, True)

            score_start_log = f"Started Step 3: Scoring{fold_log_extra}"
            update_time_log(time_log_file, score_start_log, True)
            generate_score_matrix(prediction_out_path, args.cutoff, is_monomer)
            score_end_log = f"Completed Step 3: Scoring{fold_log_extra}"
            update_time_log(time_log_file, score_end_log, True)

    if not is_batch:
        prediction_end_time = current_time_raw()
        prediction_duration = prediction_end_time - prediction_start_time
        prediction_duration_seconds = prediction_duration.total_seconds()
        program_timing_msg = f"\n-- Total time in seconds: {prediction_duration_seconds}"
        update_time_log(time_log_file, program_timing_msg, False)

