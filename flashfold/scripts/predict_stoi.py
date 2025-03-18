import os
import itertools
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Union

from tensorboard import summary

from flashfold.utils import get_valid_sequence_records_from_fasta, get_input_fasta_features, is_valid_protein_fasta, \
    manage_output_path, run_single_job, current_time_raw, update_time_log, get_files_from_path_by_extension
from .run_alphafold3 import get_af3_configurations, run_alphafold3_docker
from .generate_report import get_summary_table_rows_from_result_path, filter_columns, generate_html_table, \
    generate_csv_table, remove_na_columns


def number_to_alphabet(num: int) -> str:
    if num < 1 or num > 26:
        raise ValueError("Number must be between 1 and 26")
    return chr(num + 96).upper()


def modify_query_to_stoichiometry(query_name: str) -> str:
    query_name_split = query_name.split("_")
    stoichiometry_num = []

    for i in range(len(query_name_split)):
        if i%2 != 0:
            stoichiometry_num.append(query_name_split[i])

    stoichiometry = []
    for i in range(len(stoichiometry_num)):
        stoichiometry_str = f"{number_to_alphabet(i+1)}{stoichiometry_num[i]}"
        stoichiometry.append(stoichiometry_str)

    return "".join(stoichiometry)


def find_stoichiometries_with_highest_value(data: List[List[str]], column_index: int) -> Optional[List[str]]:
    highest_value = float('-inf')
    value_to_stoichiometries = defaultdict(list)
    for row in data:
        try:
            value = float(row[column_index])
            if value >= highest_value:
                highest_value = value
                stoichiometry = row[0]
                value_to_stoichiometries[value].append(stoichiometry)
        except ValueError:
            continue

    stoichiometries = value_to_stoichiometries[highest_value]

    if len(stoichiometries) == 0:
        return None

    return stoichiometries


def change_query_row_value(data: List[List[str]], column_index: int = 0) -> None:
    for row in data:
        try:
            value = row[column_index]
            row[column_index] = modify_query_to_stoichiometry(value)
        except ValueError:
            continue


def predict_stoichiometry(args) -> None:

    query_file_path = args.query
    if not is_valid_protein_fasta(query_file_path):
        return

    current_working_directory = os.path.dirname(os.path.abspath(__file__))
    af3_config_file_path = os.path.join(current_working_directory, "af3_config.json")

    use_af3 = args.use_af3

    config: Optional[Dict[str, str]] = None
    if use_af3:
        config = get_af3_configurations(af3_config_file_path, args.af3_image, args.af3_db, args.af3_params)
        if not config:
            return

    valid_fasta_file = os.path.realpath(query_file_path)
    input_fasta_records = get_valid_sequence_records_from_fasta(valid_fasta_file)
    query_fasta_features = get_input_fasta_features(input_fasta_records)


    total_num_of_unique_chains = len(query_fasta_features.chain_seq_hashes)
    all_chains_in_query = len(query_fasta_features.seqs)

    if total_num_of_unique_chains == 0 or all_chains_in_query != total_num_of_unique_chains:
        print(f"\n-- Error: Invalid query file. \n")
        return

    unique_chain_ids = query_fasta_features.accnrs

    for chain_id in unique_chain_ids:
        if chain_id.count("_") > 0:
            print(f"\n-- Error: Invalid chain name detected in the query file. \n"
                  f"-- Tip: Chain name should not contain underscore. \n")
            return

    is_single_chain = total_num_of_unique_chains == 1

    stoi_dict: Dict[str, int] = {}

    specific_stoichiometry: Optional[List[List[str]]] = args.specific_stoichiometry
    if specific_stoichiometry:
        for chain_id, copy_num in specific_stoichiometry:
            if chain_id in stoi_dict:
                print(f"\n-- Error: Please provide a unique chain name while using '--ss' option. "
                      f"Chain {chain_id} is repeated.\n")
                return

            if chain_id not in unique_chain_ids:
                print(f"\n-- Error: Chain {chain_id} is not found in the query file. \n")
                return

            stoi_dict[chain_id] = int(copy_num)

        if len(stoi_dict) != total_num_of_unique_chains:
            print(f"\n-- Error: Invalid number of specific stoichiometry input has been detected. \n"
                  f"-- Tip: Chain number detected in the query file ({total_num_of_unique_chains}) is not equal to "
                  f"the number assigned in the specific stoichiometry (-ss) input ({len(stoi_dict)}). \n")
            return

        if len(stoi_dict) == 1 and len(stoi_dict) == total_num_of_unique_chains:
            if not sum(list([int(stoi_dict[chain]) for chain in stoi_dict])) >= 2:
                print(f"\n-- Error: Invalid number has been assigned for specific stoichiometry (-ss)."
                      f"\n-- Tip: For the provided query file, the copy_number should be at least 2.\n")
                return


    global_stoichiometry: Optional[int] = args.global_stoichiometry
    if global_stoichiometry:
        if global_stoichiometry < 2:
            print(f"\n-- Error: Invalid number has been assigned for global stoichiometry (-gs) . "
                  f"Minimum requirement is >=2 \n")
            return
        for chain_id in unique_chain_ids:
            stoi_dict[chain_id] = global_stoichiometry

    if len(stoi_dict) == 0:
        print(f"\n-- Error: No stoichiometry input has been detected. \n"
              f"-- Tip: Please provide either specific stoichiometry (-ss) or global stoichiometry (-gs) input. \n")
        return

    list_of_chain_copy_num_list = []
    for key in stoi_dict:
        chain_copy_num_list = []
        for i in range(1, stoi_dict[key] + 1):
            stoi_str = f"{key}_{i}"
            chain_copy_num_list.append(stoi_str)
        list_of_chain_copy_num_list.append(chain_copy_num_list)


    stoi_to_analyze = list(itertools.product(*list_of_chain_copy_num_list))

    output_dir = manage_output_path(args.output, args.overwrite_existing_results)

    intermediate_file_dir = os.path.join(output_dir, "intermediate_files")
    os.makedirs(intermediate_file_dir)

    fasta_dir = os.path.join(intermediate_file_dir, "fasta_files")
    os.makedirs(fasta_dir)

    msa_dir = os.path.join(intermediate_file_dir, "msa_files")
    predicted_model_dir = os.path.join(intermediate_file_dir, "model_files")

    for chains in stoi_to_analyze:
        list_of_chains = list(chains)
        fasta_file_name = "_".join(list_of_chains) + ".fasta"
        fasta_file_path = os.path.join(fasta_dir, fasta_file_name)
        with open(fasta_file_path, "w") as fasta_file:
            for chain in list_of_chains:
                chain_id = "_".join(chain.split("_")[0:-1])
                chain_num = int(chain.split("_")[-1])
                chain_seq = query_fasta_features.seqs[unique_chain_ids.index(chain_id)]
                for i in range(1, chain_num+1):
                    fasta_file.write(f">{chain_id}_{i}\n{chain_seq}\n")


    # run flashfold msa prediction

    msa_command = f"flashfold fold -q {fasta_dir} -d {args.database} -o {msa_dir} -t {args.threads} --batch --only_msa"
    run_single_job(msa_command, "Building MSA")


    # structure prediction
    if use_af3:

        json_file_path = os.path.join(intermediate_file_dir, "json_files")
        json_command = f"flashfold fold -q {msa_dir} -o {json_file_path} --batch --only_json"
        run_single_job(json_command, "Generating JSON files")

        query_files = get_files_from_path_by_extension(json_file_path, "json")

        # Run the AlphaFold3 task(s)
        if not config:
            return

        os.makedirs(predicted_model_dir)

        run_alphafold3_docker(config, query_files, predicted_model_dir, None)

        summary_table_headers, summary_table_rows = get_summary_table_rows_from_result_path(predicted_model_dir,
                                                                                            True)
        filtered_table_headers,  filtered_table_rows = filter_columns(summary_table_headers, summary_table_rows,
                                                                      [0, 3, 4, 5, 7, 11, 12, 14])

    else:
        model_command = f"flashfold fold -q {msa_dir} -o {predicted_model_dir} --batch --calc_extra_ptm"
        run_single_job(model_command, "Predicting structures")

        summary_table_headers, summary_table_rows = get_summary_table_rows_from_result_path(predicted_model_dir,
                                                                                            False)
        filtered_table_headers, filtered_table_rows = filter_columns(summary_table_headers, summary_table_rows,
                                                                     [0, 3, 4, 5, 7, 11, 12, 14])

    if not filtered_table_rows or len(filtered_table_rows) == 0:
        print(f"\n-- Error: No results have been generated. Please try again.\n")
        return

    change_query_row_value(filtered_table_rows)

    clean_sum_tab_headers, clean_sum_tab_rows = remove_na_columns(filtered_table_headers, filtered_table_rows)

    if is_single_chain:
        ptm_index = 2
        stoi_with_higest_score = find_stoichiometries_with_highest_value(filtered_table_rows, ptm_index)
        score_criteria = filtered_table_headers[ptm_index]
    else:
        iptm_index = 3
        stoi_with_higest_score = find_stoichiometries_with_highest_value(filtered_table_rows, iptm_index)
        score_criteria = filtered_table_headers[iptm_index]

    if not stoi_with_higest_score or not score_criteria:
        print(f"\n-- Error: No results have been generated. Please try again.\n")
        return

    str_stoichiometry = ", ".join(stoi_with_higest_score)
    footer_text = f"Based on {score_criteria} scoring, predicted best stoichiometry: {str_stoichiometry}"

    print(f"\n-- Result: ")
    print(f"-- {footer_text}")

    output_html_file_path = os.path.join(output_dir, 'summary.html')
    output_csv_file_path = os.path.join(output_dir, 'summary.csv')
    generate_html_table(clean_sum_tab_headers, clean_sum_tab_rows, output_html_file_path, footer_text)
    generate_csv_table(clean_sum_tab_headers, clean_sum_tab_rows, output_csv_file_path)
    print(f"-- HTML report has been generated: {output_html_file_path}\n")

    return