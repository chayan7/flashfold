import os
import shutil
import subprocess
from typing import List, Dict, Optional
from pathlib import Path
from flashfold.utils import run_jobs_in_parallel, run_single_job, load_json_file, write_dict_to_json_as_file, \
    get_files_from_path_by_extension, current_time_raw, current_time, manage_output_path, update_time_log
from flashfold.tools import run_alphafold3_docker


def get_af3_configurations(af3_config_file_path: str, af3_image: str, af3_db: Optional[Path],
                           af3_params: Optional[Path]) -> Optional[Dict[str, str]]:
    """
    Get the AlphaFold3 configurations from the configuration file or command line arguments.
    Args:
        af3_config_file_path: A path to the AlphaFold3 configuration file.
        af3_image: Name of docker image for AlphaFold3.
        af3_db: Af3_database path.
        af3_params: Af3_parameters path.

    Returns:
        A dictionary of AlphaFold3 configurations or None if the configurations are invalid.

    """

    old_config = {}
    new_config = {}

    if os.path.exists(af3_config_file_path):
        old_config.update(load_json_file(af3_config_file_path))

    # Check if the AlphaFold3 docker image is valid
    if af3_image != '' and af3_image != old_config.get('image', ''):
        af3_docker_check_com = f"docker run {af3_image} python run_alphafold.py --help"
        af3_check = subprocess.run(af3_docker_check_com, shell=True, check=False,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if not af3_check.stdout.decode().startswith("AlphaFold 3"):
            print(f"\t-- Error: The AlphaFold3 docker image is invalid.\n")
            return None

        new_config['image'] = af3_image


    # Check if the AlphaFold3 database is valid
    if af3_db:
        af3_db_str = str(os.path.realpath(af3_db))
        if af3_db_str != old_config.get('db', ''):
            if af3_db.is_file():
                print(f"\t-- Error: The path to AlphaFold3 database is not valid.\n")
                return None
            else:
                new_config['db'] = af3_db_str

    # Check if the AlphaFold3 database and parameters are valid
    if af3_params:
        af3_params_str = str(os.path.realpath(af3_params))
        if af3_params_str != old_config.get('params', ''):
            if af3_params.is_file():
                print(f"\t-- Error: The path to AlphaFold3 database is not valid.\n")
                return None
            else:
                new_config['params'] = af3_params_str

    is_new_config_same_as_old = old_config == new_config

    if len(new_config) == 3 and not is_new_config_same_as_old:
        # Save the configuration to a file
        print(f"\t-- Update: Saving new AlphaFold3 configurations to: \n\t'{af3_config_file_path}'\n")
        write_dict_to_json_as_file(new_config, af3_config_file_path)
        return new_config

    return old_config


def get_query_files(query_path: str, is_batch: bool) -> Optional[List[str]]:

    query_files = []
    if is_batch:
        if not os.path.exists(query_path) or not os.path.isdir(query_path):
            print(f"\t-- Error: The path to the query directory is not valid. "
                  f"Or, remove '--batch' if provided path is a single query. \n")
            return None

        potential_query_files = get_files_from_path_by_extension(query_path, "json")

        if len(potential_query_files) == 0:
            print(f"\t-- Error: No JSON files detected in the query directory.\n")
            return None

        return potential_query_files
    else:
        if not os.path.exists(query_path) or not os.path.isfile(query_path):
            print(f"\t-- Error: The path to the query file is not valid. "
                  f"Or, use '--batch' if provided path is a directory.\n")
            return None

        query_files.append(query_path)

    return query_files


def process_alphafold3_task(args) -> None:

    prediction_start_time = current_time_raw()

    # Check the AlphaFold3 configurations
    print(f"\n-- {current_time()} > Checking configuration for AlphaFold3 ...\n")

    current_working_directory = os.path.dirname(os.path.abspath(__file__))
    af3_config_file_path = os.path.join(current_working_directory, "af3_config.json")

    config = get_af3_configurations(af3_config_file_path, args.af3_image, args.af3_db, args.af3_params)

    if not config:
        return

    print(f"-- {current_time()} > Checking configuration for AlphaFold3 is complete. \n")


    # Check the query files
    print(f"\n-- {current_time()} > Searching JSON file(s) for AlphaFold3 ...\n")

    is_batch = args.batch
    query_files = get_query_files(args.query, is_batch)

    if not query_files or len(query_files) == 0:
        return

    print(f"-- {current_time()} > Searching JSON file(s) for AlphaFold3 is complete [Found= {len(query_files)}]. \n")


    # Check and create the output directory if it does not exist
    overwrite = args.overwrite_existing_results
    output_dir = manage_output_path(args.output, overwrite)

    # create a time log file
    time_log_file = os.path.join(output_dir, "timings.txt")
    initial_msg = (f"-- AlphaFold3 structure prediction (with FlashFold generated JSON) analysis began at: "
                   f"{prediction_start_time}")
    update_time_log(time_log_file, initial_msg, False)

    # Run the AlphaFold3 task(s)
    run_alphafold3_docker(config, query_files, output_dir, time_log_file)

    prediction_end_time = current_time_raw()
    prediction_duration = prediction_end_time - prediction_start_time
    prediction_duration_seconds = prediction_duration.total_seconds()
    program_timing_msg = f"\n-- Total time in seconds: {prediction_duration_seconds}"
    update_time_log(time_log_file, program_timing_msg, False)
    print(f"-- {current_time()} > AlphaFold3 task(s) has been completed. Check output at: '{output_dir}'.\n")

    return