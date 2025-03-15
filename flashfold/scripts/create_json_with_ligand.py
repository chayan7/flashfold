import os
import sys
import shutil
from typing import List, Dict, Optional, Literal
from collections import namedtuple

from flashfold.tools import run_msa_to_json, modjson
from flashfold.utils import is_valid_protein_a3m, get_files_from_path_by_extension, current_time, \
                             join_list_elements_by_character, manage_output_path, is_valid_path, load_json_file


Valid_Ligand_Input = namedtuple('Valid_Ligand_Input', ['file_ext', 'file_paths'])


def check_userccd_path(userccd_file_paths: Optional[List[str]]) -> Optional[List[str]]:
    if userccd_file_paths is None:
        return None

    valid_userccd_files = []

    for userccd_file in userccd_file_paths:
        if not os.path.isfile(userccd_file):
            print(f"\n-- Error: Invalid path: '{userccd_file}'. Please provide a valid path to a userCCD file.\n")
            sys.exit()
        valid_userccd_files.append(userccd_file)

    return valid_userccd_files


def get_valid_ligand_input_files_with_type(dir_or_file: str, is_batch: bool) -> Valid_Ligand_Input:
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
    file_ext: Literal["json", "a3m"] = "json"
    if is_batch:
        if not os.path.isdir(dir_or_file):
            print(f"\n-- Error: For --batch option, the input query should be a directory with valid JSON/A3M files. "
                  f"The input query is not a directory. Please try again with a valid path.\n")
            sys.exit()
        else:
            files_with_json_ext = get_files_from_path_by_extension(dir_or_file, ".json")
            files_with_a3m_ext = get_files_from_path_by_extension(dir_or_file, ".a3m")
            if len(files_with_json_ext) >= 1 and len(files_with_a3m_ext) >= 1:
                print(f"\n-- Error: Provided path contains both .json and .a3m files. Please use either JSON or A3M.")
                sys.exit()
            if len(files_with_json_ext) >= 1:
                file_paths.extend(files_with_json_ext)
            elif len(files_with_a3m_ext) >= 1:
                file_ext = "a3m"
                file_paths.extend(files_with_a3m_ext)
            else:
                print(f"\n-- Error: No files detected with .json or .a3m extension.")
                sys.exit()
    else:
        if os.path.isfile(dir_or_file):
            if dir_or_file.endswith(".json"):
                file_paths.append(dir_or_file)
            elif dir_or_file.endswith(".a3m"):
                file_ext = "a3m"
                file_paths.append(dir_or_file)
            else:
                print(f"\n-- Error: No files detected with .json or .a3m extension.")
                sys.exit()
        else:
            print(f"\n-- Error: The input query is not a file. "
                  f"Please provide a valid JSON/A3M file path. Or, use '--batch' if provided path is a directory of "
                  f"valid files. \n")
            sys.exit()

    if is_batch:
        print(f"\n-- {current_time()} > Validating {len(file_paths)} input file(s) ...")

    # checking if file is valid
    invalid_files: List[str] = []
    if file_ext == "json":
        for json_file in file_paths:
            if not len(load_json_file(json_file)) > 0:
                invalid_files.append(json_file)
    elif file_ext == "a3m":
        for a3m_file in file_paths:
            if not is_valid_protein_a3m(a3m_file):
                invalid_files.append(a3m_file)

    valid_file_paths = [file_path for file_path in file_paths if file_path not in invalid_files]
    format_invalids = join_list_elements_by_character(invalid_files, "\n")
    if is_batch:
        if len(invalid_files) == len(file_paths):
            print(f"\n-- Error: No valid {file_ext.upper()} files detected in the input directory. "
                  f"Please provide valid input files and try again.\n")
            sys.exit()
        elif len(invalid_files) > 0:
            print(f"\n-- Warning: Invalid {file_ext.upper()} files: \n{format_invalids}"
                  f"\n\n-- Warning: Continuing with {len(valid_file_paths)} valid file(s).")
            valid_files = Valid_Ligand_Input(file_ext, valid_file_paths)
            return valid_files
    else:
        if len(invalid_files) != 0:
            format_invalids = join_list_elements_by_character(invalid_files, "\n")
            print(f"\n-- Error: The following file(s) are not valid {file_ext.upper()} files: \n{format_invalids}\n"
                  f"\n\n-- Error: Please provide valid input file(s) and try again.\n")
            sys.exit()

    valid_files = Valid_Ligand_Input(file_ext, valid_file_paths)
    return valid_files


def make_json_with_ligand(args) -> None:
    input_files = get_valid_ligand_input_files_with_type(args.query, args.batch)
    is_a3m = input_files.file_ext == "a3m"
    is_batch = args.batch

    ligand_to_be_added: Optional[List[List[str]]] = args.add_ligand
    is_purge: bool = args.purge_ligands
    ccdcodes_to_be_removed: Optional[List[str]] = args.remove_ccdcodes
    path_to_userccd_file: Optional[List[str]] = check_userccd_path(args.add_userccd)
    name: Optional[str] = args.name

    # Create output directory
    out_dir_path = manage_output_path(args.output, args.overwrite_existing_results)

    if is_a3m:
        temp_dir = os.path.join(out_dir_path, "temp")
        os.makedirs(temp_dir)
        # process a3m files
        for a3m_file_path in input_files.file_paths:
            a3m_file_basename = os.path.basename(a3m_file_path)
            json_file_name = f"{os.path.splitext(a3m_file_basename)[0]}.json"
            temp_json_file_path = os.path.join(temp_dir, json_file_name)
            ligand_json_file_path = os.path.join(out_dir_path, json_file_name)
            run_msa_to_json(a3m_file_path, temp_json_file_path)
            modjson(temp_json_file_path, ligand_json_file_path, ligand_to_be_added, is_purge,
                    ccdcodes_to_be_removed, name, path_to_userccd_file, is_batch)

        shutil.rmtree(temp_dir)

    else:
        # process json files
        for json_file_path in input_files.file_paths:
            ligand_json_file_path = os.path.join(out_dir_path, os.path.basename(json_file_path))
            modjson(json_file_path, ligand_json_file_path, ligand_to_be_added, is_purge,
                    ccdcodes_to_be_removed, name, path_to_userccd_file, is_batch)

    print(f"\n-- {current_time()} > Process completed. Check output at: '{out_dir_path}'.\n")

    return None

