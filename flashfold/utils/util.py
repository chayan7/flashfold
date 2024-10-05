import glob
import hashlib
import os
import re
import shutil
import argparse
from datetime import datetime
from collections import defaultdict
from typing import Dict, Set, List, Literal


def calculate_md5_hash(input_type: Literal["path", "prot"], content: str) -> str:
    # Calculate MD5 hash based on input type
    if input_type == "path":
        # Read the file contents and calculate hash
        with open(content, 'rb') as file:
            md5_hash = hashlib.md5()
            while chunk := file.read(4096):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()
    elif input_type == "prot":
        # Calculate hash directly from the provided content string
        string_bytes = content.encode('utf-8')
        md5_hash = hashlib.md5()
        md5_hash.update(string_bytes)
        return md5_hash.hexdigest()
    else:
        raise ValueError("Input type must be either 'path' or 'prot'")


def is_valid_path(directory_path: str) -> bool:
    if os.path.exists(directory_path):
        return True
    else:
        return False


def create_new_directory(directory_path: str) -> None:
    if not is_valid_path(directory_path):
        os.makedirs(directory_path)
    else:
        pass


def remove_all_contents_in_directory(directory: str):
    """
    Removes all contents in the specified directory, including files, directories, and symbolic links.

    Parameters:
    directory (str): The path to the directory to clear.

    Raises:
    ValueError: If the provided path is not a directory.
    """
    if not os.path.isdir(directory):
        raise ValueError(f"The path '{directory}' is not a directory or does not exist.")

    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            shutil.rmtree(file_path)  # Removes files, directories, and symlinks recursively
        except OSError:
            os.remove(file_path)


def has_desired_file_extensions(file: str, desired_extensions: List) -> bool:
    file_extension = os.path.splitext(file)[1]
    if file_extension in desired_extensions:
        return True
    else:
        return False


def get_filename_to_path_set_by_directory(input_directory: str, extensions: List) -> Dict[str, Set[str]]:
    """
    Maps filenames to their absolute paths in the given directory.

    Parameters:
    input_directory (str): The directory to search for files.

    Returns:
    Dict[str, Set[str]]: A dictionary where keys are filenames and values are sets of their absolute paths.
    """
    absolute_path = os.path.abspath(input_directory)
    subdirectories = glob.glob(os.path.join(absolute_path, "*"))
    filename_to_path: Dict[str, Set[str]] = defaultdict(set)
    for potential_file in subdirectories:
        if os.path.isfile(potential_file):
            if has_desired_file_extensions(potential_file, extensions):
                filename_to_path[os.path.basename(potential_file)].add(potential_file)
        elif os.path.isdir(potential_file):
            potential_database_dir = glob.glob(os.path.join(potential_file, "*"))
            for potential_db_file in potential_database_dir:
                if has_desired_file_extensions(potential_db_file, extensions):
                    filename_to_path[os.path.basename(potential_db_file)].add(potential_db_file)
    return filename_to_path


def get_hash_to_files_with_extensions_from_dir(directory: str, extensions: List) -> dict:
    hash_to_desired_files = {}

    if is_valid_path(directory):
        absolute_path = os.path.abspath(directory)
        subdirectories = glob.glob(os.path.join(absolute_path, "*"))
        count = 0
        for potential_file in subdirectories:
            if os.path.isfile(potential_file):
                if has_desired_file_extensions(potential_file, extensions):
                    count += 1
                    print(f"-- indexing genbank file {count}: {potential_file}")
                    hash_to_desired_files[calculate_md5_hash("path", potential_file)] = potential_file
            elif os.path.isdir(potential_file):
                potential_assemblies = glob.glob(os.path.join(potential_file, "*"))
                for file in potential_assemblies:
                    if has_desired_file_extensions(file, extensions):
                        count += 1
                        print(f"-- indexing genbank file {count}: {file}")
                        hash_to_desired_files[calculate_md5_hash("path", file)] = file
    return hash_to_desired_files


def current_time() -> str:
    # Get the current time
    present_time = datetime.now()

    # Format the current time in a human-readable way
    return present_time.strftime("%d-%m-%y %H:%M:%S")


def current_time_raw() -> datetime:
    return datetime.now()


def is_zero_or_pos_int(input_str: str):
    try:
        n = int(input_str)
        if n < 0:
            raise argparse.ArgumentTypeError("Input must be a 0 or a positive integer (input >= 0)")
        else:
            return n
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid input")


def is_pos_int(input_str: str):
    try:
        n = int(input_str)
        if n < 1:
            raise argparse.ArgumentTypeError("Input must be a positive integer (input >= 1)")
        else:
            return n
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid input")


def join_list_elements_by_character(input_list: list,
                                    delimiter: Literal["_", "-", ",", ":", "\t", "\n", "", " ", "null"]):
    """
    Joins elements of the list into a single string with a specified delimiter.

    Parameters:
    lst (list): The list of elements to join.
    delimiter (str): The character to use as the delimiter.

    Returns:
    str: The joined string.
    """
    if delimiter == "null":
        null_character = chr(0)
        return null_character.join(map(str, input_list))

    return delimiter.join(map(str, input_list))


def get_new_elements_from_second_list(first_list: List, second_list: List) -> List:
    set1 = set(first_list)
    set2 = set(second_list)
    # Get the difference
    return list(set2 - set1)


def is_installed(tool_name: str) -> bool:
    """
    Check if a tool is installed by verifying its presence in the system's PATH.

    Parameters:
    tool_name (str): The name of the tool to check.

    Returns:
    bool: True if the tool is installed, False otherwise.
    """
    # Check if the tool is available in the PATH
    tool_path = shutil.which(tool_name)
    if tool_path:
        return True
    else:
        return False


def get_files_from_path_by_extension(path: str, extension: str) -> list:
    absolute_path = os.path.abspath(path)
    all_contents_in_path = glob.glob(os.path.join(absolute_path, "*"))
    file_list = []
    for potential_file in all_contents_in_path:
        if os.path.isfile(potential_file) and has_desired_file_extensions(potential_file, [extension]):
            file_list.append(potential_file)
    return file_list


def get_filename_without_extension(file_path: str) -> str:
    """
    Get the file name without its extension from a given file path.

    Parameters:
    file_path (str): The full path to the file.

    Returns:
    str: The file name without its extension.
    """
    base_name = os.path.basename(file_path)  # Extracts the file name from the path
    file_name_without_ext = os.path.splitext(base_name)[0]  # Removes the extension
    return file_name_without_ext


def get_file_path_by_filename(file_path_list: list, file_name) -> str:
    file_path_to_return = ""
    for file_path in file_path_list:
        base_name = os.path.basename(file_path)
        if base_name == file_name:
            file_path_to_return = file_path
        else:
            pass
    return file_path_to_return


def is_pattern_matched(pattern: str, query: str) -> bool:
    return bool(re.findall(pattern, query))


def get_directory_from_file_path(file_path: str) -> str:
    return os.path.dirname(file_path)


def replace_char_from_string(input_string: str, replace_with: Literal["_", "-", ":"]) -> str:
    return re.sub("[^a-zA-Z0-9._]", replace_with, input_string)


def update_time_log(log_file: str, message: str, use_timestamp: bool) -> None: 
    with open(log_file, "a") as f:
        if use_timestamp:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"[{timestamp}] {message}\n")
        else:
            f.write(f"{message}\n")
