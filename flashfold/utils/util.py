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
    """
    Calculate the MD5 hash of a file or a string.
    Args:   
        input_type (str): The type of input to calculate the hash for.
        content (str): The content to calculate the hash for.
    Returns:
        str: The MD5 hash of the input content.
    Raises:
        ValueError: If the input type is not 'path' or 'prot'.
    """
    if input_type == "path":
        # Read the file contents and calculate hash
        with open(content, 'rb') as file:
            md5_hash = hashlib.md5()
            while chunk := file.read(4096):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()
    if input_type == "prot":
        # Calculate hash directly from the provided content string
        string_bytes = content.encode('utf-8')
        md5_hash = hashlib.md5()
        md5_hash.update(string_bytes)
        return md5_hash.hexdigest()
    raise ValueError("Invalid input type. Please provide 'path' or 'prot'.")


def is_valid_path(directory_path: str) -> bool:
    """
    Args:
        directory_path (str): The path to the directory to check.
    Returns:
        bool: True if the directory exists, False otherwise
    """
    if os.path.exists(directory_path):
        return True
    else:
        return False


def create_new_directory(directory_path: str) -> None:
    """
    Args:
        directory_path (str): The path to the directory to create.
    """
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
    path_to_directory = os.path.realpath(directory)
    if not os.path.isdir(path_to_directory):
        raise ValueError(f"The path '{path_to_directory}' is not a directory or does not exist.")

    for filename in os.listdir(path_to_directory):
        file_path = os.path.join(path_to_directory, filename)
        try:
            shutil.rmtree(file_path)  # Removes files, directories, and symlinks recursively
        except OSError:
            os.remove(file_path)


def has_desired_file_extensions(file: str, desired_extensions: List) -> bool:
    """
    Checks if the file has one of the specified extensions.

    Parameters:
    file (str): The file to check.
    desired_extensions (List): The list of desired file extensions.

    Returns:
    bool: True if the file has one of the desired extensions, False otherwise.
    """
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
    extensions (List): The list of desired file extensions.

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


def get_hash_to_files_with_extensions_from_dir(directory: str, extensions: List) -> Dict[str, str]:
    """
    Indexes the files in the specified directory by their MD5 hash.

    Parameters:
    directory (str): The directory to index.
    extensions (List): The list of desired file extensions.

    Returns:
    Dict: A dictionary where keys are MD5 hashes and values are file paths.
    """
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
    """
    Get the current time in a human-readable format.
    """
    present_time = datetime.now()
    return present_time.strftime("%d-%m-%y %H:%M:%S")


def current_time_raw() -> datetime:
    """
    Get the current time in a raw format.
    """
    return datetime.now()


def is_zero_or_pos_int(input_str: str) -> int:
    """
    Check if the input is a 0 or a positive integer.
    """
    try:
        n = int(input_str)
        if n < 0:
            raise argparse.ArgumentTypeError("Input must be a 0 or a positive integer (input >= 0)")
        return n
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid input")


def is_pos_int(input_str: str) -> int:
    """
    Check if the input is a positive integer.
    """
    try:
        n = int(input_str)
        if n < 1:
            raise argparse.ArgumentTypeError("Input must be a positive integer (input >= 1)")
        return n
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid input")


def join_list_elements_by_character(input_list: List,
                                    delimiter: Literal["_", "-", ",", ":", "\t", "\n", "", " ", "null"]) -> str:
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
        raise FileNotFoundError(f"{tool_name} is not installed. Please install it.")


def get_files_from_path_by_extension(path: str, extension: str) -> list:
    """
    Get a list of files with a specified extension from a given path.
    """
    real_path = os.path.realpath(path)
    all_contents_in_path = glob.glob(os.path.join(real_path, "*"))
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
    """
    Get the file path from a list of file paths based on the file name.
    Args:
        file_path_list
        file_name

    Returns:
        str: The file path corresponding to the file name.
    """
    file_path_to_return = ""
    for file_path in file_path_list:
        base_name = os.path.basename(file_path)
        if base_name == file_name:
            file_path_to_return = file_path
        else:
            pass
    return file_path_to_return


def is_pattern_matched(pattern: str, query: str) -> bool:
    """
    Check if a pattern is matched in a query string.
    Args:
        pattern: regex pattern to match
        query: string to search for the pattern

    Returns:
        bool: True if the pattern is matched in the query, False otherwise.
    """
    return bool(re.findall(pattern, query))


def get_directory_from_file_path(file_path: str) -> str:
    """
    Get the directory path from a given file path.
    Args:
        file_path:  The path to the file.

    Returns:
        str: The directory path.
    """
    return os.path.dirname(file_path)


def replace_char_from_string(input_string: str, replace_with: Literal["_", "-", ":"]) -> str:
    """
    Replace special characters in a string with a specified character.
    Args:
        input_string:       The input string.
        replace_with:       The character to replace the special characters with.

    Returns:
        str: The string with special characters replaced.
    """
    return re.sub("[^a-zA-Z0-9._]", replace_with, input_string)


def update_time_log(log_file: str, message: str, use_timestamp: bool) -> None:
    """
    Update the time log with a message.
    Args:
        log_file:    The path to the log file.
        message:    The message to write to the log file.
        use_timestamp:  A boolean value to determine if the message should be timestamped.

    Returns:
        None
    """
    with open(log_file, "a") as f:
        if use_timestamp:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"[{timestamp}] {message}\n")
        else:
            f.write(f"{message}\n")


def get_contents_by_column(file_path: str, column: int = 1) -> List[str]:
    """
    Get a list of contents from a file.
    Args:
        file_path: The path to the file.
        column: The column number to extract the contents from. Columns are "not 0-indexed".

    Returns:
        List: A list of contents from the file.
    """
    if column < 1:
        raise ValueError("Column number should be greater than 0.")
    contents = []
    with open(file_path, "r") as in_file:
        for line in in_file:
            processed_line = line.strip().split('\t')[column - 1]
            contents.append(processed_line)
    return contents


def get_sum_of_all_file_sizes_by_path(path: str) -> int:
    """
    Check the size of all files in a given path.

    :param path: Path to the directory
    Returns: total size of all files in the directory
    """
    total_size = 0
    path_to_directory = os.path.realpath(path)
    if not os.path.isdir(path_to_directory):
        print(f"The path '{path_to_directory}' is not a directory or does not exist.")
        return total_size

    for root, dirs, files in os.walk(path):
        for file in files:
            file_path = os.path.join(root, file)
            file_size = os.path.getsize(file_path)
            total_size += file_size
    return total_size


def file_has_content(file_path: str) -> bool:
    """
    Check if a file has content.
    Args:
        file_path: The path to the file.

    Returns: True if it has content
    """
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0
