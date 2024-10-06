import os
import os.path
import json
import ujson
from typing import Dict


# Function to load stored history info from file
def load_json_file(file_path: str) -> Dict:
    """
    Load a JSON file and return its content as a dictionary
    Args:
        file_path:  The path to the JSON file to load.

    Returns:
        Dict: The content of the JSON file as a dictionary.
    """
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            return ujson.load(file)
    else:
        return {}


def write_dict_to_json_as_file(recent_info: Dict, file_path: str) -> None:
    """
    Save a dictionary to a JSON file.
    Args:
        recent_info (Dict): The dictionary to save.
        file_path (str): The path to the file where the dictionary will be saved.
    """
    with open(file_path, 'w') as file:
        # noinspection PyTypeChecker
        json.dump(recent_info, file, indent=4)
    return None


def write_dict_of_set_to_json_as_file(data: Dict, file_path: str) -> None:
    """
    Write a dictionary of sets to a JSON file.
    Args:
        data:   The dictionary of sets to write.
        file_path:  The path to the file where the dictionary will be saved.

    Returns:
        None
    """
    regular_dict = {key: list(value) for key, value in data.items()}
    write_dict_to_json_as_file(regular_dict, file_path)
