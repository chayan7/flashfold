# To know more about LMDB, visit: https://lmdb.readthedocs.io/en/release/

import lmdb
import ujson
from .util import remove_all_contents_in_directory, is_valid_path
from typing import Dict, Set


def lmdb_to_dict(lmdb_path: str) -> Dict[str, Set[str]]:
    """
    Convert an LMDB database to a dictionary.
    Args:
        lmdb_path: Path to the LMDB database.

    Returns:
        Dict[str, set[str]]
    """
    env = lmdb.open(lmdb_path, readonly=True, lock=False)
    data = {}
    with env.begin() as txn:
        with txn.cursor() as cursor:
            for key, value in cursor:
                key_decode: str = key.decode()
                item_list = ujson.loads(value)
                data[key_decode] = set(str(item) for item in item_list)
    env.close()
    return data


def estimate_lmdb_map_size(data: Dict[str, Set[str]], overhead_factor: float = 3) -> int:
    """
    Estimate the size of an LMDB map based on the data.
    Args:
        data: A dictionary where keys are strings and values are set of string.
        overhead_factor: Set the overhead factor for LMDB map size estimation. Default is 3.
        User can adjust this based on their needs.

    Returns: Map size in bytes as an integer.
    """
    total_bytes = 0
    for k, v in data.items():
        key_bytes = len(str(k).encode())
        value_bytes = len(ujson.dumps(list(v)).encode())
        total_bytes += key_bytes + value_bytes
    return int(total_bytes * overhead_factor)


def convert_dict_to_lmdb(data_dict: Dict[str, Set[str]], lmdb_file_path: str) -> None:
    """
    Convert a dictionary to an LMDB database.
    Args:
        data_dict: A dictionary where keys are strings and values are set of string.
        lmdb_file_path: Path where the LMDB database will be created.

    Returns: None

    """
    if is_valid_path(lmdb_file_path):
        remove_all_contents_in_directory(lmdb_file_path)

    estimated_map_size = estimate_lmdb_map_size(data_dict)
    print(f"-- Estimated LMDB map size: {estimated_map_size} bytes")

    # Create LMDB environment
    env = lmdb.open(lmdb_file_path, map_size=estimated_map_size)  # Adjust map_size as needed

    # Write data to LMDB
    with env.begin(write=True) as txn:
        for key, value in data_dict.items():
            txn.put(key.encode(), ujson.dumps(list(value)).encode())

    env.close()
    print(f"-- Created LMDB at: {lmdb_file_path}")


def extract_values_from_lmdb(lmdb_path: str, key: str) -> Set[str]:
    """
    Extract values from an LMDB database for a given key.
    Args:
        lmdb_path: Path to the LMDB database.
        key: The key to extract values for.

    Returns:
        Set[str]: A set of values associated with the key.
    """
    env = lmdb.open(lmdb_path, readonly=True, lock=False)
    with env.begin() as txn:
        value = txn.get(key.encode())
        if value is not None:
            return set(ujson.loads(value))
        else:
            print(f"-- Warning: No values found for the {key} in the LMDB database.")
            return set()
