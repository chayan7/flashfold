# To know more about LMDB, visit: https://lmdb.readthedocs.io/en/release/

import lmdb
from .util import remove_all_contents_in_directory, is_valid_path
from typing import Dict, Set, List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed


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
                value_decode = set(value.decode().split(","))
                data[key_decode] = value_decode
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
        value_bytes = len(",".join(list(v)))
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
            str_value = ",".join(list(value))
            txn.put(key.encode(), str_value.encode())

    env.close()
    print(f"-- Created LMDB at: {lmdb_file_path}")


def extract_values_from_lmdb(lmdb_path: str, keys: Set[str], max_workers: int) -> Dict[str, List[str]]:
    key_to_values = {}
    env = lmdb.open(lmdb_path, readonly=True, lock=False, max_readers=512)
    with env.begin() as txn:

        # Function to fetch values for a given key
        def fetching_lmdb(key: str) -> Tuple[str, List[str]]:
            value = txn.get(key.encode(), )
            if value is not None:
                return key, value.decode().split(",")
            return key, []

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(fetching_lmdb, key): key for key in keys}
            for future in as_completed(futures):
                key, values = future.result()
                if values is not None:
                    key_to_values[key] = values
    env.close()
    return key_to_values
