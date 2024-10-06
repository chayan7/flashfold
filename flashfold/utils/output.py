import os
import sys
from .util import is_valid_path, remove_all_contents_in_directory


def manage_output_path(out_path: str, overwrite: bool) -> str:
    """
    Manage the output path for the foldflash results.
    Args:
        out_path: The output path for the foldflash results.
        overwrite: A boolean value indicating if the existing results should be overwritten or not.

    Returns:
        str: The absolute path to the output directory.

    """
    if not is_valid_path(out_path):
        os.makedirs(out_path)
        print(f"Created new output directory: {os.path.abspath(out_path)}")
        return os.path.abspath(out_path)
    else:
        if not os.path.isdir(out_path):
            raise ValueError(f"Use a valid path for foldflash output. Current input: '{out_path}' is not a directory.")
        else:
            # List the directory contents; True if the directory is empty, False otherwise.
            if not any(os.scandir(out_path)):
                return os.path.abspath(out_path)
            else:
                if not overwrite:
                    print(f"The provided output directory: {os.path.abspath(out_path)} is not empty. "
                          f"Use '--overwrite_existing_results true' option to run again.")
                    sys.exit()
                else:
                    # Maybe later remove selective contents; for example foldflash_msa can be kept
                    remove_all_contents_in_directory(out_path)
                    return os.path.abspath(out_path)
