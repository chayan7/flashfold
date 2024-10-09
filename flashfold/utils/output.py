import os
import sys
from .util import is_valid_path, remove_all_contents_in_directory


def manage_output_path(path_to_output: str, overwrite: bool) -> str:
    """
    Manage the output path for the flashfold results.
    Args:
        path_to_output: The output path for the flashfold results.
        overwrite: A boolean value indicating if the existing results should be overwritten or not.

    Returns:
        str: The absolute path to the output directory.

    """
    out_path = os.path.realpath(path_to_output)
    if not is_valid_path(out_path):
        os.makedirs(out_path)
        return out_path
    else:
        if not os.path.isdir(out_path):
            raise ValueError(f"Use a valid path for flashfold output. Current input: '{out_path}' is not a directory.")
        else:
            # List the directory contents; True if the directory is empty, False otherwise.
            if not any(os.scandir(out_path)):
                return out_path
            else:
                if not overwrite:
                    print(f"The provided output directory: {out_path} is not empty. "
                          f"Use '--overwrite_existing_results true' option to run again.")
                    sys.exit()
                else:
                    # Maybe later remove selective contents; for example flashfold_msa can be kept
                    remove_all_contents_in_directory(out_path)
                    return out_path
