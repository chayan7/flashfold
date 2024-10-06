# Author: Chayan Kumar Saha

import os
from flashfold.utils import is_valid_path, wget_file_from_url


# Need to update the rep_database_link link later.

rep_database_link: str = ""


def download_representative_db(output_path: str) -> None:
    if not (is_valid_path(output_path)):
        os.makedirs(output_path)
    else:
        wget_file_from_url(rep_database_link, output_path)
