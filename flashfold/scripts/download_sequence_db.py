# Author: Chayan Kumar Saha

import os
from flashfold.utils import is_valid_path, wget_file_from_url, load_json_file


def download_database_from_cloud(args) -> None:
    output_path = args.output
    json_file_path = args.input
    if not (is_valid_path(output_path)):
        os.makedirs(output_path)
    database_dict = load_json_file(json_file_path)
    print(f"-- Total {len(database_dict)} database(s) will be downloaded. "
          f"Please wait until the process is finished.\n")
    for db_name in database_dict:
        database_link = database_dict[db_name]
        wget_file_from_url(database_link, os.path.abspath(output_path))
