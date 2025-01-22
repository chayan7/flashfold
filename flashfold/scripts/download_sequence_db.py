# Author: Chayan Kumar Saha

import os
import shutil
from flashfold.utils import is_valid_path, wget_file_from_url, load_json_file


def download_database_from_cloud(args) -> None:
    output_path = os.path.realpath(args.output)
    json_file_path = args.input
    if not (is_valid_path(output_path)):
        os.makedirs(output_path)
    else:
        if any(os.scandir(output_path)):
            print(f"Provided output directory '{output_path}' is not empty, please use a different output directory.")
            return

    database_dict = load_json_file(json_file_path)
    print(f"-- Total {len(database_dict)} database(s) will be downloaded. "
          f"Please wait until the process is finished.\n")

    for db_name in database_dict:
        print(f"-- Downloading {db_name} ...")
        database_link = database_dict[db_name]
        wget_file_from_url(database_link, output_path)
        zip_file = f"{output_path}/{db_name}.zip"
        if os.path.exists(f"{output_path}/{db_name}"):
            pass
        else:
            print(f"-- Extracting {db_name} ...")
            shutil.unpack_archive(zip_file, f"{output_path}", "zip")
            os.remove(zip_file)
            print(f"-- Extraction complete for {db_name}")
