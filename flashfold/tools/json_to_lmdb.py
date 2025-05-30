import argparse
from flashfold.utils import convert_dict_to_lmdb, load_json_file


def main():
    parser = argparse.ArgumentParser(description="Convert a JSON file to an LMDB database.")
    parser.add_argument("-j", "--json_file", help="Path to the input JSON file")
    parser.add_argument("-l", "--lmdb_path", help="Path to the output LMDB directory")
    args = parser.parse_args()

    # Load JSON as dict[str, set[str]]

    data = load_json_file(args.json_file)
    data_dict = {k: set(map(str, v)) for k, v in data.items()}
    convert_dict_to_lmdb(data_dict, args.lmdb_path)


if __name__ == "__main__":
    main()
