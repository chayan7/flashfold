from collections import defaultdict, namedtuple
from .util import (is_valid_path, get_filename_to_path_set_by_directory, get_files_from_path_by_extension,
                   get_filename_without_extension, is_pattern_matched, get_file_dir_by_file_path)
from .json import write_dict_to_json_as_file
from .lmdb import extract_values_from_lmdb
from typing import Dict, List, Set

Db_Content = namedtuple('Db_Content', ['protein_hash', 'is_new_protein',
                                       'new_accessions', 'new_gbks', 'new_fasta'])

files_to_be_in_database = ["prot_hash_to_accession.json", "sequence_db.fasta", "data.mdb", "lock.mdb"]


def is_valid_database_file_count(db_file_list: List[str], query_dict: Dict[str, Set[str]]) -> bool:
    """
    Checks if the database file count is valid.

    Args:
        db_file_list (List[str]): List of database file names.
        query_dict (Dict[str, Set[str]]): Dictionary mapping file names to sets of file paths.

    Returns:
        bool: True if each file name has exactly one file path, False otherwise.
    """
    for filename in db_file_list:
        count = len(query_dict[filename])
        # Database should have 1 file path per file name
        if count != 1:
            return False
    return True


def is_valid_database_dir(database_dir: str) -> bool:
    """
    Checks if the database directory is valid.

    Args:
        database_dir (str): Path to the database directory.

    Returns:
        bool: True if the database directory is valid, False otherwise.
    """
    if is_valid_path(database_dir):
        filename_to_path = get_filename_to_path_set_by_directory(database_dir, [".json", ".fasta", ".mdb"])
        if not is_valid_database_file_count(files_to_be_in_database, filename_to_path):
            print(f"Invalid sequence database detected, check: {database_dir} "
                  f"\n- Download latest version of database using flashfold featured download_db subcommand, or "
                  f"\n- Create database using the create_db subcommand.")
            return False
        else:
            return True
    else:
        print(f"Invalid sequence database detected, check: {database_dir} "
              f"\n- Download latest version of database using flashfold featured download_db subcommand, or "
              f"\n- Create database using the create_db subcommand.")
        return False


class Database:
    def __init__(self, path: str) -> None:
        if not is_valid_database_dir(path):
            raise ValueError(f"Invalid database directory: {path}")
        self.database_path = path
        self.database_files = get_filename_to_path_set_by_directory(self.database_path, [".fasta", ".json", ".mdb"])
        self.fasta_db = self._sequence_db()
        self.protein_to_gbks = self._protein_to_gbks()

    def _sequence_db(self) -> str:
        return list(self.database_files["sequence_db.fasta"])[0]

    def _protein_to_gbks(self) -> str:
        return get_file_dir_by_file_path(list(self.database_files["data.mdb"])[0])

    def process_homology_search_output(self, path_to_alignment: str, query_seq_hashes: List[str],
                                       json_out_file: str, threads: int) -> None:
        a3m_files = get_files_from_path_by_extension(path_to_alignment, ".a3m")
        if len(a3m_files) == 0:
            raise FileNotFoundError(f"No .a3m files found in the directory: {path_to_alignment}")
        gbk_to_hits: Dict[str, List[str]] = defaultdict(list)

        query_seq_hash_to_a3m_file: Dict[str, str] = \
            {get_filename_without_extension(a3m_file): a3m_file for a3m_file in a3m_files}

        hit_hash_keys = []
        query_hash_colon_hit_accessions = []
        for query_seq_hash in query_seq_hashes:
            a3m_file_path = query_seq_hash_to_a3m_file.get(query_seq_hash)
            if not a3m_file_path:
                continue
            with (open(a3m_file_path, "r", encoding="utf-8") as a3m_in):
                for line in a3m_in:
                    if not line.startswith(">"):
                        continue
                    split_line = line.strip().split("\t")
                    if not len(split_line) > 1:
                        continue
                    hit_accession = split_line[0][1:]
                    hit_hash_key = split_line[1]
                    slash_digit_to_digit_pattern = r'/\d+-\d+'
                    if not is_pattern_matched(slash_digit_to_digit_pattern, hit_accession):
                        continue
                    # print(hit_accession, hit_hash_key)
                    query_hash_colon_hit_accession = f"{query_seq_hash}:{hit_accession}"
                    hit_hash_keys.append(hit_hash_key)
                    query_hash_colon_hit_accessions.append(query_hash_colon_hit_accession)

        hit_hash_keys_to_gbk = extract_values_from_lmdb(self.protein_to_gbks, set(hit_hash_keys), threads)

        for query_hash_colon_hit_accession, hit_hash_key in zip(query_hash_colon_hit_accessions, hit_hash_keys):
            for gbk in hit_hash_keys_to_gbk.get(hit_hash_key, []):
                if query_hash_colon_hit_accession not in gbk_to_hits[gbk]:
                    gbk_to_hits[gbk].append(query_hash_colon_hit_accession)

        write_dict_to_json_as_file(gbk_to_hits, json_out_file)

        return None


class CreateDbContent:
    def __init__(self, protein_hash: str, is_new_protein: bool, new_accessions: List[str],
                 new_gbks: List[str], new_fasta: str) -> None:
        self.protein_hash = protein_hash
        self.is_new_protein = is_new_protein
        self.new_accessions = new_accessions
        self.new_gbks = new_gbks
        self.new_fasta = new_fasta

    def get_formatted_content(self) -> Db_Content:
        return Db_Content(self.protein_hash, self.is_new_protein, self.new_accessions, self.new_gbks, self.new_fasta)
