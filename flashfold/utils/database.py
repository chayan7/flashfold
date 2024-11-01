import threading
from collections import defaultdict, namedtuple
from .util import (is_valid_path, get_filename_to_path_set_by_directory, get_files_from_path_by_extension,
                   get_filename_without_extension, is_pattern_matched)
from .json import load_json_file, write_dict_to_json_as_file
from typing import Dict, List, Set


Db_Content = namedtuple('Db_Content', ['protein_hash', 'is_new_protein',
                                       'new_accessions', 'new_gbks', 'new_fasta'])


files_to_be_in_database = ["prot_hash_to_accession.json", "protein_to_gbks.json", "sequence_db.fasta"]


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
        filename_to_path = get_filename_to_path_set_by_directory(database_dir, [".json", ".fasta"])
        if not is_valid_database_file_count(files_to_be_in_database, filename_to_path):
            print(f"Invalid sequence database detected, check: {database_dir} "
                  f"\nTo create database please use the create_db command provided with flashfold.")
            return False
        else:
            return True
    else:
        print(f"Invalid sequence database detected, check: {database_dir} "
              f"\nTo create database please use the create_db command provided with flashfold.")
        return False


class Database:
    def __init__(self, path: str) -> None:
        if not is_valid_database_dir(path):
            raise ValueError(f"Invalid database directory: {path}")
        self.database_path = path
        self.database_files = get_filename_to_path_set_by_directory(self.database_path, [".fasta", ".json"])
        self.fasta_db = self._sequence_db()
        self.prot_hash_to_accession = self._prot_hash_to_accession()
        self.protein_to_gbks = self._protein_to_gbks()
        self._protein_to_gbks_loaded = None
        self._protein_to_gbks_thread = threading.Thread(target=self._load_protein_to_gbks_in_background)
        self._protein_to_gbks_thread.start()

    def _sequence_db(self) -> str:
        seq_db_fasta_path_set = self.database_files["sequence_db.fasta"]
        return list(seq_db_fasta_path_set)[0]

    def _prot_hash_to_accession(self) -> str:
        prot_hash_to_accession_path_set = self.database_files["prot_hash_to_accession.json"]
        return list(prot_hash_to_accession_path_set)[0]

    def _protein_to_gbks(self) -> str:
        protein_to_gbks_json_set = self.database_files["protein_to_gbks.json"]
        return list(protein_to_gbks_json_set)[0]

    def _load_protein_to_gbks_in_background(self) -> None:
        self._protein_to_gbks_loaded = load_json_file(self.protein_to_gbks)

    @property
    def load_protein_to_gbks(self) -> Dict[str, List[str]]:
        # Ensure the background thread has completed
        self._protein_to_gbks_thread.join()
        return self._protein_to_gbks_loaded

    def process_homology_search_output(self, path_to_alignment: str, query_seq_hashes: List[str],
                                       json_out_file: str) -> None:
        sto_files = get_files_from_path_by_extension(path_to_alignment, ".sto")
        if len(sto_files) == 0:
            raise FileNotFoundError(f"No .sto files found in the directory: {path_to_alignment}")
        gbk_to_hits: Dict[str, List[str]] = defaultdict(list)
        for query_seq_hash in query_seq_hashes:
            for sto_file_path in sto_files:
                query_hash = get_filename_without_extension(sto_file_path)
                if query_seq_hash == query_hash:
                    with open(sto_file_path, "r", encoding="utf-8") as sto_in:
                        for line in sto_in:
                            if line.startswith("#=GS"):
                                hit_accession = line.split(":")[0].split(" ")[1]
                                hit_hash_key = line.split(":")[0].split(" ")[-1]
                                slash_digit_to_digit_pattern = r'/\d+-\d+'
                                if is_pattern_matched(slash_digit_to_digit_pattern, hit_accession):
                                    query_hash_colon_hit_accession = f"{query_seq_hash}:{hit_accession}"
                                    for gbk in self.load_protein_to_gbks[hit_hash_key]:
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
