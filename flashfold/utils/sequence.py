# noinspection PyPackageRequirements
from Bio import SeqIO
# noinspection PyPackageRequirements
from Bio.Seq import Seq
# noinspection PyPackageRequirements
from Bio.SeqRecord import SeqRecord
from typing import List, Dict
from collections import namedtuple
from .util import calculate_md5_hash, join_list_elements_by_character, replace_char_from_string, is_valid_path

# Define the named tuple
Sequence = namedtuple('Sequence', ['hash', 'fasta'])
Infile_feats = namedtuple('Infile_feats', ['accnrs', 'seqs', 'chain_accnrs', 'chain_seqs',
                                           'chain_seq_hashes', 'a3m_header', 'hash_to_fasta', 'empty_subunits'])
Chain = namedtuple('Subunit', ['subunits', 'frequency', 'sequence'])
Fasta_record = namedtuple('Fasta_record', ['accession', 'fasta'])


def extract_protein_sequences(gbff_file_path: str) -> List:
    """
    Extract protein sequences from a GenBank file. The protein sequences are extracted from the CDS features in the
    GenBank file.

    Args:
        gbff_file_path (str): Path to the GenBank file.

    Returns:
        List: A list of tuples containing accession, gene name, product, protein sequence, and hash value.

    """
    protein_sequences = []
    with open(gbff_file_path, "r") as file:
        for record in SeqIO.parse(file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    accession = feature.qualifiers["protein_id"][0]
                    gene_name = feature.qualifiers.get("gene", [""])[0]
                    product = feature.qualifiers.get("product", [""])[0]
                    protein_seq = feature.qualifiers.get("translation", [""])[0]
                    hash_value = calculate_md5_hash("prot", protein_seq)
                    protein_sequences.append((accession, gene_name, product, protein_seq, hash_value))
    return protein_sequences


def create_fasta_for_db(accession: str, gene: str, desc: str, hash_str: str, sequence: str) -> str:
    """
    Converts a sequence to a FASTA format.
    Args:
        accession
        gene
        desc
        hash_str
        sequence

    Returns:
        str: The sequence in FASTA format

    """
    seq = Seq(sequence)
    desc_with_hash = f"{hash_str}: {desc}"
    seq_record = SeqRecord(seq, id=accession, name=gene, description=desc_with_hash)
    return seq_record.format("fasta")


def is_protein_sequence(sequence: str, is_seq_from_msa: bool = False) -> bool:
    """
    Checks if a given sequence is a valid protein sequence.

    Parameters:
    sequence (str): The sequence to be checked.

    Returns:
    bool: True if the sequence is a valid protein sequence, False otherwise.
    """

    # According to https://wiki.thegpm.org/wiki/Amino_acid_symbols, two additional amino acid code is added to the
    # end of the valid list: U=Selenocysteine; O=Pyrrolysine

    valid_amino_acids = set("ARNDCEQGHILKMFPSTWYVUO")

    if is_seq_from_msa:
        extended_aa = set("BJXZ")
        valid_amino_acids.update(extended_aa)

    sequence = sequence.upper()  # Convert sequence to uppercase

    # Remove the DNA/RNA characters from the sequence
    dna_rna_chars = set("ATGCU")
    sequence_without_dna_rna_char = "".join([char for char in sequence if char not in dna_rna_chars])

    if sequence_without_dna_rna_char == "":
        print("Invalid sequence: input sequence is either DNA or RNA.")
        return False

    for char in sequence_without_dna_rna_char:
        if char not in valid_amino_acids:
            print(f"Invalid character found in sequence: {char}")

    return all(char in valid_amino_acids for char in sequence)


def is_valid_protein_fasta(file_path: str) -> bool:
    """
    Validates a protein FASTA file.

    Parameters:
    - file_path: Path to the FASTA file to be validated.

    Returns:
    - bool: True if the file is a valid protein FASTA file, False otherwise.
    """
    if not is_valid_path(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    try:
        with open(file_path, 'r') as fasta_file:
            lines = fasta_file.readlines()

        if not lines or not lines[0].startswith('>'):
            raise ValueError(f"File does not start with a '>' character.\nCheck: {file_path}")

        id_count = 0
        for i, line in enumerate(lines):
            line = line.strip()
            if line.startswith('>'):
                if len(line) > 1:
                    id_count += 1
            elif not is_protein_sequence(line):
                raise ValueError(
                    f"Invalid sequence character(s) found in line {i + 1} of the input protein FASTA file."
                    f"\nCheck: {file_path}")

        if id_count == 0:
            raise ValueError(f"No valid protein IDs found in the file.\nCheck: {file_path}")

        return True

    except Exception as e:
        print(f"\n-- An error occurred: {e}")
        return False


def is_valid_protein_a3m(file_path: str) -> bool:
    """
    Validates a protein A3M file.

    Parameters:
    - file_path: Path to the A3M file to be validated.

    Returns:
    - bool: True if the file is a valid protein A3M file, False otherwise.
    """
    null_character = chr(0)
    if not is_valid_path(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    try:
        with open(file_path, 'r') as a3m_file:
            lines = a3m_file.readlines()

        if not lines or not lines[0].startswith('#'):
            raise ValueError(f"File does not start with a '#' character.\nCheck: {file_path}")

        a3m_id_count = 0
        for i, raw_line in enumerate(lines[1:]):
            line = raw_line.strip()
            if line.startswith('>'):
                if len(line) > 1:
                    a3m_id_count += 1
            elif not is_protein_sequence(line.replace("-", "").replace(null_character, "").upper(), True):
                raise ValueError(
                    f"Invalid sequence character(s) found in line {i + 1} of the input protein a3m file."
                    f"\nCheck: {file_path}")

        if a3m_id_count == 0:
            raise ValueError(f"No valid protein IDs found in the file.\nCheck: {file_path}")

        return True

    except Exception as e:
        print(f"\n-- An error occurred: {e}")
        return False


def is_a3m_monomer(a3m_file_path: str) -> bool:
    with open(a3m_file_path, 'r') as a3m_file:
        for line in a3m_file:
            if line.startswith("#"):
                split_line = line.rstrip().split("\t")
                if split_line[1] == "1":
                    return True
    return False


def get_valid_sequence_records_from_fasta(fasta_file: str) -> List[Dict]:
    """
    Get valid sequence records from a FASTA file. The sequence records are stored in a list of dictionaries.
    Args:
        fasta_file:     Path to the FASTA file.

    Returns:
        List: A list of dictionaries containing the accession, sequence, FASTA format, and hash value of the sequence.
    """
    parsed_sequence = SeqIO.parse(fasta_file, "fasta")
    seq_records = []
    count = 0
    for record in parsed_sequence:
        accession = str(record.id)
        sequence = str(record.seq)
        count += 1
        accession_without_bad_char = replace_char_from_string(accession, "_")
        changed_accession = f"S{count}_{accession_without_bad_char}"
        record.id = changed_accession
        seq_hash = calculate_md5_hash("prot", str(sequence))
        record_dict: Dict[str, str] = {"accession": changed_accession,
                                       "sequence": sequence,
                                       "fasta": record.format("fasta"),
                                       "seq_hash": seq_hash}
        seq_records.append(record_dict)
    return seq_records


def get_chain_from_sequence(sequence_list: List[str]) -> Chain:
    """
    Get a chain of subunits from a list of sequences.
    Args:
        sequence_list:  A list of sequences.

    Returns:
        Chain: A chain of subunits.
    """
    uniq_sequence_set = set()
    abundances = []
    seq_to_abundance = []
    for seq in sequence_list:
        if seq not in uniq_sequence_set:
            abundance = sequence_list.count(seq)
            abundances.append(abundance)
            seq_to_abundance.append(
                [seq, abundance]
            )
            uniq_sequence_set.add(seq)

    min_abundance = min(abundances)

    subunits = []
    for seq, abundance in seq_to_abundance:
        rational_abundance = round(abundance/min_abundance)
        for i in range(rational_abundance):
            subunits.append(seq)

    subunit_sequence = join_list_elements_by_character(subunits, "")
    total_sequence = join_list_elements_by_character(sequence_list, "")
    frequency = total_sequence.count(subunit_sequence)
    if frequency == 0:
        return Chain(sequence_list, 1, total_sequence)
    return Chain(subunits, frequency, subunit_sequence)


def get_input_fasta_features(sequence_records: List[Dict]) -> Infile_feats:
    """
    Get input FASTA features.
    Args:
        sequence_records:  A list of dictionaries containing the accession, sequence,
        FASTA format, and hash value of the sequence.

    Returns:
        Infile_feats: A named tuple containing the input FASTA features.
    """
    accessions = []
    sequences = []
    seq_hashes = []
    hash_to_fasta = {}
    for record in sequence_records:
        accessions.append(record["accession"])
        sequences.append(record["sequence"])
        seq_hashes.append(record["seq_hash"])
        hash_to_fasta[record["seq_hash"]] = record["fasta"]

    # Make single chain with subúnits from input sequences
    # It generates subunits as a list, frequency as integer, subunit_sequence as string
    chain = get_chain_from_sequence(sequences)

    chain_accnrs = []
    chain_seq_hashes = []
    chain_sequences = chain.subunits
    for seq in chain_sequences:
        chain_accnrs.append(accessions[sequences.index(seq)])
        chain_seq_hashes.append(seq_hashes[sequences.index(seq)])

    # Make unpaired alignment
    empty_subunits = [""] * len(chain_seq_hashes)

    for i in range(len(chain_sequences)):
        empty_seq = "-" * len(chain_sequences[i])
        empty_subunits[i] = empty_seq

    if chain.frequency == 1:
        uniq_seq_lengths = []
        uniq_seq_units = []
        for uniq_seq in chain_sequences:
            uniq_seq_lengths.append(len(uniq_seq))
            uniq_seq_units.append(chain.frequency)
        joined_uniq_seq_lengths = join_list_elements_by_character(uniq_seq_lengths, ",")
        joined_uniq_seq_units = join_list_elements_by_character(uniq_seq_units, ",")
        a3m_header = f"#{joined_uniq_seq_lengths}\t{joined_uniq_seq_units}"
    else:
        a3m_header = f"#{len(chain.sequence)}\t{chain.frequency}"

    return Infile_feats(accessions, sequences, chain_accnrs, chain_sequences, chain_seq_hashes, a3m_header,
                        hash_to_fasta, empty_subunits)


def get_records(file_path: str) -> Dict[str, str]:
    """
    Get alignment records from a file. The records are stored in a dictionary.
    Args:
        file_path: Path to the FASTA/A3M file.

    Returns:
        Dict: A dictionary containing the accession and sequence
    """
    parsed_sequence = SeqIO.parse(file_path, "fasta")
    records = {}
    for record in parsed_sequence:
        accession = str(record.id)
        sequence = str(record.seq)
        # Accession not empty validity check
        if accession == "":
            pass
        elif sequence == "":
            pass
        else:
            records[accession] = sequence
    return records


def get_sequence_length_from_single_fasta(fasta_file: str) -> int:
    """
    Get the length of a sequence from a single FASTA file.
    Args:
        fasta_file: Path to the FASTA file.

    Returns:
        int: The length of the sequence.
    """
    sequence_record = get_records(fasta_file)
    if len(sequence_record) > 1:
        raise ValueError("More than one sequence found in the file.")
    sequence_length = len(list(sequence_record.values())[0])
    return sequence_length


def combine_sequences(accessions: List, sequences: List) -> Sequence:
    """
    Combine sequences from a list of accessions and sequences
    Args:
        accessions:     A list of accessions.
        sequences:    A list of sequences.

    Returns:
        Sequence: A named tuple containing the hash value and the combined sequences.
    """
    accession_combo = join_list_elements_by_character(accessions, "\t")
    joined_combo_sequence = join_list_elements_by_character(sequences, "null")
    hash_combo_sequence = calculate_md5_hash("prot", joined_combo_sequence)
    sequences_formatted = f">{accession_combo}\n{joined_combo_sequence}\n"
    return Sequence(hash_combo_sequence, sequences_formatted)


def make_hash_fasta_sequence(accession: str, sequence: str) -> Sequence:
    """
    Make a FASTA sequence from an accession and a sequence.
    Args:
        accession:
        sequence:

    Returns:
        Sequence: A named tuple containing the hash value and the FASTA sequence.
    """
    hash_combo_sequence = calculate_md5_hash("prot", sequence)
    sequences_formatted = f">{accession}\n{sequence}\n"
    return Sequence(hash_combo_sequence, sequences_formatted)


def combine_gappy_sequences(accessions: List, sequences: List) -> Sequence:
    """
    Combine gappy sequences from a list of accessions and sequences.
    Args:
        accessions:     A list of accessions.
        sequences:  A list of sequences.

    Returns:
        Sequence: A named tuple containing the hash value and the combined gappy sequences
    """
    accession_index = 0
    for i in range(len(sequences)):
        seq = sequences[i]
        mod_seq = seq.replace("-", "")
        if mod_seq != "":
            accession_index = i

    accession = accessions[accession_index]
    joined_combo_sequence = join_list_elements_by_character(sequences, "null")
    hash_combo_sequence = calculate_md5_hash("prot", joined_combo_sequence)
    sequences_formatted = f">{accession}\n{joined_combo_sequence}\n"
    return Sequence(hash_combo_sequence, sequences_formatted)


class SequenceDbFasta:
    def __init__(self, fasta_path: str) -> None:
        protein_hash_to_record = {}
        for record in SeqIO.parse(fasta_path, "fasta"):
            protein_hash = calculate_md5_hash("prot", str(record.seq))
            fasta_record = Fasta_record(record.id, record.format("fasta"))
            protein_hash_to_record[protein_hash] = fasta_record
        self.sequence_records = protein_hash_to_record

    def get_record_by_protein_hash(self, protein_hash: str) -> Fasta_record:
        if protein_hash in self.sequence_records:
            return self.sequence_records[protein_hash]
        raise ValueError("Warning: File does not contain sequence.")

