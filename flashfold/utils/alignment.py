import os
import itertools
from typing import List, Dict, Optional
from .util import get_filename_without_extension, join_list_elements_by_character, current_time
from .sequence import Sequence, Infile_feats, get_records, combine_sequences, \
    combine_gappy_sequences, make_hash_fasta_sequence
from .json import load_json_file
from collections import namedtuple


A3M_Records = namedtuple("A3M_records", ["query_hash_to_seq", "hits_per_query"])


def has_good_coverage(sequence: str, coverage: float = 0.5) -> bool:
    """
    Check if the sequence has good coverage.

    Args:
        sequence (str): The sequence to check.
        coverage (float): The coverage threshold. Default is 0.5.

    Returns:
        bool: True if the sequence has good coverage, False otherwise.
    """
    sequence_length_with_gaps = len(sequence)
    desired_sequence_length_without_gaps = sequence_length_with_gaps * coverage
    gap_count = sequence.count("-")
    sequence_length_without_gaps = sequence_length_with_gaps - gap_count
    if sequence_length_without_gaps >= desired_sequence_length_without_gaps:
        return True
    else:
        return False


def get_combined_a3m_records(a3m_files: List[str], tsv_files: List[str],
                             query_hashes: List[str]) -> Optional[A3M_Records]:
    """
    Get combined A3M records for all subunits.

    Args:
        a3m_files (List[str]): List of A3M file paths.
        tsv_files (List[str]): List of TSV file paths.
        query_hashes (List[str]): List of query hashes.

    Returns:
        A3M_Records or None: A named tuple of query hash to sequence and hits per query or None.
    """

    query_hash_to_tsv_records = dict()
    for tsv_file in tsv_files:
        tsv_file_without_extension = get_filename_without_extension(tsv_file)
        tsv_query_hash = tsv_file_without_extension.split("_", 1)[0]
        tsv_record_set = {f"{tsv_query_hash}:{line.rstrip().split()[0].strip()}" for line in open(tsv_file, "r")}
        query_hash_to_tsv_records[tsv_query_hash] = tsv_record_set

    query_hash_colon_hit_to_a3m = {}
    hit_count_per_query = []
    for query_hash in query_hashes:
        count = 0
        for a3m_file in a3m_files:
            a3m_file_without_extension = get_filename_without_extension(a3m_file)
            a3m_query_hash = a3m_file_without_extension.split("_", 1)[0]
            if query_hash == a3m_query_hash:
                tsv_records = query_hash_to_tsv_records.get(query_hash, None)
                if tsv_records is None:
                    print(f"-- Warning: No TSV records found for query hash: {query_hash}")
                    continue
                a3m_records = get_records(a3m_file)
                for hit_accession in a3m_records:
                    query_hash_colon_hit_accession = f"{query_hash}:{hit_accession}"
                    query_hash_colon_hit_to_a3m[query_hash_colon_hit_accession] = {
                        "is_ignorable": query_hash_colon_hit_accession not in tsv_records,
                        "a3m_sequence": a3m_records[hit_accession]
                    }
                    count += 1
        hit_count_per_query.append(count)

    if sum(hit_count_per_query) == 0:
        return None

    return A3M_Records(query_hash_to_seq=query_hash_colon_hit_to_a3m, hits_per_query=hit_count_per_query)


def make_alignment_pair(query_colon_hits: List[str], input_query_feats: Infile_feats,
                        a3m_seq_records: Dict[str, Dict]) -> List[Sequence]:
    """
    Create alignment pairs from query hits and input query features.

    Args:
        query_colon_hits (List[str]): List of query hits in the format 'query:hit'.
        input_query_feats (Infile_feats): Input query features.
        a3m_seq_records (Dict[str, Dict]): Dictionary of A3M records.

    Returns:
        List[Sequence]: A list of paired sequences
    """
    query_hashes = input_query_feats.chain_seq_hashes
    gappy_seq_limit = len(query_hashes) - 1

    list_of_top_hit_list = []
    for query_hash in query_hashes:
        top_hit_list = []
        for query_colon_hit in query_colon_hits:
            query_hash_value, hit_accession = query_colon_hit.split(":", 1)
            if query_hash == query_hash_value:
                if len(top_hit_list) == 0:
                    top_hit_list.append(query_colon_hit)
        if len(top_hit_list) == 1:
            list_of_top_hit_list.append(top_hit_list)
        else:
            list_of_top_hit_list.append(["-"])

    pairs_of_best_hits = list(itertools.product(*list_of_top_hit_list))

    list_of_alignment = []
    for potential_combination in pairs_of_best_hits:
        hit_combo = list(potential_combination)
        if hit_combo.count("-") < gappy_seq_limit:
            accession_combo = [""] * len(hit_combo)
            seq_combo = [""] * len(hit_combo)
            dummy_counter = 0
            for hit_index in range(len(hit_combo)):
                subunit_based_hit = hit_combo[hit_index]
                length = len(input_query_feats.chain_seqs[hit_index])
                if subunit_based_hit == "-":
                    dummy_counter += 1
                    created_sequence = "-" * length
                    accession_combo[hit_index] = "DUMMY"
                    seq_combo[hit_index] = created_sequence
                else:
                    hit_accession = subunit_based_hit.split(":", 1)[1]
                    sequence = a3m_seq_records[subunit_based_hit].get("a3m_sequence", "")
                    if has_good_coverage(sequence):
                        accession_combo[hit_index] = hit_accession
                        seq_combo[hit_index] = sequence
                    else:
                        dummy_counter += 1
                        created_sequence = "-" * length
                        accession_combo[hit_index] = "DUMMY"
                        hit_combo[hit_index] = "-"
                        seq_combo[hit_index] = created_sequence

            if dummy_counter < gappy_seq_limit:
                list_of_alignment.append(combine_sequences(accession_combo, seq_combo))

    return list_of_alignment


def introduce_gap_in_subunit(subunit_seq: List[str]) -> List[List[str]]:
    """
    Introduce gaps in subunit sequences.

    Args:
        subunit_seq (List[str]): List of subunit sequences.

    Returns:
        List[List[str]]: List of subunit sequences with gaps introduced.
    """
    placeholder = "-"
    # Generate combinations
    combinations = []
    for i in range(len(subunit_seq)):
        combination = subunit_seq.copy()
        for j in range(len(subunit_seq)):
            if i == j:
                pass
            else:
                combination[j] = placeholder * len(subunit_seq[j])
        combinations.append(combination)
    return combinations


def create_a3m_for_folding(summary_json: str, a3m_records: A3M_Records,
                           query_fasta: Infile_feats, output_a3m_file: str) -> None:
    """
    Create A3M files for folding.

    Args:
        summary_json (str): Path to the summary JSON file.
        a3m_records (A3M_Records): A named tuple of query hash to sequence and hits per query.
        query_fasta (Infile_feats): Input query features.
        output_a3m_file (str): Output file.

    Returns:
        None
    """
    chain_subunit_accessions = query_fasta.chain_accnrs
    a3m_sequence_records = a3m_records.query_hash_to_seq

    chain_subunit_mod_accessions = []
    for i in range(len(chain_subunit_accessions)):
        mod_accession = 101 + i
        chain_subunit_mod_accessions.append(mod_accession)

    chain_subunit_sequences = query_fasta.chain_seqs

    a3m_alignments = []
    a3m_alignment_hashes = set()

    processed_query = combine_sequences(chain_subunit_mod_accessions, chain_subunit_sequences)
    if processed_query.hash not in a3m_alignment_hashes:
        a3m_alignments.append(processed_query.fasta)
        a3m_alignment_hashes.add(processed_query.hash)

    gbk_to_query_colon_hits: Dict[str, List[str]] = load_json_file(summary_json)

    chain_seq_hashes = query_fasta.chain_seq_hashes
    is_query_a_hetero_complex = len(query_fasta.chain_seq_hashes) > 1

    gappy_seq_dict = {}
    if is_query_a_hetero_complex:
        for gbk in gbk_to_query_colon_hits:
            query_colon_hits = gbk_to_query_colon_hits[gbk]
            paired_alignment_lists = make_alignment_pair(query_colon_hits, query_fasta, a3m_sequence_records)
            for a3m_alignment in paired_alignment_lists:
                if a3m_alignment.hash not in a3m_alignment_hashes:
                    a3m_alignments.append(a3m_alignment.fasta)
                    a3m_alignment_hashes.add(a3m_alignment.hash)

        g = 0
        for gappy_subunit_sequences in introduce_gap_in_subunit(chain_subunit_sequences):
            combo_gappy = combine_gappy_sequences(chain_subunit_mod_accessions, gappy_subunit_sequences)
            gappy_seq_dict[g] = combo_gappy
            g += 1

        for j in range(len(chain_seq_hashes)):
            empty_chain_subunit_seq_list = query_fasta.empty_subunits.copy()
            for query_hash_colon_hit in a3m_sequence_records:
                query_hash_value, hit_accession = query_hash_colon_hit.split(":", 1)
                is_ignorable = a3m_sequence_records[query_hash_colon_hit].get("is_ignorable", False)
                if query_hash_value == chain_seq_hashes[j] and not is_ignorable:
                    empty_chain_subunit_seq_list[j] = a3m_sequence_records[query_hash_colon_hit].get("a3m_sequence", "")
                    combo_unpaired = combine_sequences([hit_accession], empty_chain_subunit_seq_list)
                    if gappy_seq_dict[j].hash not in a3m_alignment_hashes:
                        a3m_alignments.append(gappy_seq_dict[j].fasta)
                        a3m_alignment_hashes.add(gappy_seq_dict[j].hash)
                    if combo_unpaired.hash not in a3m_alignment_hashes:
                        a3m_alignments.append(combo_unpaired.fasta)
                        a3m_alignment_hashes.add(combo_unpaired.hash)

    else:
        for query_hash_colon_hit in a3m_sequence_records:
            query_hash_value, hit_accession = query_hash_colon_hit.split(":", 1)
            is_ignorable = a3m_sequence_records[query_hash_colon_hit].get("is_ignorable", False)
            if query_hash_value == chain_seq_hashes[0] and not is_ignorable:
                sequence = a3m_sequence_records[query_hash_colon_hit].get("a3m_sequence", "")
                fasta_sequence = make_hash_fasta_sequence(hit_accession, sequence)
                if fasta_sequence.hash not in a3m_alignment_hashes:
                    a3m_alignments.append(fasta_sequence.fasta)
                    a3m_alignment_hashes.add(fasta_sequence.hash)

    with open(output_a3m_file, "w") as con_out:
        # noinspection PyTypeChecker
        print(query_fasta.a3m_header, file=con_out)
        for seq_aln in a3m_alignments:
            # noinspection PyTypeChecker
            print(seq_aln.rstrip(), file=con_out)
    end_time = current_time()
    print(f"-- {end_time} > Completed filtering MSA, check output at: \n\t'{output_a3m_file}'\n")
    return None
