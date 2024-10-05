import os
import itertools
from typing import List, Dict, Set
from .util import get_filename_without_extension, join_list_elements_by_character, current_time
from .sequence import Sequence, Infile_feats, get_alignment_records_from_a3m_file, combine_sequences, \
    combine_gappy_sequences, make_fasta_sequence
from .json import load_json_file
from .execute import run_jobs_in_parallel


def run_jackhmmer(fasta_files: list, database_fasta: str, cpu: int, out_path: str, perl_path: str) -> None:
    jackhmmer_commands = []
    sto_to_a3m_commands = []
    n_threads = max(1, round(cpu / len(fasta_files)))
    for fasta_file in fasta_files:
        fasta_basename = os.path.basename(fasta_file)
        sto_file_name = f"{os.path.splitext(fasta_basename)[0]}.sto"
        sto_file_path = os.path.join(os.path.abspath(out_path), sto_file_name)
        a3m_file_name = f"{os.path.splitext(fasta_basename)[0]}.a3m"
        a3m_file_path = os.path.join(os.path.abspath(out_path), a3m_file_name)
        sto_command = ("jackhmmer --noali --F1 0.0005 --F2 0.00005 --F3 0.0000005 --incE 0.0001 -E 0.0001 -N 1 "
                       "-o /dev/null --cpu %s -A %s %s %s" % (n_threads, sto_file_path, fasta_file, database_fasta))
        jackhmmer_commands.append(sto_command)
        sto_to_a3m_command = ("%s extra/reformat.pl sto a3m %s %s > /dev/null" % (perl_path, sto_file_path,
                                                                                  a3m_file_path))
        sto_to_a3m_commands.append(sto_to_a3m_command)

    run_jobs_in_parallel(cpu, jackhmmer_commands, "Homology searching")
    run_jobs_in_parallel(cpu, sto_to_a3m_commands, "Alignment formatting")


def has_good_coverage(sequence: str, coverage: float = 0.5) -> bool:
    sequence_length = len(sequence)
    coverage_threshold = sequence_length * coverage
    gap_count = sequence.count("-")
    if gap_count < coverage_threshold:
        return True
    else:
        return False


def get_query_to_a3m_records(a3m_files: List[str], query_hashes: List[str]) -> Dict[str, str]:
    query_hash_colon_hit_to_a3m = {}
    for query_hash in query_hashes:
        for a3m_file in a3m_files:
            a3m_file_without_extension = get_filename_without_extension(a3m_file)
            a3m_query_hash = a3m_file_without_extension.split("_", 1)[0]
            if query_hash == a3m_query_hash:
                a3m_records = get_alignment_records_from_a3m_file(a3m_file)
                for hit_accession in a3m_records:
                    query_hash_colon_hit_accession = f"{query_hash}:{hit_accession}"
                    query_hash_colon_hit_to_a3m[query_hash_colon_hit_accession] = a3m_records[hit_accession]
    return query_hash_colon_hit_to_a3m


def make_alignment_pair(query_colon_hits: List[str], input_query_feats: Infile_feats,
                               a3m_records: Dict[str, str]) -> tuple[List[Sequence], List[str]]:

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
    list_of_a3m_records_keys_paired = []
    for potential_combination in pairs_of_best_hits:
        hit_combo = list(potential_combination)
        if hit_combo.count("-") < gappy_seq_limit:
            accession_combo = [""] * len(hit_combo)
            seq_combo = [""] * len(hit_combo)
            dummy_counter = 0
            for hit_index in range(len(hit_combo)):
                subunit_based_hit = hit_combo[hit_index]
                if subunit_based_hit == "-":
                    dummy_counter += 1
                    length = len(input_query_feats.chain_seqs[hit_index])
                    created_sequence = "-" * length
                    accession_combo[hit_index] = "DUMMY"
                    seq_combo[hit_index] = created_sequence
                else:
                    hit_accession = subunit_based_hit.split(":", 1)[1]
                    sequence = a3m_records[subunit_based_hit]
                    if has_good_coverage(sequence):
                        accession_combo[hit_index] = hit_accession
                        seq_combo[hit_index] = sequence
                    else:
                        print("Warning: Sequence coverage is less than 50% for {}".format(hit_accession))
                        dummy_counter += 1
                        created_sequence = "-" * len(sequence)
                        accession_combo[hit_index] = "DUMMY"
                        hit_combo[hit_index] = "-"
                        seq_combo[hit_index] = created_sequence

            if dummy_counter < gappy_seq_limit:
                list_of_a3m_records_keys_paired.extend(hit_combo)
                list_of_alignment.append(combine_sequences(accession_combo, seq_combo))

    return list_of_alignment, list_of_a3m_records_keys_paired


def introduce_gap_in_subunit(subunit_seq: List[str]) -> List[List[str]]:
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


def create_a3m_for_folding(summary_json: str, a3m_records: Dict[str, str],
                           query_fasta: Infile_feats, out_path: str) -> None:
    chain_subunit_accessions = query_fasta.chain_accnrs

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

    set_of_a3m_records_keys_paired: Set[str] = set()
    gappy_seq_dict = {}

    if is_query_a_hetero_complex:
        for gbk in gbk_to_query_colon_hits:
            query_colon_hits = gbk_to_query_colon_hits[gbk]
            paired_alignment_lists, list_of_a3m_records_keys_paired \
                = make_alignment_pair(query_colon_hits, query_fasta, a3m_records)
            for a3m_alignment in paired_alignment_lists:
                if a3m_alignment.hash not in a3m_alignment_hashes:
                    a3m_alignments.append(a3m_alignment.fasta)
                    a3m_alignment_hashes.add(a3m_alignment.hash)
                    set_of_a3m_records_keys_paired.update(list_of_a3m_records_keys_paired)

        g = 0
        for gappy_subunit_sequences in introduce_gap_in_subunit(chain_subunit_sequences):
            combo_gappy = combine_gappy_sequences(chain_subunit_mod_accessions, gappy_subunit_sequences)
            gappy_seq_dict[g] = combo_gappy
            g += 1

        for j in range(len(chain_seq_hashes)):
            empty_chain_subunit_seq_list = query_fasta.empty_subunits.copy()
            for query_hash_colon_hit in a3m_records:
                if query_hash_colon_hit not in set_of_a3m_records_keys_paired:
                    query_hash_value, hit_accession = query_hash_colon_hit.split(":", 1)
                    if query_hash_value == chain_seq_hashes[j]:
                        empty_chain_subunit_seq_list[j] = a3m_records[query_hash_colon_hit]
                        combo_unpaired = combine_sequences([hit_accession], empty_chain_subunit_seq_list)
                        if gappy_seq_dict[j].hash not in a3m_alignment_hashes:
                            a3m_alignments.append(gappy_seq_dict[j].fasta)
                            a3m_alignment_hashes.add(gappy_seq_dict[j].hash)
                        if combo_unpaired.hash not in a3m_alignment_hashes:
                            a3m_alignments.append(combo_unpaired.fasta)
                            a3m_alignment_hashes.add(combo_unpaired.hash)
    else:
        for query_hash_colon_hit in a3m_records:
            query_hash_value, hit_accession = query_hash_colon_hit.split(":", 1)
            if query_hash_value == chain_seq_hashes[0]:
                sequence = a3m_records[query_hash_colon_hit]
                fasta_sequence = make_fasta_sequence(hit_accession, sequence)
                if fasta_sequence.hash not in a3m_alignment_hashes:
                    a3m_alignments.append(fasta_sequence.fasta)
                    a3m_alignment_hashes.add(fasta_sequence.hash)

    concatenated_filename = join_list_elements_by_character(query_fasta.accnrs, "-")
    concat_filepath = os.path.join(out_path, f"{concatenated_filename}.a3m")
    with open(concat_filepath, "w") as con_out:
        # noinspection PyTypeChecker
        print(query_fasta.a3m_header, file=con_out)
        for seq_aln in a3m_alignments:
            # noinspection PyTypeChecker
            print(seq_aln.rstrip(), file=con_out)
    end_time = current_time()
    print(f"-- {end_time} > Finished preparation for structure prediction, check: \n\t{concat_filepath}\n")
    return None
