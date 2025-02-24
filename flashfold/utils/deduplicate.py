import argparse
import os
import shutil
import subprocess
from typing import Dict
from Bio import SeqIO


Memory_CD_hit = 2000


def get_fasta_records(file_path: str) -> Dict[str, str]:
    """
    Get Fasta records from a file.
    Args:
        file_path: Path to the FASTA file.

    Returns:
        Dict: A dictionary containing the accession and sequence
    """
    parsed_fasta = SeqIO.parse(file_path, "fasta")
    records = {}
    for record in parsed_fasta:
        accession = str(record.description).strip()
        sequence = str(record.seq)
        # Accession not empty validity check
        if accession == "" or sequence == "":
            pass
        else:
            records[accession] = sequence
    return records


def get_min_num_of_diverse_hits(fasta_path: str, query_length: int, outfile: str, min_sequences: int = 1000,
                                threads: int = 8) -> None:
    """
    Get the most diverse hits from JackHmmer parsed FASTA output.
    :param fasta_path: Path to the FASTA file.
    :param query_length: Length of the query sequence.
    :param outfile: Output file path.
    :param min_sequences: Minimum number of sequences to include.
    :param threads: Number of threads.
    :return: None
    """

    init_seq_records = get_fasta_records(fasta_path)

    outfile_base_name = os.path.basename(outfile)

    if len(init_seq_records) <= min_sequences:
        with open(outfile, 'w') as out_file:
            for record in init_seq_records:
                out_file.write(f"{record}\n")
        return

    # Get the most diverse hits
    root_out_dir = os.path.dirname(outfile)

    for i in range(1, 10):
        identity = f"{i/10:.2f}"
        out_dir = os.path.join(root_out_dir, f"{outfile_base_name[:-4]}_cd_hit_{identity}")
        out_fas = os.path.join(out_dir, f"{outfile_base_name[:-4]}.fas")
        os.makedirs(out_dir, exist_ok=True)
        cd_hit_command = ("cd-hit -i %s -o %s -c %s -T %s -M %s" %
                          (fasta_path, out_fas, identity, threads, Memory_CD_hit))
        try:
            subprocess.run(cd_hit_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as _:
            continue

        seq_records = get_fasta_records(out_fas)

        if len(seq_records) >= min_sequences:
            with open(outfile, 'w') as out_file:
                for record in seq_records:
                    out_file.write(f"{record}\n")
            return
        elif i == 9:
            with open(outfile, 'w') as out_file:
                s = -1
                for record in init_seq_records:
                    s += 1
                    if s < 1000:
                        out_file.write(f"{record}\n")
                    else:
                        sequence_length = len(init_seq_records[record])
                        # Check if the sequence length is at least 50% of the query length
                        if round(sequence_length/query_length, 2) >= 0.5:
                            out_file.write(f"{record}\n")
            return
        else:
            continue
    return


def main():
    parser = argparse.ArgumentParser(description='Gets diverse hits from JackHmmer parsed FASTA output.')
    parser.add_argument('fasta_path', type=str, help='Path to the FASTA file.')
    parser.add_argument('-ql', '--query_length', type=int, required=True,
                        help='Length of the query sequence.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output file.')
    parser.add_argument('-m', '--min_sequences', type=int, default=1000,
                        help='Minimum number of sequences to include.')
    parser.add_argument('-t', '--threads', type=int, default=8,
                        help='Number of threads.')
    args = parser.parse_args()

    get_min_num_of_diverse_hits(
        fasta_path = args.fasta_path,
        query_length = args.query_length,
        outfile = args.output,
        min_sequences = args.min_sequences,
        threads = args.threads
    )


if __name__ == '__main__':
    main()
