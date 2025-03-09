import argparse


def align_sequence_to_gapless_query(sequence: str, query_sequence: str) -> str:
    """
    Aligns a sequence to a gapless query sequence by removing gaps.
    Args:
        sequence (str): The sequence to align.
        query_sequence (str): The gapless query sequence.
    Returns:
        str: The aligned sequence with gaps removed.
    """
    if len(sequence) != len(query_sequence):
        raise ValueError(f"The sequence ({len(sequence)}) and the query "
                         f"sequence ({len(query_sequence)}) don't have the same length.")

    output = []
    for residue_index in range(len(sequence)):
        query_residue = query_sequence[residue_index]
        residue = sequence[residue_index]
        if query_residue != '-':
            output.append(residue)
        elif residue == '-':
            continue
        else:
            output.append(residue.lower())

    return ''.join(output)


def convert_stockholm(
    stockholm_path: str,
    a3m_file: str,
    fas_file: str | None = None,
    max_sequences: int | None = None,
    remove_first_row_gaps: bool = True,
    line_width: int | None = None,
) -> None:
    """Converts MSA in Stockholm format to the A3M format."""
    descriptions = {}
    sequences = {}

    if line_width is not None and line_width <= 0:
        raise ValueError('line_width must be > 0 or None')

    stockholm = open(stockholm_path, 'r')

    for line in stockholm:
        #reached_max_sequences = max_sequences and len(sequences) >= max_sequences
        line = line.strip()
        # Ignore blank lines, markup and end symbols - remainder are alignment
        # sequence parts.
        if not line or line.startswith(('#', '//')):
            continue
        seq_name, aligned_seq = line.split(maxsplit=1)
        if seq_name not in sequences:
            sequences[seq_name] = ''
        sequences[seq_name] += aligned_seq

    if not sequences:
        print(f"\n-- Warning: No sequences found in '{stockholm_path}'\n")
        return None

    stockholm.seek(0)
    for line in stockholm:
        line = line.strip()
        if line[:4] == '#=GS':
            # Description row - example format is:
            # #=GS UniRef90_Q9H5Z4/4-78            DE [subseq from] cDNA: FLJ22755 ...
            columns = line.split(maxsplit=3)
            seq_name, feature = columns[1:3]
            value = columns[3] if len(columns) == 4 else ''
            if feature != 'DE':
                continue
            #if reached_max_sequences and seq_name not in sequences:
            #    continue
            descriptions[seq_name] = value.split(":")[0].split(" ")[-1]
            if len(descriptions) == len(sequences):
                break

    assert len(descriptions) <= len(sequences)

    # Convert sto format to a3m line by line
    a3m_sequences = {}
    # query_sequence is assumed to be the first sequence
    query_sequence = next(iter(sequences.values()))

    seq_hash_set = set()
    for seq_name, sto_sequence in sequences.items():
        seq_hash = descriptions.get(seq_name, '')
        if remove_first_row_gaps:
            if seq_hash in seq_hash_set:
                continue
            a3m_sequences[seq_name] = align_sequence_to_gapless_query(
                sequence=sto_sequence, query_sequence=query_sequence).replace('.', '')
            seq_hash_set.add(seq_hash)
        else:
            if seq_hash in seq_hash_set:
                continue
            a3m_sequences[seq_name] = sto_sequence.replace('.', '')
            seq_hash_set.add(seq_hash)

    a3m_seq_chunks = []
    fasta_seq_chunks = []

    seq_count = 0
    for seq_name, a3m_seq in a3m_sequences.items():
        seq_count += 1
        fasta_seq = a3m_seq.replace('-', '')
        if max_sequences and seq_count > max_sequences:
            break
        a3m_seq_chunks.append(f'>{seq_name.strip()}\t{descriptions.get(seq_name, "").strip()}')
        fasta_seq_chunks.append(f'>{seq_name.strip()}\t{descriptions.get(seq_name, "").strip()}')
        if line_width:
            a3m_seq_chunks.extend(
                a3m_seq[i: line_width + i] for i in range(0, len(a3m_seq), line_width)
            )
            a3m_seq_chunks.extend(
                fasta_seq[i: line_width + i] for i in range(0, len(a3m_seq), line_width)
            )
        else:
            a3m_seq_chunks.append(a3m_seq)
            fasta_seq_chunks.append(fasta_seq)

    with open(a3m_file, 'w') as a3m_out:
        a3m_out.write('\n'.join(a3m_seq_chunks) + '\n')

    if fas_file is None:
        return

    with open(fas_file, 'w') as fas_out:
        fas_out.write('\n'.join(fasta_seq_chunks) + '\n')

    return


def main():
    parser = argparse.ArgumentParser(description='Convert Stockholm format to A3M and FASTA format.')
    parser.add_argument('stockholm_path', type=str, help='Path to the Stockholm file.')
    parser.add_argument('-oa', '--output_a3m', type=str, required=True,
                        help='Output A3M file path.')
    parser.add_argument('-of', '--output_fas', type=str,
                        help='Output Fasta file path.')
    parser.add_argument('-m', '--max_sequences', type=int, default=None,
                        help='Maximum number of sequences to include.')
    parser.add_argument('-r', '--remove_first_row_gaps', type=bool, default=True,
                        help='Remove gaps in the first row.')
    parser.add_argument('-l', '--line_width', type=int, default=None,
                        help='Line width for the output A3M file.')

    args = parser.parse_args()

    convert_stockholm(
        stockholm_path=args.stockholm_path,
        a3m_file=args.output_a3m,
        fas_file=args.output_fas,
        max_sequences=args.max_sequences,
        remove_first_row_gaps=args.remove_first_row_gaps,
        line_width=args.line_width
    )


if __name__ == '__main__':
    main()
