import os
from flashfold.utils import run_jobs_in_parallel, run_single_job, get_sequence_length_from_single_fasta


min_jackhmmer_hits = 1000


def run_jackhmmer(fasta_files: list, database_fasta: str, provided_cpu: int, out_path: str) -> None:
    """
    Run jackhmmer for homology searching.

    Args:
        fasta_files (list): List of paths to FASTA files.
        database_fasta (str): Path to the database FASTA file.
        provided_cpu (int): Number of CPU cores to use.
        out_path (str): Output directory path.

    Returns:
        None
    """
    current_file_dir = os.path.dirname(os.path.abspath(__file__))
    reformat_script = os.path.join(current_file_dir, "reformat_msa.py")
    deduplicate_script = os.path.join(current_file_dir, "deduplicate_msa.py")
    combined_commands = []
    cpu_per_job = min(8, provided_cpu)
    for fasta_file in fasta_files:
        fasta_seq_length = get_sequence_length_from_single_fasta(fasta_file)
        fasta_basename = os.path.basename(fasta_file)
        sto_file_name = f"{os.path.splitext(fasta_basename)[0]}.sto"
        a3m_file_name = f"{os.path.splitext(fasta_basename)[0]}.a3m"
        fas_file_name = f"{os.path.splitext(fasta_basename)[0]}.fas"
        diverse_hits_file_name = f"{os.path.splitext(fasta_basename)[0]}.tsv"
        sto_file_path = os.path.join(os.path.abspath(out_path), sto_file_name)
        a3m_file_path = os.path.join(os.path.abspath(out_path), a3m_file_name)
        fas_file_path = os.path.join(os.path.abspath(out_path), fas_file_name)
        diverse_hits_file_path = os.path.join(os.path.abspath(out_path), diverse_hits_file_name)
        sto_command = ("jackhmmer --noali --F1 0.0005 --F2 0.00005 --F3 0.0000005 --incE 0.0001 -E 0.0001 -N 1 "
                       "-o /dev/null --cpu %s -A %s %s %s" % (cpu_per_job, sto_file_path, fasta_file, database_fasta))
        reformat_sto_command = ("python3 %s %s --output_a3m %s --output_fas %s"
                                % (reformat_script, sto_file_path, a3m_file_path, fas_file_path))
        deduplicate_command = ("python3 %s %s --query_length %s --min_sequences %s --output %s --threads %s"
                               % (deduplicate_script, fas_file_path, fasta_seq_length, min_jackhmmer_hits,
                                   diverse_hits_file_path, cpu_per_job))
        combined_command = f"{sto_command} && {reformat_sto_command} && {deduplicate_command}"
        combined_commands.append(combined_command)
    run_jobs_in_parallel(provided_cpu, cpu_per_job, combined_commands, "Homology searching")
    return


def run_af3tools(a3m_file_path: str, json_file_path: str) -> None:
    """
    Run af3tools for creating json file from msa for AF3 input.

    Args:
        a3m_file_path (str): Path to the A3M file.
        json_file_path (str): Path to the JSON file.

    Returns:
        None
    """
    current_file_dir = os.path.dirname(os.path.abspath(__file__))
    af3tools_script = os.path.join(current_file_dir, "msa_to_json.py")
    af3tools_command = f"python3 {af3tools_script} -i {a3m_file_path} -o {json_file_path}"
    run_single_job(af3tools_command, "JSON file creation")
    return


