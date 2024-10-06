import os
import re
import sys
from collections import namedtuple
import numpy as np
import pandas as pd
from typing import List, Optional, Tuple

# noinspection PyPackageRequirements
from Bio.PDB.Structure import Structure

# noinspection PyPackageRequirements
from Bio.PDB import PDBParser

from .util import current_time, get_filename_without_extension
from .json import load_json_file


PdockQ_Calc_Out = namedtuple('PdockQ_Calc_Out', ["num_res_list", "if_plddt_list", "plddt_list",
                                                 "num_diso_list", "pdockq_list", "length_list"])


def extract_plddt(log_line: str) -> Optional[float]:
    """
    Extracts the pLDDT value from a log line.

    :param log_line: A string containing the log line.
    :return: The pLDDT value as a float if found, otherwise None.
    """
    pattern = r"pLDDT=(\d+(\.\d+)?)"
    match = re.search(pattern, log_line)
    if match:
        return float(match.group(1))
    else:
        return None


def extract_ptm(log_line: str) -> Optional[float]:
    """
    Extracts the pTM value from a log line.

    :param log_line: A string containing the log line.
    :return: The pTM value as a float if found, otherwise None.
    """
    pattern = r"pTM=(\d+(\.\d+)?)"
    match = re.search(pattern, log_line)
    if match:
        return float(match.group(1))
    else:
        return None


def extract_plddt_ptm(log_line: str) -> Optional[Tuple[str, float, float]]:
    """
    Extracts the name, pLDDT, and pTM values from a log line.

    :param log_line: A string containing the log line.
    :return: A tuple containing the name, pLDDT, and pTM values if found, otherwise None.
    """
    if "rank_" not in log_line:
        return None

    split_line = log_line.rstrip().split()
    name = split_line[2]
    plddt = extract_plddt(log_line)
    ptm = extract_ptm(log_line)
    if plddt and ptm:
        return name, plddt, ptm
    else:
        return None


def load_structure(pdb_file_path: str) -> Structure:
    """
    Loads a PDB structure from a file and returns a Structure object.

    :param pdb_file_path: Path to the PDB file.
    :return: A Structure object representing the PDB structure.
    """
    structure_id = get_filename_without_extension(pdb_file_path)
    bio_parser = PDBParser()
    return bio_parser.get_structure(structure_id, pdb_file_path)


def get_popt_from_cutoff(cutoff: float) -> List[float]:
    """
    Returns the popt values based on the cutoff value.

    :param cutoff: The cutoff value.
    :return: A list of popt values.
    """
    if cutoff <= 5:
        return [6.96234405e-01, 2.35483775e+02, 2.25322970e-02, 2.88445245e-02]
    elif cutoff <= 6:
        return [7.02605033e-01, 2.91749822e+02, 2.70621128e-02, 2.25416051e-02]
    elif cutoff <= 7:
        return [7.06385097e-01, 3.32456259e+02, 2.97005237e-02, 2.24488132e-02]
    elif cutoff <= 8:
        return [7.18442739e-01, 3.60791204e+02, 3.01635944e-02, 2.04076969e-02]
    elif cutoff <= 9:
        return [7.23328534e-01, 3.80036094e+02, 3.06316084e-02, 1.98471192e-02]
    elif cutoff <= 10:
        return [7.20293782e-01, 3.95627723e+02, 3.15235037e-02, 2.37304238e-02]
    elif cutoff <= 11:
        return [7.22015998e-01, 4.09095024e+02, 3.11905555e-02, 2.59467513e-02]
    else:
        return [7.20555781e-01, 4.21033584e+02, 3.09024241e-02, 2.88659629e-02]


def calc_pdockq(structure_bio: Structure, cutoff_score: float, diso_cut_1=50, diso_cut_2=70,
                diso_cut_3=90) -> PdockQ_Calc_Out:
    """
    Calculates the pDockQ score for a given structure.

    :param structure_bio: A Structure object representing the PDB structure.
    :param cutoff_score: The cutoff score for calculating the pDockQ score.
    :param diso_cut_1: The first disorder cutoff value.
    :param diso_cut_2: The second disorder cutoff value.
    :param diso_cut_3: The third disorder cutoff value.
    :return: A namedtuple containing the pDockQ calculation output.
    """
    popt = get_popt_from_cutoff(cutoff_score)
    if2 = None
    i = 0
    tiny = 1.e-20
    chains = []
    num_res_list = []
    if_plddt_list = []
    plddt_list = []
    num_diso_list = []
    pdockq_list = []
    length_list = []
    max_len = 10000
    for chain in structure_bio.get_chains():
        chains += [chain]
        i += 1

    for c in range(len(chains)):
        interface_residues1 = []
        interface_residues2 = []
        k = tiny
        b = 0
        for d in range(len(chains)):
            if c == d:
                continue
            for res1 in chains[c]:
                for res2 in chains[d]:
                    test = False
                    for i in res1:
                        if test:
                            break
                        for j in res2:
                            dist = np.linalg.norm(i.coord - j.coord)
                            if dist < cutoff_score:
                                interface_residues1.append(res1.id[1])
                                interface_residues2.append(res2.id[1] + max_len * d)
                                test = True
                                break
                            elif dist > 2 * cutoff_score:
                                test = True
                                break
            if2 = np.unique(interface_residues2)
            for res in chains[d]:
                if res.id[1] in if2:
                    b += res['CA'].get_bfactor()
                    k += 1
        if1 = np.unique(interface_residues1)
        b1 = 0
        i1 = 0
        num_diso1 = [0, 0, 0, 0]
        for res in chains[c]:
            b1 += res['CA'].get_bfactor()
            i1 += 1
            if res['CA'].get_bfactor() > diso_cut_3:
                num_diso1[0] += 1
            elif res['CA'].get_bfactor() > diso_cut_2:
                num_diso1[1] += 1
            elif res['CA'].get_bfactor() > diso_cut_1:
                num_diso1[2] += 1
            else:
                num_diso1[3] += 1
            if res.id[1] in if1:
                b += res['CA'].get_bfactor()
                k += 1
        length_list += [i1]
        num_res_list += [if1.shape[0] + if2.shape[0]]
        if_plddt_list += [b / k]
        plddt_list += [b1 / i1]
        num_diso_list += [num_diso1]
        pdockq_list += [sigmoid(np.log(num_res_list[c] + tiny) * if_plddt_list[c], *popt)]

    return PdockQ_Calc_Out(num_res_list, if_plddt_list, plddt_list, num_diso_list, pdockq_list, length_list)


def sigmoid(x: float, l_cap: float, x0: float, k: float, b: float) -> float:
    """
    Computes the sigmoid function.

    :param x: The input value.
    :param l_cap: The maximum value of the sigmoid function.
    :param x0: The x-value of the sigmoid's midpoint.
    :param k: The steepness of the sigmoid curve.
    :param b: The y-offset of the sigmoid function.
    :return: The computed sigmoid value.
    """
    y = l_cap / (1 + np.exp(-k * (x - x0))) + b
    return y


def generate_score_matrix(path_to_predicted_structure: str, cutoff: float, is_monomer: bool) -> None:
    """
    Generates a score matrix for the predicted structure.

    :param path_to_predicted_structure: Path to the directory containing the predicted structure files.
    :param cutoff: Cutoff value for calculating scores.
    :param is_monomer: Boolean indicating if the structure is a monomer.
    :return: None
    """
    score_matrix_filepath = os.path.join(path_to_predicted_structure, "score.tsv")

    if is_monomer:
        log_file = os.path.join(path_to_predicted_structure, "log.txt")
        with open(log_file, "r") as log_in, open(score_matrix_filepath, "w") as score_out:
            print("name", "pLDDT", "pTM", sep="\t", file=score_out)
            for line in log_in:
                scores = extract_plddt_ptm(line)
                if scores:
                    name, plddt, ptm = scores
                    print(name, plddt, ptm, sep="\t", file=score_out)
                else:
                    pass
        return None

    all_files = os.listdir(path_to_predicted_structure)
    pdb_files = [f for f in all_files if f.endswith(".pdb") and "unrelaxed" in f]

    names = []
    mean_pdockq_scores = []
    min_pdockq_scores = []
    all_pdockq2 = []
    iptm = []
    iptm_ptm = []
    plddt_score = []
    for predicted_pdb in pdb_files:
        pdb_file = predicted_pdb.split("/")[-1]
        if 'rank' not in pdb_file:
            print(f"Skipping {pdb_file} as it does not contain 'rank' in the file name.")
            return

        json_file_name = predicted_pdb.replace(".pdb", ".json").replace("unrelaxed", "scores")
        json_file_path = os.path.join(path_to_predicted_structure, json_file_name)
        names.append(predicted_pdb)
        json_file = load_json_file(json_file_path)

        iptm_score = json_file['iptm']
        iptm.append(iptm_score)
        iptm_ptm_score = 0.8 * json_file["iptm"] + 0.2 * json_file["ptm"]
        iptm_ptm.append(iptm_ptm_score)

        mean_plddt = np.mean(json_file["plddt"])
        plddt_score.append(mean_plddt)

        pdb_file_path = os.path.join(path_to_predicted_structure, predicted_pdb)
        structure = load_structure(pdb_file_path)
        num_res, if_plddt, plddt, num_diso, pdockq, length = calc_pdockq(structure, cutoff)

        chain_number = len(length)
        if chain_number < 2:
            print(f'Warning: pDockQ2 score is not calculated because predicted pdb file contains {chain_number} '
                  f'chain. FoldFlash currently offers calculation when the pdb file contains at least two chains.\n')
            sys.exit()

        mean_pdockq = sum(pdockq) / len(pdockq)
        min_pdockq = min(pdockq)
        mean_pdockq_scores.append(mean_pdockq)
        min_pdockq_scores.append(min_pdockq)

        chain_counter = 0
        pdockq_scores_per_model = []
        for chain_unit in structure.get_chains():
            chain_id = chain_unit.get_id()
            pdockq_score_with_chain = f"{chain_id}: {pdockq[chain_counter]}"
            pdockq_scores_per_model.append(pdockq_score_with_chain)
            chain_counter += 1

        all_pdockq2.append(pdockq_scores_per_model)

    score_dict = {'name': names, 'mean_pDockQ2': mean_pdockq_scores, 'min_pDockQ2': min_pdockq_scores,
                  'ipTM': iptm, 'ipTM+pTM': iptm_ptm, 'mplDDT': plddt_score, 'pDockQ2_per_chain': all_pdockq2}
    collect_data = pd.DataFrame(score_dict)
    sorted_collected_data = collect_data.sort_values(by='name')
    sorted_collected_data.to_csv(score_matrix_filepath, sep='\t', index=False)
    print(f"-- {current_time()} > Completed pDockQ2 score calculation for predicted complex\n")
    return None

