# pDockQ2 Code is collected from: https://gitlab.com/ElofssonLab/huintaf2/-/blob/main/bin/pDockQ2.py
# And changed accordingly to fit the current project.

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


def generate_score_matrix(path_to_predicted_structure: str, cutoff: Optional[float], is_monomer: bool,
                          is_actifptm: bool) -> None:
    """
    Generates a score matrix for the predicted structure.

    :param path_to_predicted_structure: Path to the directory containing the predicted structure files.
    :param cutoff: Cutoff value for calculating scores.
    :param is_monomer: Boolean indicating if the structure is a monomer.
    :param is_actifptm: Boolean indicating if the structure contains active site PTM.
    :return: None
    """
    score_matrix_filepath = os.path.join(path_to_predicted_structure, "score.tsv")

    all_files = os.listdir(path_to_predicted_structure)
    pdb_files = [f for f in all_files if f.endswith(".pdb") and "_unrelaxed_rank_" in f]

    if is_monomer:
        if cutoff:
            print(f"\n-- Warning: Cutoff {cutoff} is not applicable for monomers, "
                  f"but will be used for pDockQ2 scoring for multimers (if provided any).\n")
            pass

        names = []
        plddt_score = []
        ptm = []

        for predicted_pdb in pdb_files:
            json_file_name = (predicted_pdb.replace(".pdb", ".json").
                              replace("_unrelaxed_", "_scores_"))
            json_file_path = os.path.join(path_to_predicted_structure, json_file_name)
            names.append(predicted_pdb)
            json_file = load_json_file(json_file_path)
            ptm.append(json_file['ptm'])
            mean_plddt = np.mean(json_file["plddt"])
            plddt_score.append(mean_plddt)

        score_dict = {'name': names, 'pLDDT': plddt_score, 'pTM': ptm}
        collect_data = pd.DataFrame(score_dict)
        sorted_collected_data = collect_data.sort_values(by='name')
        sorted_collected_data.to_csv(score_matrix_filepath, sep='\t', index=False)
        return None

    if not is_monomer:
        cutoff_score = cutoff if cutoff else 10.0

        names = []
        mean_pdockq_scores = []
        min_pdockq_scores = []
        all_pdockq2 = []
        ptm = []
        iptm = []
        actifptm = []
        iptm_ptm = []
        actifptm_ptm = []
        plddt_score = []
        pairwise_actifptm = []
        pairwise_iptm = []
        per_chain_ptm = []

        for predicted_pdb in pdb_files:
            json_file_name = (predicted_pdb.replace(".pdb", ".json").
                              replace("_unrelaxed_", "_scores_"))
            json_file_path = os.path.join(path_to_predicted_structure, json_file_name)

            names.append(predicted_pdb)
            json_file = load_json_file(json_file_path)

            ptm.append(json_file.get('ptm', 0))
            iptm.append(json_file.get('iptm', 0))
            actifptm.append(json_file.get('actifptm', 0))
            iptm_ptm.append(0.8 * json_file.get("iptm", 0) + 0.2 * json_file.get("ptm", 0))
            actifptm_ptm.append(0.8 * json_file.get("actifptm", 0) + 0.2 * json_file.get("ptm", 0))
            plddt_score.append(np.mean(json_file["plddt"]))
            pairwise_actifptm.append(json_file.get("pairwise_actifptm", 0))
            pairwise_iptm.append(json_file.get("pairwise_iptm", 0))
            per_chain_ptm.append(json_file.get("per_chain_ptm", 0))

            pdb_file_path = os.path.join(path_to_predicted_structure, predicted_pdb)
            structure = load_structure(pdb_file_path)
            num_res, if_plddt, plddt, num_diso, pdockq, length = calc_pdockq(structure, cutoff_score)

            chain_number = len(length)
            if chain_number < 2:
                print(f'Warning: pDockQ2 score is not calculated because predicted pdb file contains {chain_number} '
                      f'chain. FlashFold currently offers calculation when a pdb file contains at least two chains.\n')
                sys.exit()

            mean_pdockq = sum(pdockq) / len(pdockq)
            min_pdockq = min(pdockq)
            mean_pdockq_scores.append(mean_pdockq)
            min_pdockq_scores.append(min_pdockq)

            chain_counter = 0
            pdockq_scores_per_model = {}
            for chain_unit in structure.get_chains():
                chain_id = chain_unit.get_id()
                pdockq_score_with_chain = pdockq[chain_counter]
                pdockq_scores_per_model[chain_id] = pdockq_score_with_chain
                chain_counter += 1

            all_pdockq2.append(pdockq_scores_per_model)

        score_dict = {
            'name': names,
            'pLDDT': plddt_score,
            'pTM': ptm,
            'ipTM': iptm,
            'ipTM+pTM': iptm_ptm,
            'min_pDockQ2': min_pdockq_scores,
            'mean_pDockQ2': mean_pdockq_scores,
            'pDockQ2_per_chain': all_pdockq2
        }

        score_dict_with_actifptm = {
            'name': names,
            'pLDDT': plddt_score,
            'pTM': ptm,
            'ipTM': iptm,
            'ipTM+pTM': iptm_ptm,
            'actifpTM': actifptm,
            'actifpTM+pTM': actifptm_ptm,
            'min_pDockQ2': min_pdockq_scores,
            'mean_pDockQ2': mean_pdockq_scores,
            'pDockQ2_per_chain': all_pdockq2,
            'pairwise_actifPTM': pairwise_actifptm,
            'pairwise_ipTM': pairwise_iptm,
            'per_chain_pTM': per_chain_ptm
        }

        selected_score_dict = score_dict_with_actifptm if is_actifptm else score_dict

        collect_data = pd.DataFrame(selected_score_dict)
        sorted_collected_data = collect_data.sort_values(by='name')
        sorted_collected_data.to_csv(score_matrix_filepath, sep='\t', index=False)
        print(f"-- {current_time()} > Completed pDockQ2 score calculation for predicted complex\n")
        return None

