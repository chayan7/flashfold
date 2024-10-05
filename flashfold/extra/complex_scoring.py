import os
import argparse
from collections import namedtuple
import numpy as np
import pandas as pd
from typing import List, Dict, Union
import ujson
from datetime import datetime

# noinspection PyPackageRequirements
from Bio.PDB.Structure import Structure

# noinspection PyPackageRequirements
from Bio.PDB import PDBParser


def current_time():
    # Get the current time
    present_time = datetime.now()

    # Format the current time in a human-readable way
    return present_time.strftime("%d-%m-%y %H:%M:%S")


# Function to load stored history info from file
def load_json_file(file_path: str) -> Dict:
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            return ujson.load(file)
    else:
        return {}


def get_filename_without_extension(file_path: str) -> str:
    """
    Get the file name without its extension from a given file path.

    Parameters:
    file_path (str): The full path to the file.

    Returns:
    str: The file name without its extension.
    """
    base_name = os.path.basename(file_path)  # Extracts the file name from the path
    file_name_without_ext = os.path.splitext(base_name)[0]  # Removes the extension
    return file_name_without_ext


PdockQ_Calc_Out = namedtuple('PdockQ_Calc_Out', ["num_res_list", "if_plddt_list", "plddt_list",
                                                 "num_diso_list", "pdockq_list", "length_list"])


def load_structure(pdb_file_path: str) -> Structure:
    """
    Loads a PDB structure from a file and returns a Structure object.

    :param pdb_file_path: Path to the PDB file.
    :return: A Structure object representing the PDB structure.
    """
    # Extract the structure ID from the file name (remove the '.pdb' extension)
    structure_id = get_filename_without_extension(pdb_file_path)
    # Initialize the PDB parser
    bio_parser = PDBParser()
    return bio_parser.get_structure(structure_id, pdb_file_path)


def get_popt_from_cutoff(cutoff: Union[int, float]) -> List[float]:
    # These popt values are collected from https://gitlab.com/ElofssonLab/huintaf2/-/raw/main/bin/pDockQ2.py.
    if cutoff <= 5:
        # 5 SpearmanrResult(correlation=0.7647585237390458, pvalue=4.749030057232305e-280)
        return [6.96234405e-01, 2.35483775e+02, 2.25322970e-02, 2.88445245e-02]
        # 0.7805034405869632
    elif cutoff <= 6:
        # 6 SpearmanrResult(correlation=0.7708834427476546, pvalue=2.707297682746201e-287)
        return [7.02605033e-01, 2.91749822e+02, 2.70621128e-02, 2.25416051e-02]
        # 0.7871982094514278
    elif cutoff <= 7:
        # 7 SpearmanrResult(correlation=0.7709518988131879, pvalue=2.2402500804327052e-287)
        return [7.06385097e-01, 3.32456259e+02, 2.97005237e-02, 2.24488132e-02]
        # 0.7859609807320201
    elif cutoff <= 8:
        # 8 SpearmanrResult(correlation=0.7632969367380509, pvalue=2.3583905451705336e-278)
        return [7.18442739e-01, 3.60791204e+02, 3.01635944e-02, 2.04076969e-02]
        # 0.7764648775754815
    elif cutoff <= 9:
        # 9 SpearmanrResult(correlation=0.7496303495195178, pvalue=4.539049646719674e-263)
        return [7.23328534e-01, 3.80036094e+02, 3.06316084e-02, 1.98471192e-02]
        # 0.7608417399783565
    elif cutoff <= 10:
        # 10 SpearmanrResult(correlation=0.7330653937901442, pvalue=7.988440779428826e-246)
        return [7.20293782e-01, 3.95627723e+02, 3.15235037e-02, 2.37304238e-02]
        # 0.7431426093979494
    elif cutoff <= 11:
        # 11 SpearmanrResult(correlation=0.71288058226417, pvalue=1.7542846392453894e-226)
        return [7.22015998e-01, 4.09095024e+02, 3.11905555e-02, 2.59467513e-02]
        # 0.7219615906164123
    else:
        # 12 SpearmanrResult(correlation=0.6938911161134763, pvalue=9.284495013784153e-210)
        return [7.20555781e-01, 4.21033584e+02, 3.09024241e-02, 2.88659629e-02]
        # 0.7023000652310362


# This function is collected from https://gitlab.com/ElofssonLab/huintaf2/-/raw/main/bin/pDockQ2.py and modified.

def calc_pdockq(structure_bio: Structure, cutoff_score: Union[int, float], diso_cut_1=50, diso_cut_2=70,
                diso_cut_3=90) -> PdockQ_Calc_Out:
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
                            elif dist > 2 * cutoff_score:  # To speed up things
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
            if res['CA'].get_bfactor() > diso_cut_3:  # >90
                num_diso1[0] += 1
            elif res['CA'].get_bfactor() > diso_cut_2:  # 70-90
                num_diso1[1] += 1
            elif res['CA'].get_bfactor() > diso_cut_1:  # 50-70
                num_diso1[2] += 1
            else:  # <50
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


# This function is collected from https://gitlab.com/ElofssonLab/huintaf2/-/raw/main/bin/pDockQ2.py and modified.
def sigmoid(x: float, l_cap: float, x0: float, k: float, b: float) -> float:
    y = l_cap / (1 + np.exp(-k * (x - x0))) + b
    return y


def generate_score_matrix(path_to_predicted_structure: str, cutoff: Union[int, float]) -> None:

    # List all files in the directory
    all_files = os.listdir(path_to_predicted_structure)

    # Filter for only .pdb files
    pdb_files = [f for f in all_files if f.endswith(".pdb")]

    names = []
    pdockq_scores_per_chain = []
    pdockq_scores_with_chain_name = []
    iptm = []
    iptm_ptm = []
    plddt_score = []
    for predicted_pdb in pdb_files:

        # Extracting information from the pdb_file_path
        pdb_file = predicted_pdb.split("/")[-1]

        # Check if the file name contains 'rank'
        if 'rank' not in pdb_file:
            print(f"Skipping {pdb_file} as it does not contain 'rank' in the file name.")
            return

        json_file_name = predicted_pdb.replace(".pdb", ".json").replace("unrelaxed",
                                                                             "scores")
        json_file_path = os.path.join(path_to_predicted_structure, json_file_name)
        names.append(predicted_pdb)
        json_file = load_json_file(json_file_path)

        # Obtain ipTM
        iptm_score = json_file['iptm']
        iptm.append(iptm_score)
        iptm_ptm_score = 0.8 * json_file["iptm"] + 0.2 * json_file["ptm"]
        iptm_ptm.append(iptm_ptm_score)

        # Mean plDDT
        mean_plddt = np.mean(json_file["plddt"])
        plddt_score.append(mean_plddt)

        # Calculate pdockq2 score
        pdb_file_path = os.path.join(path_to_predicted_structure, predicted_pdb)
        structure = load_structure(pdb_file_path)

        num_res, if_plddt, plddt, num_diso, pdockq, length = calc_pdockq(structure, cutoff)
        mean_pdockq = sum(pdockq) / len(pdockq)
        pdockq_scores_per_chain.append(mean_pdockq)

        chain_counter = 0
        pdockq_scores_per_model = []
        for chain_unit in structure.get_chains():
            chain_id = chain_unit.get_id()
            pdockq_score_with_chain = f"{chain_id}: {pdockq[chain_counter]}"
            pdockq_scores_per_model.append(pdockq_score_with_chain)
            chain_counter += 1

        pdockq_scores_with_chain_name.append(pdockq_scores_per_model)

    score_dict = {'name': names, 'pDockQ': pdockq_scores_per_chain, 'ipTM': iptm, 'ipTM+pTM': iptm_ptm,
                  'mplDDT': plddt_score, 'pDockQ_per_chain': pdockq_scores_with_chain_name}
    collect_data = pd.DataFrame(score_dict)
    sorted_collected_data = collect_data.sort_values(by='name')
    score_matrix_filepath_for_complex = os.path.join(path_to_predicted_structure, "score_new.tsv")
    sorted_collected_data.to_csv(score_matrix_filepath_for_complex, sep='\t', index=False)
    print(f"-- {current_time()} > Completed pDockQ score calculation for predicted complex\n")


def main():
    parser = argparse.ArgumentParser(description="Calculates pDockQ from NumRes and IF_plDDT")
    parser.add_argument("-d", "--directory", required=True,
                        help="Path to the predicted structure directory.")
    parser.add_argument("-c", "--cutoff", type=float, required=False, default=10.0,
                        help="Cutoff for defining distances")
    args = parser.parse_args()

    generate_score_matrix(args.directory, args.cutoff)


if __name__ == "__main__":
    main()

