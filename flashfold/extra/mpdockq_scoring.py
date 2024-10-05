import os
import argparse
import sys
from collections import defaultdict
from collections import namedtuple
import numpy as np
import pandas as pd
from typing import Dict, Union
import ujson
from datetime import datetime
import math

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


Chain_info = namedtuple('Chain_info', ['coordinates', 'plddt'])


def parse_atm_record(line: str) -> Dict:
    # Get the atm record
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record


def read_pdb(pdb_file: str) -> Chain_info:
    # Read a pdb file predicted with AF and rewritten to contain all chains
    chain_coordinates = {}
    chain_plddt = {}
    with open(pdb_file, 'r') as file:
        for line in file:
            if not line.startswith('ATOM'):
                continue
            record = parse_atm_record(line)
            # Get CB - CA for GLY
            if record['atm_name'] == 'CB' or (record['atm_name'] == 'CA' and record['res_name'] == 'GLY'):
                if record['chain'] in [*chain_coordinates.keys()]:
                    chain_coordinates[record['chain']].append([record['x'], record['y'], record['z']])
                    chain_plddt[record['chain']].append(record['B'])
                else:
                    chain_coordinates[record['chain']] = [[record['x'], record['y'], record['z']]]
                    chain_plddt[record['chain']] = [record['B']]

    # Convert to arrays
    for chain in chain_coordinates:
        chain_coordinates[chain] = np.array(chain_coordinates[chain])
        chain_plddt[chain] = np.array(chain_plddt[chain])

    return Chain_info(chain_coordinates, chain_plddt)


# alphapulldown

def read_pdb_file(pdb_file: str):
    # Read a pdb file per chain
    pdb_chains = {}
    chain_coords = {}
    chain_ca_inds = {}
    chain_cb_inds = {}

    with open(pdb_file) as file:
        for line in file:
            if 'ATOM' in line:
                record = parse_atm_record(line)
                if record['chain'] in [*pdb_chains.keys()]:
                    pdb_chains[record['chain']].append(line)
                    chain_coords[record['chain']].append([record['x'], record['y'], record['z']])
                    coord_ind += 1
                    if record['atm_name'] == 'CA':
                        chain_ca_inds[record['chain']].append(coord_ind)
                    if record['atm_name'] == 'CB' or (record['atm_name'] == 'CA' and record['res_name'] == 'GLY'):
                        chain_cb_inds[record['chain']].append(coord_ind)

                else:
                    pdb_chains[record['chain']] = [line]
                    chain_coords[record['chain']] = [[record['x'], record['y'], record['z']]]
                    chain_ca_inds[record['chain']] = []
                    chain_cb_inds[record['chain']] = []
                    # Reset coord ind
                    coord_ind = 0

    return pdb_chains, chain_coords, chain_ca_inds, chain_cb_inds


def parse_bfactor(pdb_file: str) -> np.array:
    parser = PDBParser()
    structure = parser.get_structure('PDB_structure', pdb_file)
    residue_bfactors = []
    for model in structure:
        for chain in model:
            for residue in chain:
                bfactor_sum = 0
                atom_count = 0
                for atom in residue:
                    bfactor_sum += atom.get_bfactor()
                    atom_count += 1
                if atom_count > 0:
                    avg_bfactor = bfactor_sum / atom_count
                    residue_bfactors.append(avg_bfactor)
    return np.array(residue_bfactors)


def get_best_plddt(work_dir: str):
    best_plddt = None

    # List all files in the directory
    all_files = os.listdir(work_dir)

    # Filter for only .pdb files
    pdb_files = [f for f in all_files if f.endswith(".pdb")]

    best_ranked_pdb_file = None
    for predicted_pdb in pdb_files:
        # Extracting information from the pdb_file_path
        pdb_file = predicted_pdb.split("/")[-1]

        # Check if the file name contains 'rank'
        if 'rank' in pdb_file:
            best_ranked_pdb_file = pdb_file
    try:
        if os.path.exists(os.path.join(work_dir, best_ranked_pdb_file)):
            best_plddt = parse_bfactor(os.path.join(work_dir, best_ranked_pdb_file))
            print(f"Successfully parsed plddt values")
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        print(
            f"rank not found in {work_dir}. Failed to parse information of pLDDT scores of the best"
            f" model. The programme will crash.")

    return best_plddt


def read_plddt(best_plddt, chain_ca_inds):
    # Get the plDDT for each chain
    chain_names = chain_ca_inds.keys()
    chain_lengths = dict()
    for name in chain_names:
        curr_len = len(chain_ca_inds[name])
        chain_lengths[name] = curr_len

    plddt_per_chain = dict()
    curr_len = 0
    for k, v in chain_lengths.items():
        curr_plddt = best_plddt[curr_len:curr_len + v]
        plddt_per_chain[k] = curr_plddt
        curr_len += v
    return plddt_per_chain


def score_complex(path_coords, path_cb_inds, path_plddt):
    """

    Score all interfaces in the current complex

    Modified from the score_complex() function in MoLPC repo:
    https://gitlab.com/patrickbryant1/molpc/-/blob/main/src/complex_assembly/score_entire_complex.py#L106-154

    """

    chains = [*path_coords.keys()]
    chain_inds = np.arange(len(chains))
    complex_score = 0
    # Get interfaces per chain
    for i in chain_inds:
        chain_i = chains[i]
        chain_coords = np.array(path_coords[chain_i])
        chain_cb_inds = path_cb_inds[chain_i]
        l1 = len(chain_cb_inds)
        chain_cb_coords = chain_coords[chain_cb_inds]
        chain_plddt = path_plddt[chain_i]

        for int_i in np.setdiff1d(chain_inds, i):
            int_chain = chains[int_i]
            int_chain_cb_coords = np.array(path_coords[int_chain])[path_cb_inds[int_chain]]
            int_chain_plddt = path_plddt[int_chain]
            # Calc 2-norm
            mat = np.append(chain_cb_coords, int_chain_cb_coords, axis=0)
            a_min_b = mat[:, np.newaxis, :] - mat[np.newaxis, :, :]
            dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
            contact_dists = dists[:l1, l1:]
            contacts = np.argwhere(contact_dists <= 8)
            # The first axis contains the contacts from chain 1
            # The second the contacts from chain 2
            if contacts.shape[0] > 0:
                av_if_plddt = np.concatenate((chain_plddt[contacts[:, 0]], int_chain_plddt[contacts[:, 1]])).mean()
                complex_score += np.log10(contacts.shape[0] + 1) * av_if_plddt

    return complex_score, len(chains)


def calculate_mpdockq(complex_score):
    """
    A function that returns a complex's mpDockQ score after
    calculating complex_score
    """
    cap_l = 0.827
    x_0 = 261.398
    k = 0.036
    b = 0.221
    return cap_l / (1 + math.exp(-1 * k * (complex_score - x_0))) + b

# alphapulldown #


def calc_pdockq(chain_coordinates: Dict, chain_plddt: Dict, t: int) -> Union[int, float]:
    """
    Calculate the pDockQ scores
    pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
    L= 0.724 x0= 152.611 k= 0.052 and b= 0.018
    Args:
        chain_coordinates: Dict
        chain_plddt: Dict
        t: int

    Returns:
        pdockq: float
    """

    # Get coordinates and plddt per chain

    # noinspection PyTupleAssignmentBalance
    ch1, ch2 = [*chain_coordinates.keys()]
    coords1, coords2 = chain_coordinates[ch1], chain_coordinates[ch2]
    plddt1, plddt2 = chain_plddt[ch1], chain_plddt[ch2]

    # Calc 2-norm
    mat = np.append(coords1, coords2, axis=0)
    a_min_b = mat[:, np.newaxis, :] - mat[np.newaxis, :, :]
    dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
    l1 = len(coords1)
    contact_dists = dists[:l1, l1:]  # upper triangular --> first dim = chain 1
    contacts = np.argwhere(contact_dists <= t)

    if contacts.shape[0] < 1:
        pdockq = 0
    else:
        # Get the average interface plDDT
        avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:, 0])],
                                                  plddt2[np.unique(contacts[:, 1])]]))
        # Get the number of interface contacts
        n_if_contacts = contacts.shape[0]
        x = avg_if_plddt*np.log10(n_if_contacts)
        pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018

    return pdockq


def generate_score_matrix(path_to_predicted_structure: str, is_monomer: bool) -> None:

    if is_monomer:
        return None

    # List all files in the directory
    all_files = os.listdir(path_to_predicted_structure)

    # Filter for only .pdb files
    pdb_files = [f for f in all_files if f.endswith(".pdb")]

    names = []
    pdockq = []
    mpdockq = []
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

        # Read chains
        pdb_file_path = os.path.join(path_to_predicted_structure, predicted_pdb)
        pdb = read_pdb(pdb_file_path)
        print(pdb_file_path)
        chain_coordinates = pdb.coordinates
        chain_plddt = pdb.plddt

        pdb_chains, chain_coords, chain_ca_inds, chain_cb_inds = read_pdb_file(pdb_file_path)
        best_plddt = parse_bfactor(pdb_file_path)
        plddt_per_chain = read_plddt(best_plddt, chain_ca_inds)
        complex_score, len_chains = score_complex(chain_coords, chain_cb_inds, plddt_per_chain)
        mpdockq_score = calculate_mpdockq(complex_score)
        mpdockq.append(mpdockq_score)

        # Check chains
        chain_number = len(chain_coordinates.keys())
        if chain_number < 2:
            print(f'Warning: pDockQ score is not calculated because predicted pdb file contains {chain_coordinates} '
                  f'chain. FoldFlash currently offers calculation when the pdb file contains at least two chains.\n')
            sys.exit()

        # Calculate pdockq
        t = 8   # Distance threshold, set to 8 Å
        compute_pdockq = calc_pdockq(chain_coordinates, chain_plddt, t)
        pdockq.append(compute_pdockq)

    score_dict = {'name': names, 'mpDockQ': mpdockq, 'pDockQ': pdockq,
                  'ipTM': iptm, 'ipTM+pTM': iptm_ptm, 'mplDDT': plddt_score}
    collect_data = pd.DataFrame(score_dict)
    sorted_collected_data = collect_data.sort_values(by='name')
    score_matrix_filepath_for_complex = os.path.join(path_to_predicted_structure, "score.tsv")
    sorted_collected_data.to_csv(score_matrix_filepath_for_complex, sep='\t', index=False)
    print(f"-- {current_time()} > Completed pDockQ score calculation for predicted complex\n")


def main():
    parser = argparse.ArgumentParser(description="Calculate mpDockQ score for whole complex.")
    parser.add_argument("-d", "--directory", required=True, help="Path to the predicted structure directory.")
    args = parser.parse_args()

    generate_score_matrix(args.directory, False)


if __name__ == "__main__":
    main()
