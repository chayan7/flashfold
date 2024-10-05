# Following script is collected from https://gitlab.com/ElofssonLab/huintaf2/-/raw/main/bin/pDockQ2.py and modified.

import numpy as np
import argparse
import re
from typing import Union
from collections import namedtuple

# noinspection PyPackageRequirements
from Bio.PDB.PDBParser import PDBParser
# noinspection PyPackageRequirements
from Bio.PDB.Structure import Structure
# noinspection PyPackageRequirements
from Bio.PDB.MMCIFParser import MMCIFParser


PdockQ_Calc_Out = namedtuple('PdockQ_Calc_Out', ["num_res_list", "if_plddt_list", "plddt_list",
                                                 "num_diso_list", "pdockq_list", "length_list"])


def calc_pdockq(structure_bio: Structure, cutoff_score: Union[int, float], diso_cut_1=50, diso_cut_2=70,
                diso_cut_3=90) -> PdockQ_Calc_Out:
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


def sigmoid(x: float, l_cap: float, x0: float, k: float, b: float) -> float:
    y = l_cap / (1 + np.exp(-k * (x - x0))) + b
    return y


# popt=[7.07140240e-01, 3.88062162e+02, 3.14767156e-02, 3.13182907e-02]

arg_parser = argparse.ArgumentParser(description="Calculates pDockQ from NumRes and IF_plDDT")
group = arg_parser.add_mutually_exclusive_group(required=True)
group.add_argument("-p", "--pdb", type=argparse.FileType('r'), help="Input pdb file pLddt values in "
                                                                   "b_factor columns")
group.add_argument("-c", "--cif", type=argparse.FileType('r'), help="Input cif file plddt values in "
                                                                   "b_factor columns")
arg_parser.add_argument("-v", "--verbose", action='store_true', required=False, help="Verbose output")
arg_parser.add_argument("-C", "--cutoff", type=float, required=False, default=10.0,
                        help="Cutoff for defining distances")
args = arg_parser.parse_args()

cutoff = args.cutoff
cutoff2 = 3

structure = None
file_name = None
if args.cif:
    file_name = args.cif.name
    bio_parser = MMCIFParser()
    structure_file = args.cif
    structure_id = args.cif.name[:-4]
    structure = bio_parser.get_structure(structure_id, structure_file)
elif args.pdb:
    file_name = args.pdb.name
    bio_parser = PDBParser()
    structure_file = args.pdb
    structure_id = args.pdb.name[:-4]
    structure = bio_parser.get_structure(structure_id, structure_file)


if cutoff <= 5:
    # 5 SpearmanrResult(correlation=0.7647585237390458, pvalue=4.749030057232305e-280)
    popt = [6.96234405e-01, 2.35483775e+02, 2.25322970e-02, 2.88445245e-02]
    # 0.7805034405869632
elif cutoff <= 6:
    # 6 SpearmanrResult(correlation=0.7708834427476546, pvalue=2.707297682746201e-287)
    popt = [7.02605033e-01, 2.91749822e+02, 2.70621128e-02, 2.25416051e-02]
    # 0.7871982094514278
elif cutoff <= 7:
    # 7 SpearmanrResult(correlation=0.7709518988131879, pvalue=2.2402500804327052e-287)
    popt = [7.06385097e-01, 3.32456259e+02, 2.97005237e-02, 2.24488132e-02]
    # 0.7859609807320201
elif cutoff <= 8:
    # 8 SpearmanrResult(correlation=0.7632969367380509, pvalue=2.3583905451705336e-278)
    popt = [7.18442739e-01, 3.60791204e+02, 3.01635944e-02, 2.04076969e-02]
    # 0.7764648775754815
elif cutoff <= 9:
    # 9 SpearmanrResult(correlation=0.7496303495195178, pvalue=4.539049646719674e-263)
    popt = [7.23328534e-01, 3.80036094e+02, 3.06316084e-02, 1.98471192e-02]
    # 0.7608417399783565
elif cutoff <= 10:
    # 10 SpearmanrResult(correlation=0.7330653937901442, pvalue=7.988440779428826e-246)
    popt = [7.20293782e-01, 3.95627723e+02, 3.15235037e-02, 2.37304238e-02]
    # 0.7431426093979494
elif cutoff <= 11:
    # 11 SpearmanrResult(correlation=0.71288058226417, pvalue=1.7542846392453894e-226)
    popt = [7.22015998e-01, 4.09095024e+02, 3.11905555e-02, 2.59467513e-02]
    # 0.7219615906164123
else:
    # 12 SpearmanrResult(correlation=0.6938911161134763, pvalue=9.284495013784153e-210)
    popt = [7.20555781e-01, 4.21033584e+02, 3.09024241e-02, 2.88659629e-02]
    # 0.7023000652310362

num_res, if_plddt, plddt, num_diso, pdockq, length = calc_pdockq(structure, cutoff)
num_res_overlap, if_plddt_overlap, plddt_overlap, num_diso_overlap, pdockq_overlap, length_overlap = calc_pdockq(
    structure, cutoff2)

file_name = re.sub(r'.*/', '', file_name)
file_name = re.sub(r'.pdb$', '', file_name)
file_name = re.sub(r'.cif$', '', file_name)

# print (NumRes,tiny,IF_plDDT, popt)

# print (NumRes,IF_plDDT,plDDT,NumDiso,pdockq,NumResOverlap)

if args.verbose:
    chain_counter = 0
    print("Name,id,Chain,pDockQ,NumRes,IF_plDDT,plDDTNumDiso1+90,NumDiso1-70-90,NumDiso1-50-70,NumDiso1-50,Length")
    for chain_unit in structure.get_chains():
        chain_id = chain_unit.get_id()
        print("%s,%d,%s,%f,%d,%f,%f,%d,%d,%d,%d,%d,%d" % (
            file_name, chain_counter, chain_id, pdockq[chain_counter], num_res[chain_counter], if_plddt[chain_counter],
            plddt[chain_counter], num_diso[chain_counter][0], num_diso[chain_counter][1], num_diso[chain_counter][2],
            num_diso[chain_counter][3], num_res_overlap[chain_counter], length[chain_counter]))
        chain_counter += 1
else:
    chain_counter = 0
    for chain_unit in structure.get_chains():
        chain_id = chain_unit.get_id()
        print(chain_counter, chain_id, pdockq[chain_counter])
        chain_counter += 1
