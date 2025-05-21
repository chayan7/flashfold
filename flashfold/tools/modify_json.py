# This code is collected from: https://github.com/cddlab/alphafold3_tools/blob/main/alphafold3tools/modjson.py
# And changed accordingly to fit the project requirements.


import copy
from typing import Literal, cast, List, Dict, Optional
from flashfold.tools.util import int_id_to_str_id
from flashfold.utils import load_json_file, write_dict_to_json_as_file, current_time


def remove_ccdcodes_from_data(data: Dict, ccdcodes_to_remove: list[str]) -> Dict:
    """Removes ligand entities from AlphaFold3 json data.
    Args:
        data (Dict): AlphaFold3 json data.
        ccdcodes_to_remove (list[str]): ccdcodes to remove.
    Returns:
        dict (Dict): AlphaFold3 json data with ccdcodes removed.
    """
    new_data = copy.deepcopy(data)
    sequence_contents = new_data["sequences"]
    new_sequence_contents = []

    is_removed = False
    for sequence_content in sequence_contents:
        if "ligand" in sequence_content:
            if "ccdCodes" in sequence_content["ligand"]:
                ccd_codes = sequence_content["ligand"]["ccdCodes"]
                if any(ligand in ccd_codes for ligand in ccdcodes_to_remove):
                    print(
                        f"\t ∞ Removing ligand: {sequence_content['ligand']['ccdCodes']}"
                    )
                    is_removed = True
                else:
                    new_sequence_contents.append(sequence_content)
            else:
                # Keeping other ligands (e.g., smiles) if present.
                new_sequence_contents.append(sequence_content)
        else:
            new_sequence_contents.append(sequence_content)

    if not is_removed:
        print(
            f"--Warning: No ligand with ccdCodes: {ccdcodes_to_remove} found "
            "in the input JSON file."
        )
    new_data["sequences"] = new_sequence_contents
    return new_data


def purge_all_ligands(data: Dict) -> Dict:
    """Purges all ligand entities from AlphaFold3 json data.
    Args:
        data (Dict): AlphaFold3 json data.
    Returns:
        dict (Dict): AlphaFold3 json data with ligands removed.
    """
    new_data = copy.deepcopy(data)
    sequence_contents = new_data["sequences"]
    new_sequence_contents = []

    for sequence_content in sequence_contents:
        if "ligand" in sequence_content:
            if sequence_content["ligand"].get("smiles"):
                print(f"\t ∞ Purging smiles: {sequence_content['ligand']['smiles']}")
            elif sequence_content["ligand"].get("ccdCodes"):
                print(f"\t ∞ Purging ccdCodes: {sequence_content['ligand']['ccdCodes']}")
            else:
                new_sequence_contents.append(sequence_content)
        else:
            new_sequence_contents.append(sequence_content)
    new_data["sequences"] = new_sequence_contents

    data_name = data["name"]
    new_data_name = data_name.split("-")[0]
    if new_data_name != data_name:
        print(f"\t ∞ Resetting the job name to: {new_data_name} (previous: {data_name})")
        new_data["name"] = new_data_name
    return new_data


def add_ligand_to_data(
    data: Dict,
    ligand_type: Literal["smiles", "ccdCodes"],
    ligand_name: str,
    num_ligand: int,
) -> Dict:
    """Adds ligand entities to AlphaFold3 json data.
    Args:
        data (dict): AlphaFold3 json data.
        ligand_type (Literal["smiles", "ccdCodes"]): Type of ligand to add.
        ligand_name (str): Ligand string to add.
        num_ligand (int): Number of ligand molecules to add.
    Returns:
        dict (Dict): AlphaFold3 json data with ligands added.
    """
    print(f"\t ∞ Adding {num_ligand} ligand: {ligand_name} as {ligand_type}")
    new_data = copy.deepcopy(data)
    sequence_contents = new_data["sequences"]

    num_ids = [int_id_to_str_id(num) for num in range(1, num_ligand + 1)]
    if ligand_type == "smiles":
        sequence_contents.append(
            {
                "ligand": {
                    "id": num_ids,
                    "smiles": ligand_name,
                }
            }
        )
    elif ligand_type == "ccdCodes":
        sequence_contents.append(
            {
                "ligand": {
                    "id": num_ids,
                    "ccdCodes": [ligand_name],
                }
            }
        )
    return new_data


def fix_sequence_ids(data: Dict) -> Dict:
    """Fixes the sequence IDs in the AlphaFold3 JSON data.

    This function updates the IDs in the "sequences" field of the provided
    AlphaFold3 JSON data. It ensures that each ID is unique and follows a
    sequential order using the `int_id_to_str_id` function to convert integers
    to string IDs in a reverse spreadsheet style naming
    (e.g., 1 = A, 2 = B, ..., 27 = AA).

    Args:
        data (Dict): The AlphaFold3 JSON data containing sequences with IDs to be fixed.

    Returns:
        dict (dict): A new dictionary with the updated sequence IDs.
    """
    new_data = copy.deepcopy(data)
    sequence_contents = new_data["sequences"]

    id_counter = 1  # 1-based indexing.
    for sequence_content in sequence_contents:
        for key in sequence_content:
            if "id" in sequence_content[key]:
                if isinstance(sequence_content[key]["id"], list):
                    new_ids = []
                    for _ in sequence_content[key]["id"]:
                        new_ids.append(int_id_to_str_id(id_counter))
                        id_counter += 1
                    sequence_content[key]["id"] = new_ids
                elif isinstance(sequence_content[key]["id"], str):
                    new_id = int_id_to_str_id(id_counter)
                    sequence_content[key]["id"] = new_id
                    id_counter += 1

    return new_data


def modify_name(data: Dict, new_name: str, is_batch: bool) -> Dict:
    """Modifies the job name in the AlphaFold3 JSON data.

    Args:
        data (dict): The AlphaFold3 JSON data.
        new_name (str): The new job name to set.
        is_batch (bool): A boolean flag to determine if the job is a batch job.
    Returns:
        dict (Dict): A new dictionary with the updated prediction name.
    """
    new_data = copy.deepcopy(data)

    modified_name = new_name
    if is_batch:
        # If the job is a batch job, append the new name to the previous name.
        previous_name = new_data.get("name", None)
        if previous_name:
            modified_name = f"{previous_name}-{new_name}"

    print(f"\t ∞ Setting the job name to: {modified_name}")
    new_data["name"] = modified_name
    return new_data


def add_userccd(data: Dict, userccd_files: List[str]) -> Dict:
    """Adds user provided ccdCodes to the AlphaFold3 JSON data.

    Args:
        data (Dict): The AlphaFold3 JSON data.
        userccd_files (List[str]): The path to the user provided ccdCodes file.
        Multiple files can be provided.
    Returns:
        dict (Dict): A new dictionary with the updated ccdCodes.
    """
    new_data = copy.deepcopy(data)
    userccd_as_string = ""
    for userccd_file in userccd_files:
        with open(userccd_file, "r") as file:
            userccd_as_string += file.read()
            userccd_as_string += "## \n"

    new_data["userCCD"] = userccd_as_string
    return new_data


def modjson(
    input_path: str,
    output_path: str,
    add_ligand: Optional[List[List[str]]],
    purge_ligands: bool,
    remove_ccdcodes: Optional[List[str]],
    name: Optional[str],
    userccdfiles: Optional[List[str]],
    is_batch: bool,
) -> None:
    """Modifies AlphaFold3 JSON file.
    Args:
        input_path (str): Input AlphaFold3 JSON file.
        output_path (str): Output JSON file.
        add_ligand (Optional[List[List[str]]]): Add ligand to the input JSON file.
        purge_ligands (bool): Purge all ligands from the input JSON file at first.
        remove_ccdcodes (Optional[List[str]]): Remove ligands with ccdcodes from the input JSON file.
        name (Optional[str]): Set the job name in the input JSON file.
        userccdfiles (Optional[List[str]]): Add user provided ccdCodes to the input JSON file.
        is_batch (bool): A boolean flag to determine if the job is a batch job.
    """
    print(f'\n-- {current_time()} > JSON file modification is being processed')
    data = load_json_file(input_path)
    if purge_ligands:
        print("\t ∞ Purging current ligand entities from the input JSON file.")
        data = purge_all_ligands(data)
    if remove_ccdcodes:
        data = remove_ccdcodes_from_data(data, remove_ccdcodes)
    if add_ligand:
        for ligand_info in add_ligand:
            ligand_type, ligand_name, num_ligand = ligand_info
            if ligand_type not in ["smiles", "ccdCodes"]:
                raise ValueError(
                    f"Invalid ligand type: {ligand_type}. "
                    "The ligand type must be either 'smiles' or 'ccdCodes'."
                )
            ligand_type_literal = cast(Literal["smiles", "ccdCodes"], ligand_type)
            data = add_ligand_to_data(data, ligand_type_literal, ligand_name, int(num_ligand))
    data = fix_sequence_ids(data)

    if name:
        data = modify_name(data, name, is_batch)

    if userccdfiles:
        data = add_userccd(data, userccdfiles)

    write_dict_to_json_as_file(data, output_path)
    print(f"\t ∞ Output file: '{output_path}'")

    print(f'-- {current_time()} > JSON file modification is being processed \n')

    return None

