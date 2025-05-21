import pytest
import os
import shutil
import subprocess
from flashfold.utils import (is_installed, get_contents_by_column, file_has_content,
                             load_json_file, get_sum_of_all_file_sizes_by_path, is_valid_protein_a3m,
                             is_valid_protein_fasta, is_valid_path)

# List of dependencies
dependencies = ["datasets", "jackhmmer", "colabfold_batch", "cd-hit"]


@pytest.fixture(scope="session", autouse=True)
def setup_and_remove_output_dir():
    output_dir = "test/output/"
    # Setup: Create the output directory
    os.makedirs(output_dir, exist_ok=True)
    yield output_dir
    # Teardown: Remove the output directory after all tests are finished
    shutil.rmtree(output_dir)


# Checking if dependencies are installed
def test_dependency():
    for dependency in dependencies:
        assert is_installed(dependency), f"{dependency} is not installed. Please install it."


# Checking if 'download_db' subcommand is working
def test_download_db():
    flashfold_sub = "download_db"
    result = subprocess.run(
        ["flashfold", flashfold_sub, "--input", "test/input/others/test_db.json",
         "--output", "test/output/cloud_database/"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr}"
    assert file_has_content("test/output/cloud_database/test-db/sequence_db.fasta"), f"Error: {flashfold_sub}"
    assert file_has_content("test/output/cloud_database/test-db/protein_to_gbks.json"), f"Error: {flashfold_sub}"
    assert file_has_content("test/output/cloud_database/test-db/prot_hash_to_accession.json"), f"Error: {flashfold_sub}"


# Checking if 'ncbi_data' subcommand is working
def test_ncbi_data():
    flashfold_sub = "ncbi_data"
    accession_file_path = "test/input/others/assembly_accessions.txt"
    result = subprocess.run(
        ["flashfold", flashfold_sub, "-i", accession_file_path,
         "-o", "test/output/ncbi_gbff/", "-f", "gbff"],
        capture_output=True,
        text=True
    )
    accessions = get_contents_by_column(accession_file_path, 1)
    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr}"
    for accession in accessions:
        assert file_has_content(f"test/output/ncbi_gbff/gbff/{accession}.gbff"), f"Error: {flashfold_sub}"


# Checking if 'create_db' subcommand is working
def test_create_db():
    flashfold_sub = "create_db"
    result = subprocess.run(
        ["flashfold", flashfold_sub, "-p", "test/output/ncbi_gbff/gbff/", "-o", "test/output/custom_database/"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr}"
    assert file_has_content("test/output/custom_database/sequence_db.fasta"), f"Error: {flashfold_sub}"
    assert file_has_content("test/output/custom_database/protein_to_gbks.json"), f"Error: {flashfold_sub}"
    assert file_has_content("test/output/custom_database/prot_hash_to_accession.json"), f"Error: {flashfold_sub}"


# Checking if 'extend_db' subcommand is working
def test_extend_db():
    flashfold_sub = "extend_db"

    # Extending main database with new database
    result_new_db = subprocess.run(
        ["flashfold", flashfold_sub, "-m", "test/output/custom_database/", "-n",
         "test/output/cloud_database/test-db/", "-y"],
        capture_output=True,
        text=True
    )
    assert result_new_db.returncode == 0, f"{flashfold_sub}-with new database failed with error: {result_new_db.stderr}"
    custom_json = load_json_file("test/output/custom_database/protein_hash_to_accession.json")
    # Extending main database with genbank files
    result_gbff = subprocess.run(
        ["flashfold", flashfold_sub, "-m", "test/output/cloud_database/test-db/", "-g",
         "test/output/ncbi_gbff/gbff/", "-y"],
        capture_output=True,
        text=True
    )
    assert result_gbff.returncode == 0, f"{flashfold_sub}-with gbff failed with error: {result_gbff.stderr}"
    test_json = load_json_file("test/output/cloud_database/test-db/protein_hash_to_accession.json")
    database_ext_with_new_db_size = get_sum_of_all_file_sizes_by_path("test/output/custom_database/")
    database_ext_with_gbff_size = get_sum_of_all_file_sizes_by_path("test/output/cloud_database/test-db/")
    assert database_ext_with_new_db_size > 0, f"Error: {flashfold_sub}"
    assert database_ext_with_gbff_size > 0, f"Error: {flashfold_sub}"
    assert set(custom_json.keys()) == set(test_json.keys()), f"Error: {flashfold_sub}"


# Checking with DNA monomer sequence
def test_monomer_dna():
    flashfold_sub = "fold"
    # testing fold for DNA monomer
    result = subprocess.run(
        ["flashfold", flashfold_sub,
         "-q", "test/input/fasta/monomer_dna.fasta",
         "-d", "test/output/custom_database/",
         "-o", "test/output/fold-monomer-dna/",
         "-t", "1"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr}"
    assert not is_valid_path("test/output/fold-monomer-dna/"), f"Error: {flashfold_sub} with DNA monomer"


# Checking if 'fold' subcommand is working
def test_monomer_msa():
    # testing fold only_msa for monomer
    flashfold_sub = "fold"
    result = subprocess.run(
        ["flashfold", flashfold_sub,
         "-q", "test/input/fasta/monomer.fasta",
         "-d", "test/output/custom_database/",
         "-o", "test/output/fold-monomer/",
         "-t", "1",
         "--only_msa"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr}"

    a3m_file = "test/output/fold-monomer/flashfold_filtered_msa/monomer.a3m"
    assert file_has_content(a3m_file), f"Error: {flashfold_sub} with --only_msa"
    assert is_valid_protein_a3m(a3m_file), f"Error: {flashfold_sub} with --only_msa"
    assert not is_valid_protein_fasta(a3m_file), f"Error: {flashfold_sub} with --only_msa"


def test_multimer_msa():
    flashfold_sub = "fold"
    # testing fold only_msa for heterodimer
    result_msa = subprocess.run(
        ["flashfold", flashfold_sub,
         "-q", "test/input/fasta/heterodimer.fasta",
         "-d", "test/output/custom_database/",
         "-o", "test/output/fold-heterodimer/",
         "-t", "1",
         "--only_msa"],
        capture_output=True,
        text=True
    )
    assert result_msa.returncode == 0, f"{flashfold_sub} failed with error: {result_msa.stderr}"

    a3m_file = "test/output/fold-heterodimer/flashfold_filtered_msa/heterodimer.a3m"
    assert file_has_content(a3m_file), f"Error: {flashfold_sub} with --only_msa"
    assert is_valid_protein_a3m(a3m_file), f"Error: {flashfold_sub} with --only_msa"


def test_json_batch():
    flashfold_sub = "fold"
    # testing fold for batch json output
    result = subprocess.run(
        ["flashfold", flashfold_sub,
         "-q", "test/input/msa/",
         "-o", "test/output/fold-json/",
         "-t", "1",
         "--batch",
         "--only_json"],
        capture_output=True,
        text=True
    )
    monomer_json_file = "test/output/fold-json/monomer/flashfold_json/monomer.json"
    heterodimer_json_file = "test/output/fold-json/heterodimer/flashfold_json/heterodimer.json"
    loaded_monomer = load_json_file(monomer_json_file)
    monomer_paired_msa_value = loaded_monomer['sequences'][0]['protein']['pairedMsa']
    loaded_heterodimer = load_json_file(heterodimer_json_file)
    heterodimer_paired_msa_value_1 = loaded_heterodimer['sequences'][0]['protein']['pairedMsa']
    heterodimer_paired_msa_value_2 = loaded_heterodimer['sequences'][1]['protein']['pairedMsa']
    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr}"
    assert file_has_content(monomer_json_file), f"Error: {flashfold_sub} - batch with json output"
    assert file_has_content(heterodimer_json_file), f"Error: {flashfold_sub} - batch with json output"
    assert monomer_paired_msa_value == "", f"Error: {flashfold_sub} - batch with json output"
    assert heterodimer_paired_msa_value_1 != "", f"Error: {flashfold_sub} - batch with json output"
    assert heterodimer_paired_msa_value_2 != "", f"Error: {flashfold_sub} - batch with json output"


def test_ligand_batch():
    flashfold_sub = "ligand"
    # testing flashfold for ligand batch
    result = subprocess.run(
        ["flashfold", flashfold_sub,
         "-q", "test/input/msa/",
         "-o", "test/output/fold-ligand/",
         "-a", "smiles", "CCOCCC", "1",
         "-a", "ccdCodes", "PRD", "2",
         "-n", "added_smiles_ccdCodes",
         "--batch"],
        capture_output=True,
        text=True
    )
    monomer_json_file = "test/output/fold-ligand/monomer.json"
    heterodimer_json_file = "test/output/fold-ligand/heterodimer.json"
    loaded_monomer = load_json_file(monomer_json_file)
    monomer_ligand_value_1 = loaded_monomer['sequences'][1].get('ligand', None)
    monomer_ligand_value_2 = loaded_monomer['sequences'][2].get('ligand', None)
    monomer_name = loaded_monomer['name']
    loaded_heterodimer = load_json_file(heterodimer_json_file)
    heterodimer_ligand_value_1 = loaded_heterodimer['sequences'][2].get('ligand', None)
    heterodimer_ligand_value_2 = loaded_heterodimer['sequences'][3].get('ligand', None)
    heterodimer_name = loaded_heterodimer['name']

    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr} - batch with ligand"
    assert file_has_content(monomer_json_file), f"Error: {flashfold_sub} - batch with ligand for monomer"
    assert file_has_content(heterodimer_json_file), f"Error: {flashfold_sub} - batch with ligand for heterodimer"
    assert monomer_ligand_value_1, f"Error: {flashfold_sub} - batch with ligand for monomer"
    assert monomer_ligand_value_2, f"Error: {flashfold_sub} - batch with ligand for monomer"
    assert monomer_name.split("-")[-1] == "added_smiles_ccdCodes", \
        f"Error: {flashfold_sub} - batch with ligand for monomer"
    assert heterodimer_ligand_value_1, f"Error: {flashfold_sub} - batch with ligand for heterodimer"
    assert heterodimer_ligand_value_2, f"Error: {flashfold_sub} - batch with ligand for heterodimer"
    assert heterodimer_name.split("-")[-1] == "added_smiles_ccdCodes", \
        f"Error: {flashfold_sub} - batch with ligand for heterodimer"


def test_ligand_remove_ccd_code():
    flashfold_sub = "ligand"
    # testing flashfold for remove ccdCode ligand batch
    result = subprocess.run(
        ["flashfold", flashfold_sub,
         "-q", "test/output/fold-ligand/",
         "-o", "test/output/fold-removed/",
         "-r", "PRD",
         "-n", "removed_prd",
         "--batch"],
        capture_output=True,
        text=True
    )
    monomer_json_file = "test/output/fold-removed/monomer.json"
    heterodimer_json_file = "test/output/fold-removed/heterodimer.json"
    loaded_monomer = load_json_file(monomer_json_file)
    monomer_ligand_value_1 = loaded_monomer['sequences'][1].get('ligand', None)
    monomer_ligands = len(loaded_monomer['sequences']) - 1
    monomer_name = loaded_monomer['name']
    loaded_heterodimer = load_json_file(heterodimer_json_file)
    heterodimer_ligand_value_1 = loaded_heterodimer['sequences'][2].get('ligand', None)
    heterodimer_ligands = len(loaded_heterodimer['sequences']) - 2
    heterodimer_name = loaded_heterodimer['name']

    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr} - remove ccdCode for ligand"
    assert file_has_content(monomer_json_file), f"Error: {flashfold_sub} - remove ccdCode for monomer"
    assert file_has_content(heterodimer_json_file), f"Error: {flashfold_sub} - remove ccdCode for heterodimer"
    assert monomer_ligand_value_1, f"Error: {flashfold_sub} - remove ccdCode for monomer"
    assert monomer_ligands == 1, f"Error: {flashfold_sub} - remove ccdCode for monomer"
    assert monomer_name.split("-")[-1] == "removed_prd", \
        f"Error: {flashfold_sub} - remove ccdCode for monomer"
    assert heterodimer_ligand_value_1, f"Error: {flashfold_sub} - remove ccdCode for heterodimer"
    assert heterodimer_ligands == 1, f"Error: {flashfold_sub} - remove ccdCode for heterodimer"
    assert heterodimer_name.split("-")[-1] == "removed_prd", \
        f"Error: {flashfold_sub} - remove ccdCode for heterodimer"


def test_ligand_purge():
    flashfold_sub = "ligand"
    # testing flashfold for purge ligand batch
    result = subprocess.run(
        ["flashfold", flashfold_sub,
         "-q", "test/output/fold-ligand/",
         "-o", "test/output/fold-purged/",
         "-p",
         "-n", "purged",
         "--batch"],
        capture_output=True,
        text=True
    )
    monomer_json_file = "test/output/fold-purged/monomer.json"
    heterodimer_json_file = "test/output/fold-purged/heterodimer.json"
    loaded_monomer = load_json_file(monomer_json_file)
    monomer_ligands = len(loaded_monomer['sequences']) - 1
    monomer_name = loaded_monomer['name']
    loaded_heterodimer = load_json_file(heterodimer_json_file)
    heterodimer_ligands = len(loaded_heterodimer['sequences']) - 2
    heterodimer_name = loaded_heterodimer['name']

    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr} - purge ligand"
    assert file_has_content(monomer_json_file), f"Error: {flashfold_sub} - purge monomer"
    assert file_has_content(heterodimer_json_file), f"Error: {flashfold_sub} - purge heterodimer"
    assert monomer_ligands == 0, f"Error: {flashfold_sub} -purge monomer"
    assert monomer_name.split("-")[1] == "purged", \
        f"Error: {flashfold_sub} - purge monomer"
    assert heterodimer_ligands == 0, f"Error: {flashfold_sub} - purge heterodimer"
    assert heterodimer_name.split("-")[1] == "purged", \
        f"Error: {flashfold_sub} - purge heterodimer"


def test_fold_batch():
    flashfold_sub = "fold"
    # testing batch only for MSA
    result = subprocess.run(
        ["flashfold", flashfold_sub,
         "-q", "test/input/msa/",
         "--batch",
         "-o", "test/output/fold-batch/",
         "--num_models", "1",
         "--num_recycle", "1"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr}"

    monomer_score_file = "test/output/fold-batch/monomer/flashfold_structure/score.tsv"
    assert file_has_content(monomer_score_file), f"Error: {flashfold_sub}"

    m_score_name_column = get_contents_by_column(monomer_score_file, 1)
    m_score_plddt_column = get_contents_by_column(monomer_score_file, 2)
    m_score_ptm_column = get_contents_by_column(monomer_score_file, 3)
    print(m_score_ptm_column)
    assert len(m_score_name_column) == 2, f"Error: {flashfold_sub}"
    assert float(m_score_plddt_column[1]) > 0, f"Error: {flashfold_sub}"
    assert float(m_score_ptm_column[1]) > 0, f"Error: {flashfold_sub}"

    heterodimer_score_file = "test/output/fold-batch/heterodimer/flashfold_structure/score.tsv"
    assert file_has_content(heterodimer_score_file), f"Error: {flashfold_sub}"

    h_score_name_column = get_contents_by_column(heterodimer_score_file, 1)
    h_score_mean_pdockq2_column = get_contents_by_column(heterodimer_score_file, 2)
    h_score_min_pdockq2_column = get_contents_by_column(heterodimer_score_file, 3)
    print(h_score_min_pdockq2_column)
    assert len(h_score_name_column) == 2, f"Error: {flashfold_sub}"
    assert float(h_score_mean_pdockq2_column[1]) > 0, f"Error: {flashfold_sub}"
    assert float(h_score_min_pdockq2_column[1]) > 0, f"Error: {flashfold_sub}"


def test_summary_report():
    flashfold_sub = "summary"
    result = subprocess.run(
        ["flashfold", flashfold_sub,
         "-d", "test/output/fold-batch/",
         "-o", "test/output/fold-batch/"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0, f"{flashfold_sub} failed with error: {result.stderr}"
    assert file_has_content("test/output/fold-batch/summary.html"), f"Error: {flashfold_sub}"
    assert file_has_content("test/output/fold-batch/summary.csv"), f"Error: {flashfold_sub}"


if __name__ == "__main__":
    pytest.main()
