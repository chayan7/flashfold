import argparse
from flashfold.scripts import create_protein_db_from_gbk, download_database_from_cloud, predict_3d_structure, \
    extend_main_sequence_db, download_ncbi_data, parse_formats, make_summary_report, make_json_with_ligand, \
    process_alphafold3_task, predict_stoichiometry
from flashfold.utils import is_zero_or_pos_int, is_pos_int
from pathlib import Path


def main() -> None:
    # Create the main parser
    parser = argparse.ArgumentParser(description='Predict protein structure from sequence.')
    parser.add_argument("-v", "--version", action="version", version="∞∞ FlashFold v1.2.0 ∞∞")

    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', title='subcommands',
                                       description='Choose any of the following options:')

    # command create_db parser
    desc_create_db = ''' Create sequence database from Genbank files, by extracting protein sequences. 
    Files should have either '.gbff' or '.gbk' extension. '''
    create_db = subparsers.add_parser('create_db', description=desc_create_db)
    create_db.add_argument("-p", "--path", metavar="<Data_Dir>", required=True,
                           help="path to the directory that contains all '.gbff' or '.gbk' files")
    create_db.add_argument("-o", "--output", metavar="<Output_Dir>", required=True,
                           help="output directory to store sequence database")

    # command extend_db parser
    desc_extend_db = ''' Extend current sequence database with information from new genbank file(s) or another database 
    that has been created using create_db command. '''
    extend_db = subparsers.add_parser('extend_db', description=desc_extend_db)
    extend_db.add_argument("-m", "--main_db", metavar="<Database_Dir>", required=True,
                           help="path to the main database, which will be extended with new information")
    extend_db.add_argument("-y", "--yes", action="store_true", default=False,
                           help="confirm the extension of the main database (default: False)")
    # Create a mutually exclusive group for -n and -g
    group = extend_db.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--genbank_path", metavar="<Genbank_Dir>",
                       help="path to the directory that contains genbank file(s), will be used to "
                            "update the main database")
    group.add_argument("-n", "--new_db", metavar="<Database_Dir>",
                       help="path to the new database that will be added to the main database")
    # command download_db parser
    desc_download_db = ''' Download database required for FlashFold prediction. '''
    download_db = subparsers.add_parser('download_db', description=desc_download_db)
    download_db.add_argument("-i", "--input", metavar="<FILE_In>", required=True,
                             help="path to the JSON (database.json) file that contains the download links")
    download_db.add_argument("-o", "--output", metavar="<Output_Dir>", required=True,
                             help="path to store the database built from representative sequence")

    # command ncbi_data parser
    desc_ncbi_data = '''
        Description: 

            This aims to retrieve genome assembly data from NCBI. 
            Users can download the data in the following formats:
            (1) genome - Genomic sequences
            (2) protein - Amino acid sequences
            (3) cds - Nucleotide coding sequences
            (4) gff3 - General feature file (GFF3)
            (5) gtf - Gene transfer format (GTF)
            (6) gbff - GenBank flat file (GBFF)
            (7) all - All the above formats.
            For FlashFold, only the GenBank or 'gbff' format is required.
        '''

    ncbi_data = subparsers.add_parser('ncbi_data', description=desc_ncbi_data)

    # Create a mutually exclusive group for -i and -n
    group = ncbi_data.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--input", metavar="<FILE_In>",
                       help="a text file that contains one NCBI assembly accession per line.")
    group.add_argument("-n", "--name", metavar="<String>",
                       help="name of organism eg., 'Pseudomonas aeruginosa' etc.")
    group.add_argument("-r", "--reference", metavar="<String>",
                       choices=["bacteria", "archaea", "fungi", "invertebrate", "plant",
                                "protozoa", "vertebrate_mammalian", "vertebrate_other", "viral"],
                       help="download reference genomes (for viruses it is refseq annotated genomes) from NCBI RefSeq "
                            "subdirectory. Input should be any of them: bacteria/archaea/fungi/invertebrate/plant/"
                            "protozoa/vertebrate_mammalian/vertebrate_other/viral")
    ncbi_data.add_argument("-a", "--api-key", metavar="<String>", help="Specify an NCBI API key. "
                                                                       "To get this key kindly check: https://"
                                                                       "ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/"
                                                                       "new-api-keys-for-the-e-utilities/")
    ncbi_data.add_argument("-s", "--source", metavar="<String>", choices=['refseq', 'genbank', 'all'],
                           help="assembly source selection. Valid options: 'refseq', 'genbank' or 'all' "
                                "(default: 'refseq')")
    ncbi_data.add_argument("-f", "--format", metavar="<String>", type=parse_formats, required=True,
                           help="file format to be downloaded (comma-separated for multiple). "
                                "Valid options: 'genome', 'protein', 'cds', 'gff3', 'gtf', 'gbff', 'all'.")
    ncbi_data.add_argument("-o", "--output", metavar="<Output_Dir>", required=True,
                           help="path that will contain output.")

    # command fold parser
    desc_fold = ''' Predict structure from FASTA sequence or MSA file in A3M format. FlashFold by default uses 
    AlphaFold2 models for structure prediction. Optionally, it is possible to generate JSON file that can be 
    used as AlphaFold3 input. The JSON file will contain the FlashFold predicted MSA information, therefore the 
    homology detection steps in AlphaFold3 will be skipped. '''
    fold = subparsers.add_parser('fold', description=desc_fold)
    fold.add_argument("-q", "--query", metavar="<FILE_In|File_Dir>", required=True, type=str,
                      help="path to FASTA/A3M file(s)")
    fold.add_argument("-d", "--database", metavar="<Database_Dir>", type=str,
                      help="path to sequence database(s) created using create_db command")
    fold.add_argument("-o", "--output", metavar="<Output_Dir>", required=True, type=str,
                      help="path that will contain output")
    fold.add_argument("-t", "--threads", metavar="<Integer, >=1>", type=is_pos_int, default=16,
                      help="number of threads. Only utilized when query is a path to FASTA file(s) (default: 16)")
    fold.add_argument("--batch", action="store_true", default=False,
                      help="process multiple queries. If set, --query/-q should be the path to a "
                           "directory containing FASTA/A3M files. Parallel batch processing is available only for "
                           "homology searching. (default: False)")
    fold.add_argument("--only_msa", action="store_true", default=False,
                      help="does not predict structures but only produces MSA for given query (default: False)")
    fold.add_argument("--compact_msa", action="store_true", default=False,
                      help="creates compact MSA using CD-HIT, makes structure prediction faster with sacrificing "
                           "little or no accuracy (default: False)")
    fold.add_argument("--only_json", action="store_true", default=False,
                      help="does not predict structures but produces json file for AlphaFold3 input "
                           "for given query (default: False)")
    fold.add_argument("--num_models", type=int, choices=range(1, 6), default=5,
                      help="number of models to use for structure prediction. Reducing the number of models speeds "
                           "up the prediction but results in lower quality. (default: 5)")
    fold.add_argument("--num_recycles", metavar="<Integer, >=1>", type=is_pos_int, default=3,
                      help="number of prediction recycles. Increasing recycles can improve the prediction quality "
                           "but slows down the prediction. (default: 3)")
    fold.add_argument("--stop_at_score", metavar="<Integer>", type=int, default=100,
                      help="compute models until pLDDT (single chain) or pTM-score (multimer) > threshold is reached"
                           ". This speeds up prediction by running less models for easier queries. (default: 100) ")
    fold.add_argument("--num_model_relax", metavar="<Integer, >=0>", type=is_zero_or_pos_int, default=0,
                      help="specify how many of the top-ranked models to relax using OpenMM/Amber. Typically "
                           "relaxing the top-ranked prediction is enough and speeds up the runtime. (default: 0)")
    fold.add_argument("--relax_max_iterations", metavar="<Integer>", type=int, default=2000,
                      help="maximum number of iterations for the relaxation process. (default: 2000)")
    fold.add_argument("--overwrite_existing_results", metavar="<Boolean>", type=bool, default=False,
                      help="do not recompute results, if a query has already been predicted. (default: False)")
    fold.add_argument("--cutoff", metavar="<Float>", type=float,
                      help="Cutoff to define distances used for pDockQ2 score calculation of protein complex. "
                           "(default: 10.0)")
    fold.add_argument("--calc_extra_ptm", action="store_true", default=False,
                      help="Calculates extra PTM scores (default: False)")

    # command ligand parser
    desc_ligand = ''' Add or remove ligands in JSON files. Additionally, MSA files (A3M format) can also be used to add 
    ligands. The final output will be a JSON file, that can be used as input for AlphaFold3. '''
    ligand = subparsers.add_parser('ligand', description=desc_ligand)
    ligand.add_argument("-q", "--query", type=str, metavar="<FILE_In|File_Dir>", required=True,
                        help="path to the predicted A3M or JSON file(s)")
    ligand.add_argument("-o", "--output", type=str, metavar="<Output_Dir>", required=True,
                        help="path that will contain output")
    ligand.add_argument("-a", "--add_ligand", type=str, nargs=3, action="append",
                        metavar=("ligand_type", "ligand_name", "number_of_ligand"),
                        help="add ligand to the output JSON file for each query. Provide 'ligand type', 'ligand name', "
                             "and 'number of the ligand molecule'. The 'ligand type' must be either 'smiles' or "
                             "'ccdCodes'. Multiple ligands can be added. e.g. -a smiles CCOCCC 1 -a ccdCodes PRD 2")
    ligand.add_argument("-p", "--purge_ligands", action="store_true",
                        help="purge all ligands from the query JSON file at first. (default: False)")
    ligand.add_argument("-r", "--remove_ccdcodes", type=str, nargs="*", metavar="ccdcode",
                        help="remove ligands with ccdcodes from the query JSON file. "
                             "Multiple ccdcodes can be provided. e.g. -r PRD ATP")
    ligand.add_argument("-n", "--name", type=str, metavar="new prediction name",
                        help="Set the job name in the output JSON file. If --batch is enabled, the query name will be "
                             "used as a suffix, separated by a '-'. e.g. -n 'added_ccdCodes'")
    ligand.add_argument("-u", "--add_userccd", type=str, nargs="*", metavar="userccd_file",
                        help="add user provided ccdCodes to the output JSON file. Multiple files can be provided. "
                             "e.g. -u userccd1.cif userccd2.cif")
    ligand.add_argument("--batch", action="store_true", default=False,
                        help="process multiple queries but one-by-one. If set, --query/-q should be the path to a "
                             "directory containing A3M or JSON files. (default: False)")
    ligand.add_argument("--overwrite_existing_results", metavar="<Boolean>", type=bool, default=False,
                        help="do not recompute results, if a query has already been predicted. (default: False)")

    # command run_af3 parser
    desc_run_af3 = ''' Run AlphaFold3 prediction using the FlashFold generated JSON files. User should build a 
    Docker container with all the right python dependencies for AlphaFold3 as described in the AlphaFold3 documentation.
    For the first time execution, this script requires the local path to AlphaFold3 directory, including the path to 
    AlphaFold3 database and parameters. It can be provided again later when any of the above paths change.'''
    run_af3 = subparsers.add_parser('run_af3', description=desc_run_af3)
    run_af3.add_argument("-q", "--query", type=str, metavar="<FILE_In|File_Dir>", required=True,
                         help="path to the FlashFold generated JSON file(s)")
    run_af3.add_argument("-o", "--output", type=str, metavar="<Output_Dir>", required=True,
                         help="path that will contain output")
    run_af3.add_argument("--af3_image", type=str, metavar="<String>", default="alphafold3",
                         help="name of AlphaFold3 docker image. This is required for the first time unless the path is"
                              "changed in future. default: 'alphafold3'")
    run_af3.add_argument("--af3_db", type=Path, metavar="<File_Dir>",
                         help="path to the AlphaFold3 database directory. This is required for the first time unless "
                                "the path is changed in future.")
    run_af3.add_argument("--af3_params", type=Path, metavar="<File_Dir>",
                         help="path to the AlphaFold3 parameters directory. This is required for the first time unless "
                              "the path is changed in future.")
    run_af3.add_argument("--batch", action="store_true", default=False,
                         help="process multiple queries but one-by-one. If set, --query/-q should be the path to a "
                              "directory containing JSON files. (default: False)")
    run_af3.add_argument("--overwrite_existing_results", metavar="<Boolean>", type=bool, default=False,
                         help="do not recompute results, if a query has already been predicted. (default: False)")

    # command run_stoi parser
    desc_stoi = ''' Predict stoichiometry of a protein or protein complex. '''
    stoi = subparsers.add_parser('stoi', description=desc_stoi)
    stoi.add_argument("-q", "--query", metavar="<FILE_In>", type=str, required=True,
                          help="path to the FASTA file of the protein or protein complex")
    stoi.add_argument("-d", "--database", metavar="<Database_Dir>", type=Path, required=True,
                      help="path to FlashFold provided sequence database. ")
    stoi.add_argument("-o", "--output", metavar="<Output_Dir>", type=Path, required=True,
                          help="path to the output directory")
    stoi.add_argument("-ss", "--specific_stoichiometry", type=str, nargs=2, action="append",
                          metavar=("chain_name", "maximum_copy_number"),
                          help="Multiple specific stoichiometry can be added. Provide 'chain_name' and "
                               "'maximum_copy_number'. Chain name should be unique and should be same as FASTA header. "
                               "e.g. (If the query file has two sequences with IDs >A and >B) -ss A 2 -ss B 3")
    stoi.add_argument("-gs", "--global_stoichiometry", metavar="<Integer, >=2>", type=is_pos_int,
                            help="Maximum number of copies per chain. Should be greater than or equal to 2.")
    stoi.add_argument("-t", "--threads", metavar="<Integer, >=1>", type=is_pos_int, default=16,
                          help="number of threads (default: 16)")
    stoi.add_argument("--af3_image", type=str, metavar="<String>", default="alphafold3",
                         help="name of AlphaFold3 docker image. This is required for the first time unless the path is"
                              "changed in future. default: 'alphafold3'")
    stoi.add_argument("--af3_db", type=Path, metavar="<File_Dir>",
                         help="path to the AlphaFold3 database directory. This is required for the first time unless "
                                "the path is changed in future.")
    stoi.add_argument("--af3_params", type=Path, metavar="<File_Dir>",
                         help="path to the AlphaFold3 parameters directory. This is required for the first time unless "
                              "the path is changed in future.")
    stoi.add_argument("-af3", "--use_af3", action="store_true", default=False,
                          help="Run structure prediction with AlphaFold3 (default: False)")
    stoi.add_argument("--overwrite_existing_results", metavar="<Boolean>", type=bool, default=False,
                         help="do not recompute results, if a query has already been predicted. (default: False)")
    stoi.add_argument("--use_compact_msa", action="store_true", default=False,
                      help="creates and uses compact MSA for faster structure inference with sacrificing little or "
                           "no accuracy (default: False)")

    # command summary parser
    desc_summary = ''' Generates an interactive HTML report and a CSV file from FlashFold output. '''
    summary = subparsers.add_parser('summary', description=desc_summary)
    summary.add_argument("-d", "--directory", metavar="<File_Dir>", type=str, required=True,
                         help="path to FlashFold output directory")
    summary.add_argument("-f", "--filter", type=str, nargs=2, action="append",
                         metavar=("criteria", "score"),
                         help="filters summary report by provided criteria. Provide 'criteria' and 'score'. The "
                              "'criteria' must be any of the following: 'plddt', 'ptm', 'iptm', 'iptm_plus_ptm',"
                              " 'actifptm', 'actifptm_plus_ptm', 'min_pdockq2', 'mean_pdockq2'. Multiple filters can be"
                              " added. e.g. -f iptm_plus_ptm 0.6 -f mean_pdockq2 0.23")
    summary.add_argument("-o", "--output", metavar="<File_Dir>", type=str, required=True,
                         help="Path to the summary output directory.")
    summary.add_argument("-af3", "--alphafold3", action="store_true", default=False,
                         help="Results are from AlphaFold3 (default: False)")

    # Parse the arguments
    args = parser.parse_args()

    # noinspection PyUnreachableCode
    match args.command:
        case "create_db":
            create_protein_db_from_gbk(args)
        case "extend_db":
            extend_main_sequence_db(args)
        case "download_db":
            download_database_from_cloud(args)
        case "ncbi_data":
            download_ncbi_data(args)
        case "fold":
            predict_3d_structure(args)
        case "ligand":
            make_json_with_ligand(args)
        case "run_af3":
            process_alphafold3_task(args)
        case "stoi":
            predict_stoichiometry(args)
        case "summary":
            make_summary_report(args)
        case _:
            parser.print_help()


if __name__ == '__main__':
    main()
