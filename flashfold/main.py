import argparse
from flashfold.scripts import create_protein_db_from_gbk, download_database_from_cloud, predict_3d_structure, \
    extend_main_sequence_db, download_ncbi_data, parse_formats, make_summary_report
from flashfold.utils import is_zero_or_pos_int, is_pos_int


def main() -> None:
    # Create the main parser
    parser = argparse.ArgumentParser(description='Predict protein structure from sequence.')
    parser.add_argument("-v", "--version", action="version", version="∞∞ FlashFold v1.1.3 ∞∞")

    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', title='commands',
                                       description='Choose any of the following command options')

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
    desc_fold = ''' Predict structure from FASTA sequence. '''
    fold = subparsers.add_parser('fold', description=desc_fold)
    fold.add_argument("-q", "--query", metavar="<FILE_In|File_Dir>", required=True,
                      help="path to FASTA/A3M file(s)")
    fold.add_argument("-d", "--database", metavar="<Database_Dir>",
                      help="path to sequence database(s) created using create_db command")
    fold.add_argument("-o", "--output", metavar="<Output_Dir>", required=True,
                      help="path that will contain output")
    fold.add_argument("-t", "--threads", metavar="<Integer, >=1>", type=is_pos_int, default=16,
                      help="number of threads. Only utilized when query is a path to FASTA file(s) (default: 16)")
    fold.add_argument("--batch", action="store_true", default=False,
                      help="process multiple queries (default: False). If set, --query/-q should be the path to a "
                           "directory containing FASTA files.")
    fold.add_argument("--only_msa", action="store_true", default=False,
                      help="does not predict structures but only produces MSA for given query (default: False)")
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

    # command summary parser
    desc_summary = ''' Generates an interactive HTML report and a CSV file from FlashFold output. '''
    summary = subparsers.add_parser('summary', description=desc_summary)
    summary.add_argument("-d", "--directory", metavar="<File_Dir>", type=str, required=True,
                         help="path to FlashFold output directory")
    summary.add_argument("-fl", "--filter_by_plddt", metavar="<Float>", type=float,
                         help="Filter output by pLDDT score.")
    summary.add_argument("-fp", "--filter_by_ptm", metavar="<Float>", type=float,
                         help="Filter output by pTM score.")
    summary.add_argument("-fi", "--filter_by_iptm", metavar="<Float>", type=float,
                         help="Filter output by ipTM score.")
    summary.add_argument("-fip", "--filter_by_iptm_plus_ptm", metavar="<Float>", type=float,
                         help="Filter output by ipTM+pTM score.")
    summary.add_argument("-fai", "--filter_by_actifptm", metavar="<Float>", type=float,
                         help="Filter output by actifpTM score.")
    summary.add_argument("-faip", "--filter_by_actifptm_plus_ptm", metavar="<Float>", type=float,
                         help="Filter output by actifpTM+pTM score.")
    summary.add_argument("-fmd", "--filter_by_min_pdockq2", metavar="<Float>", type=float,
                         help="Filter output by minimum pDockQ2 score.")
    summary.add_argument("-fad", "--filter_by_avg_pdockq2", metavar="<Float>", type=float,
                         help="Filter output by average pDockQ2 score.")
    summary.add_argument("-o", "--output", metavar="<File_Dir>", type=str, required=True,
                         help="Path to the summary output directory.")

    # Parse the arguments
    args = parser.parse_args()

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
        case "summary":
            make_summary_report(args)
        case _:
            parser.print_help()


if __name__ == '__main__':
    main()
