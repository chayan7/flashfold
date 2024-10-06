import argparse
from flashfold.scripts import create_protein_db_from_gbk, download_representative_db, predict_3d_structure, \
    extend_main_sequence_db
from flashfold.utils import is_zero_or_pos_int, is_pos_int


def main() -> None:
    # Create the main parser
    parser = argparse.ArgumentParser(description='Predict protein structure from sequence.')

    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', title='commands',
                                       description='Choose any of the following command options')

    # command create_db parser
    desc_create_db = ''' Create database from Genbank files, by extracting protein sequences. 
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
    # # Create a mutually exclusive group for -n and -g
    group = extend_db.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--genbank_path", metavar="<Genbank_Dir>",
                       help="path to the directory that contains genbank file(s), will be used to "
                            "update the main database")
    group.add_argument("-n", "--new_db", metavar="<Database_Dir>",
                       help="path to the new database that will be added to the main database")

    # command download_db parser
    desc_download_db = ''' Download representative database. '''
    download_db = subparsers.add_parser('download_db', description=desc_download_db)
    download_db.add_argument("-p", "--path", metavar="<Output_Dir>", required=True,
                             help="path to store the database built from representative sequence")

    # command predict_3d parser
    desc_fold = ''' Predict structure from FASTA sequence. '''
    fold = subparsers.add_parser('fold', description=desc_fold)
    fold.add_argument("-q", "--query", metavar="<FILE_In>", required=True, help="path to FASTA file")
    fold.add_argument("-d", "--database", metavar="<Database_Dir>", required=True,
                         help="path to sequence database(s) created using create_db command")
    fold.add_argument("-o", "--output", metavar="<Output_Dir>", required=True,
                         help="path that will contain output")
    fold.add_argument("-t", "--threads", metavar="<Integer, >=1>", type=is_pos_int, default=16,
                         help="number of threads. (default: 16)")
    fold.add_argument("--num_models", type=int, choices=range(1, 6), default=5,
                         help="number of models to use for structure prediction. Reducing the number of models speeds "
                              "up the prediction but results in lower quality. (default: 5)")
    fold.add_argument("--num_recycle", metavar="<Integer, >=1>", type=is_pos_int, default=3,
                         help="number of prediction recycles. Increasing recycles can improve the prediction quality "
                              "but slows down the prediction. (default: 3)")
    fold.add_argument("--stop_at_score", metavar="<Integer>", type=int, default=100,
                         help="compute models until pLDDT (single chain) or pTM-score (multimer) > threshold is reached"
                              ". This speeds up prediction by running less models for easier queries. (default: 100) ")
    fold.add_argument("--num_structure_relax", metavar="<Integer, >=0>", type=is_zero_or_pos_int, default=0,
                         help="specify how many of the top-ranked structures to relax using OpenMM/Amber. Typically "
                              "relaxing the top-ranked prediction is enough and speeds up the runtime. (default: 0)")
    fold.add_argument("--relax_max_iterations", metavar="<Integer>", type=int, default=2000,
                         help="maximum number of iterations for the relaxation process. (default: 2000)")
    fold.add_argument("--overwrite_existing_results", metavar="<Boolean>", type=bool, default=False,
                         help="do not recompute results, if a query has already been predicted. (default: False)")
    fold.add_argument("-c", "--cutoff", metavar="<Float>", type=float, default=10.0,
                         help="Cutoff to define distances used for pDockQ2 score calculation of protein complex. "
                              "(default: 10.0)")

    # Parse the arguments
    args = parser.parse_args()

    if args.command == "create_db":
        create_protein_db_from_gbk(args.output, args.path)
    elif args.command == "extend_db":
        extend_main_sequence_db(args.main_db, args.genbank_path, args.new_db)
    elif args.command == "download_db":
        download_representative_db(args.path)
    elif args.command == "fold":
        predict_3d_structure(args.query, args.database, args.output, args.threads, args.num_models,
                             args.num_recycle, args.stop_at_score, args.num_structure_relax, args.relax_max_iterations,
                             args.overwrite_existing_results, args.cutoff)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
