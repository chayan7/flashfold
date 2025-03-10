from .execute import run_single_job_with_progress


def get_rank_by(is_monomer: bool) -> str:
    """
    Get the rank by parameter for the structure prediction.
    Args:
        is_monomer:     A boolean value indicating if the sequence is a monomer or not.

    Returns:
        str: The rank by parameter.
    """
    if is_monomer:
        return "auto"
    else:
        return "multimer"


def run_colabfold(is_monomer: bool, homology_search_results: str, out_path: str, num_models: int, num_recycle: int,
                  stop_at_score: int, num_relax: int, relax_max_iterations: int, calc_ext_ptm: bool) -> None:
    """
    Run ColabFold for structure prediction.
    """
    rank_by = get_rank_by(is_monomer)
    add_calc_ext_ptm = "" if not calc_ext_ptm else "--calc-extra-ptm"
    command = ("colabfold_batch %s %s --rank %s --num-models %s --num-recycle %s --stop-at-score %s --num-relax %s "
               "--relax-max-iterations %s %s" % (homology_search_results, out_path, rank_by, num_models, num_recycle,
                                                 stop_at_score, num_relax, relax_max_iterations, add_calc_ext_ptm))
    recycle = num_recycle + 1
    steps = num_models
    log_path = out_path + "/log.txt"
    run_single_job_with_progress(command, steps, recycle, "Structure prediction", log_path)
