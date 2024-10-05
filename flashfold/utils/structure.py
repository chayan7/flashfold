from .execute import run_single_job_with_progress


def get_rank_by(is_monomer: bool) -> str:
    if is_monomer:
        return "auto"
    else:
        return "multimer"


def run_colabfold(tool_path: str, is_monomer: bool, homology_search_results: str, out_path: str, num_models: int,
                  num_recycle: int, stop_at_score: int, num_top: int, relax_max_iterations: int) -> None:
    rank_by = get_rank_by(is_monomer)
    command = ("%s %s %s --rank %s --num-models %s --num-recycle %s --stop-at-score %s --num-relax %s "
               "--relax-max-iterations %s" % (tool_path, homology_search_results, out_path, rank_by, num_models,
                                              num_recycle, stop_at_score, num_top, relax_max_iterations))
    recycle = num_recycle + 1
    steps = num_models
    log_path = out_path + "/log.txt"
    run_single_job_with_progress(command, steps, recycle, "Structure prediction", log_path)
