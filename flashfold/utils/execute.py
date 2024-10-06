import threading
import subprocess
import sys
from queue import Queue, Empty
from tqdm import tqdm
from typing import List
from .util import current_time


def worker_func(q, stopped):
    """
    Worker function to process commands from a queue until a stop event is set.

    Args:
        q (queue.Queue): A thread-safe queue containing commands to be executed.
        stopped (threading.Event): An event to signal when the worker should stop processing.

    The function retrieves commands from the queue using the get_nowait() method to avoid blocking.
    It executes each command using subprocess.run() and handles any CalledProcessError exceptions.
    The function continues processing until the stop event is set and the queue is empty.
    """    
    while not stopped.is_set() or not q.empty():
        try:
            # Use the get_nowait() method for retrieving a queued item to
            # prevent the thread from blocking when the queue is empty
            com = q.get_nowait()
        except Empty:
            continue
        try:
            subprocess.run(com, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print('Error running command:', str(e))
        finally:
            q.task_done()


def run_jobs_in_parallel(thread_count: int, jobs: List, job_name: str) -> None:
    """
        Args:
            thread_count (int): The number of threads to use for running the jobs.
            jobs (list): A list of jobs (commands) to be executed.
            job_name (str): A name for the job batch, used for logging purposes.

        Returns:
            None
    """
    stopped = threading.Event()
    q = Queue()
    if job_name != "":
        print(f'\n-- {current_time()} > {job_name} is being processed with {thread_count} thread limit')

    for command in jobs:
        q.put(command)

    threads = []
    for _ in range(thread_count):
        t = threading.Thread(target=worker_func, args=(q, stopped))
        t.start()
        threads.append(t)

    q.join()  # Block until all tasks are done
    stopped.set()

    # Ensure all threads have finished
    for t in threads:
        t.join()
    if job_name != "":
        print(f'-- {current_time()} > {job_name} has been completed \n')


def run_single_job(job: str, job_name: str) -> None:
    """
    Run a single job using subprocess.run() and print the start and completion messages.
    """
    if job.strip() == "":
        print("Warning: No executable command found.")
        sys.exit()
    if job_name.strip() == "":
        print("Warning: No job name found.")
        sys.exit()

    print(f'\n-- {current_time()} > {job_name} is being processed')
    try:
        subprocess.run(job, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print('Error running command:', str(e))
    print(f'-- {current_time()} > {job_name} has been completed \n')


def get_count_from_log(log_file_path: str) -> int:
    """
    Counts occurrences of a specific pattern in the log file.
    """
    count = 0
    retry = True
    while retry:
        try:
            with open(log_file_path, 'r', encoding='utf-8') as file:
                for line in file:
                    if "recycle=" in line:
                        count += 1
            retry = False
        except FileNotFoundError:
            retry = True
    return count


def get_updated_count_from_log(log_file_path: str, prev_count: int, total: int) -> int:
    """
    Get the updated count from the log file.
    Args:
        log_file_path (str): The path to the log file.
        prev_count (int): The previous count of occurrences.
        total (int): The total number of occurrences to reach.

    Returns:
        int: The updated count of occurrences.
    """
    if prev_count == total:
        return total

    retry = True
    new_count = 0
    while retry:
        latest_count = get_count_from_log(log_file_path)
        if latest_count > prev_count:
            new_count = latest_count
            retry = False
    return new_count


def run_single_job_with_progress(job: str, steps: int, intervals: int, job_name: str, log_path: str) -> None:
    """
    Run a single job using subprocess.Popen() and display a progress bar based on the log file.
    Args:
        job: The command to run.
        steps: The number of steps to complete the job.
        intervals: The number of intervals to update the progress bar.
        job_name: The name of the job.
        log_path: The path to the log file.

    Returns:
        None

    """
    if job.strip() == "":
        print("Warning: No executable command found.")
        sys.exit()
    if job_name.strip() == "":
        print("Warning: No job name found.")
        sys.exit()

    # Use subprocess.Popen to run the command
    try:
        print(f'\n-- {current_time()} > {job_name} is being processed')
        total_count = steps*intervals
        description = "\tPrediction-recycle progress"
        progress_bar = tqdm(total=steps, desc=description, unit="%", ncols=80, delay=1,
                            bar_format="{desc}: |{bar}| {percentage:3.0f}% ...  ", colour='#D3D3D3')
        # Start the subprocess
        process = subprocess.Popen(job, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, shell=True)
        old_count = get_count_from_log(log_path)
        ret_code = process.poll()
        retry = ret_code is None  # None means process is still running
        while retry:
            if ret_code is not None:  # Process has finished
                retry = False  # Process has finished, exit loop

            new_count = get_updated_count_from_log(log_path, old_count, total_count)
            if new_count == total_count:
                retry = False

            if new_count > old_count:
                percentage_new = round(new_count/total_count*steps, 2)
                percentage_old = round(old_count/total_count*steps, 2)
                progress_bar.update(percentage_new - percentage_old)  # Update the progress bar with new occurrences
                old_count = new_count  # Update the old count
                retry = True
        progress_bar.close()

    except subprocess.CalledProcessError as e:
        print('Error running command:', str(e))

    print(f'-- {current_time()} > {job_name} has been completed \n')
