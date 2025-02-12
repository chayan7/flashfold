import wget
import sys
import shutil


def wget_file_from_url(url: str, out_dir: str) -> str:
    """
    Downloads a file from the specified URL to the given output directory.

    Args:
        url (str): The URL of the file to download.
        out_dir (str): The directory where the downloaded file will be saved.

    Raises:
        wget.WgetError: If there is an error during the download process.
        OSError: If there is an OS-related error during the download process.

    Returns:
        str: The path to the downloaded file.

    Example:
        wget_file_from_url("https://example.com/file.zip", "/path/to/output/dir")
    """
    try:
        file_path = wget.download(url, out=out_dir)
        print(f"\nFile downloaded successfully: '{file_path}'\n")
        return file_path
    except OSError as e:  # Handle any OS-related errors (e.g., file write issues)
        print(f"\nError downloading file: {e}")
        shutil.rmtree(out_dir)
        sys.exit()
    except Exception as e:  # Handle any other general errors like connection issues
        print(f"\nAn unexpected error occurred: {e}")
        shutil.rmtree(out_dir)
        sys.exit()


