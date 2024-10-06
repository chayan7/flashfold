import ftplib
import wget
from typing import List, Dict


def wget_file_from_url(url: str, out_dir: str) -> None:
    """
    Downloads a file from the specified URL to the given output directory.

    Args:
        url (str): The URL of the file to download.
        out_dir (str): The directory where the downloaded file will be saved.

    Raises:
        wget.WgetError: If there is an error during the download process.
        OSError: If there is an OS-related error during the download process.

    Example:
        wget_file_from_url("https://example.com/file.zip", "/path/to/output/dir")
    """
    try:
        filename = wget.download(url, out=out_dir)
        print(f"\nFile '{filename}' downloaded successfully to {out_dir}.")
    except (wget.WgetError, OSError) as e:
        print(f"\nError downloading file: {e}")


def get_to_be_downloaded_files(ftp_link: str, extensions: List[str]) -> Dict:
    """
    Retrieves a dictionary of files to be downloaded from an FTP server based on specified file extensions.

    Args:
        ftp_link (str): The FTP link to the directory containing the files.
        extensions (List[str]): A list of file extensions to filter the files to be downloaded.

    Returns:
        Dict: A dictionary where the keys are filenames and the values are the corresponding download links.

    Raises:
        ftplib.all_errors: If any FTP-related error occurs.
        OSError: If any OS-related error occurs.
    """
    files_to_download = {}
    try:
        ftp_split_directory = ftp_link.split('/')[3:]
        ftp_path = '/'.join(map(str, ftp_split_directory))
        ftp = ftplib.FTP('ftp.ncbi.nih.gov', 'anonymous', 'anonymous@ftp.ncbi.nih.gov')
        ftp.cwd('/' + ftp_path)
        files = ftp.nlst()
        for filename in files:
            for extension in extensions:
                if extension in filename:
                    files_to_download[filename] = f"{ftp_link}/{filename}"
        return files_to_download
    except (ftplib.all_errors, OSError) as e:
        print(f"Error: {e}")
