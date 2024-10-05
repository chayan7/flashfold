import ftplib
import wget


def wget_file_from_url(url: str, out_dir: str) -> None:
    # Download the file to the current directory
    try:
        filename = wget.download(url, out=out_dir)
        print(f"\nFile '{filename}' downloaded successfully to {out_dir}.")
    except Exception as e:
        print(f"\nError downloading file: {e}")


def get_to_be_downloaded_files(ftp_link: str, extensions: list[str]) -> dict:
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
    except Exception as e:
        print(f"Error: {e}")
