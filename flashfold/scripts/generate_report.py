# Author: Chayan Kumar Saha

import os
import sys
import csv
import py3Dmol
import re
from typing import List, Dict, Tuple, Union, Optional
from flashfold.utils import is_valid_path

query_index = 0
length_index = 1
stoichiometry_index = 2
pLDDT_index = 3
pTM_index = 4
ipTM_index = 5
ipTM_plus_pTM_index = 6
actifpTM_index = 7
actifpTM_plus_pTM_index = 8
min_pDockQ2_index = 9
mean_pDockQ2_index = 10
model_index = 11
relaxed_model_index = 12
result_path_index = 13

summary_table_headers = [''] * 14
summary_table_headers[query_index] = 'Query'
summary_table_headers[length_index] = 'Length'
summary_table_headers[stoichiometry_index] = 'Stoichiometry'
summary_table_headers[pLDDT_index] = 'pLDDT'
summary_table_headers[ipTM_index] = 'ipTM'
summary_table_headers[pTM_index] = 'pTM'
summary_table_headers[ipTM_plus_pTM_index] = 'ipTM+pTM'
summary_table_headers[actifpTM_index] = 'actifpTM'
summary_table_headers[actifpTM_plus_pTM_index] = 'actifpTM+pTM'
summary_table_headers[min_pDockQ2_index] = 'min_pDockQ2'
summary_table_headers[mean_pDockQ2_index] = 'mean_pDockQ2'
summary_table_headers[model_index] = 'Predicted model'
summary_table_headers[relaxed_model_index] = 'Predicted model (relaxed)'
summary_table_headers[result_path_index] = 'Path to result'


row_per_page_options = [1, 5, 10, 25, 50]

current_file_dir = os.path.dirname(os.path.abspath(__file__))
static_file_dir = os.path.join(os.path.dirname(current_file_dir), "utils", "static")
scripts_js = os.path.join(static_file_dir, "scripts.js")
styles_css = os.path.join(static_file_dir, "style.css")


def remove_query_prefix(query: str) -> str:
    pattern = re.compile(r'^S\d+_')
    if pattern.match(query):
        return "_".join(query.split("_")[1:])
    return query


def make_float(input_item: Union[str, float]) -> float:
    return 0 if input_item == 'n/a' else float(input_item)


def round_if_float(input_item: Union[str, float]) -> str:
    try:
        return str(round(float(input_item), 3))
    except ValueError:
        return str(input_item)


def if_none_return_zero(input_item: Optional[float]) -> float:
    return input_item if input_item is not None else 0.0


def get_best_score_from_tsv(file_path: str) -> Dict[str, str]:
    with open(file_path, newline='') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        headers = next(reader)
        rows = [row for row in reader]
    best_score: Dict[str, str] = {}
    for row in rows:
        if "_rank_001_" in row[0]:
            for i, header in enumerate(headers):
                best_score[header] = round_if_float(row[i])
    return best_score


def get_length_stoichiometry_from_a3m(file_path: str) -> Tuple[str, str]:
    length, stoichiometry = None, None

    if not os.path.isfile(file_path):
        print(f"Error: The file '{file_path}' does not exist.")
        sys.exit()

    with open(file_path, 'r') as a3m_file:
        for line in a3m_file:
            if line.startswith('#'):
                length = line[1:].split()[0].replace(',', ':')
                stoichiometry = line[1:].split()[1].replace(',', ':')

    if length is None or stoichiometry is None:
        print(f"Error: Length and stoichiometry not found in the file '{file_path}'.")
        sys.exit()

    return length, stoichiometry


def generate_3dmol_html(pdb_path: str) -> str:
    if not is_valid_path(pdb_path):
        return 'n/a'
    
    with open(pdb_path, 'r') as f:
        pdb_data = f.read()

    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(pdb_data, 'pdb')
    viewer.setStyle({'cartoon': {'color': 'spectrum'}})
    viewer.zoomTo()
    pdb_html_path = pdb_path.replace('.pdb', '.html')
    with open(pdb_html_path, 'w') as pdb_html:
        # noinspection PyProtectedMember
        pdb_html.write(viewer._make_html())
    return pdb_html_path


def get_summary_table_rows_from_result_path(path_to_results: str) -> List[List[str]]:
    result_directory = os.path.realpath(path_to_results)
    row_of_rows = []
    for root, _, files in os.walk(result_directory):
        for file in files:
            if file == 'score.tsv':
                row = [''] * len(summary_table_headers)
                row[result_path_index] = root
                tsv_file_path = os.path.join(root, file)
                best_score_from_tsv = get_best_score_from_tsv(tsv_file_path)
                model_name = best_score_from_tsv['name']
                relaxed_model_name = model_name.replace('_unrelaxed_', '_relaxed_')
                query_id = model_name.split('_unrelaxed_rank_001_')[0]
                a3m_file = os.path.join(root, f"{query_id}.a3m")
                q_length, q_stoichiometry = get_length_stoichiometry_from_a3m(a3m_file)
                row[query_index] = query_id
                row[length_index] = q_length
                row[stoichiometry_index] = q_stoichiometry
                model_path = os.path.join(root, model_name)
                relaxed_model_path = os.path.join(root, relaxed_model_name)
                row[model_index] = generate_3dmol_html(model_path)
                row[relaxed_model_index] = generate_3dmol_html(relaxed_model_path) 
                row[pLDDT_index] = best_score_from_tsv.get('pLDDT', 'n/a')
                row[ipTM_index] = best_score_from_tsv.get('ipTM', 'n/a')
                row[pTM_index] = best_score_from_tsv.get('pTM', 'n/a')
                row[ipTM_plus_pTM_index] = best_score_from_tsv.get('ipTM+pTM', 'n/a')
                row[actifpTM_index] = best_score_from_tsv.get('actifpTM', 'n/a')
                row[actifpTM_plus_pTM_index] = best_score_from_tsv.get('actifpTM+pTM', 'n/a')
                row[min_pDockQ2_index] = best_score_from_tsv.get('min_pDockQ2', 'n/a')
                row[mean_pDockQ2_index] = best_score_from_tsv.get('mean_pDockQ2', 'n/a')
                row_of_rows.append(row)
    return row_of_rows


def remove_na_columns(headers: List[str], row_of_rows: List[List[str]]) -> Tuple[List[str], List[List[str]]]:
    """
    Removes columns where all the values are 'n/a'.

    :param headers: List of headers.
    :param row_of_rows: List of rows, where each row is a list of values.
    :return: Tuple of filtered headers and rows.
    """
    # Transpose rows to columns
    columns = list(zip(*row_of_rows))

    # Identify columns where all values are 'n/a'
    columns_to_keep = [i for i, col in enumerate(columns) if not all(value == 'n/a' for value in col)]

    # Filter headers and rows
    filtered_headers = [headers[i] for i in columns_to_keep]
    filtered_rows = [[row[i] for i in columns_to_keep] for row in row_of_rows]

    return filtered_headers, filtered_rows


def generate_html_table(headers: List[str], rows: List[List[str]], output_file: str):
    rev_model_index = headers.index('Predicted model')
    rev_relaxed_model_index = headers.index('Predicted model (relaxed)') \
        if 'Predicted model (relaxed)' in headers else -1
    rev_result_path_index = headers.index('Path to result')
    with open(output_file, 'w') as f:
        f.write('<html>\n<head>\n<title>Summary Table</title>\n')
        f.write(
            f'<link rel="stylesheet" type="text/css" href="{styles_css}">\n')
        f.write(f'<script src="{scripts_js}"></script>\n')
        f.write('</head>\n<body>\n')
        f.write('<p class="header-text">FlashFold Summary Report</p>\n')
        f.write('<input type="text" id="searchInput" onkeyup="searchTable()" placeholder="Search.. ">\n')
        f.write('<label for="rowsPerPage"> Show </label>\n')
        f.write('<select id="rowsPerPage" onchange="paginateTable()">\n')
        for option in row_per_page_options:
            selected = ' selected' if option == 10 else ''
            f.write(f'<option value="{option}"{selected}>{option}</option>\n')
        f.write('</select>\n')
        f.write('<label for="rowsPerPage"> prediction(s) </label>\n')
        f.write('<table id="summaryTable">\n')
        f.write('<thead>\n<tr>\n')
        for header in headers[:rev_model_index]:
            f.write(f'<th onclick="sortTable({headers.index(header)})">{header}</th>\n')
        for rest_header in headers[rev_model_index:]:
            f.write(f'<th>{rest_header}</th>\n')
        f.write('</tr>\n</thead>\n<tbody>\n')
        for row in rows:
            f.write('<tr>\n')
            for cell in row[:rev_model_index]:
                f.write(f'<td>{cell}</td>\n')
            f.write(
                f'<td><button onclick="visualizePDB(\'{row[rev_model_index]}\')">Show structure</button></td>\n')
            if rev_relaxed_model_index != -1:
                if row[rev_relaxed_model_index] != 'n/a':
                    f.write(f'<td><button onclick="visualizePDB(\'{row[rev_relaxed_model_index]}\')">'
                            f'Show structure</button></td>\n')
                else:
                    f.write(f'<td>{row[rev_relaxed_model_index]}</td>\n')
            f.write(
                f'<td><button onclick="openInFolder(\'{row[rev_result_path_index]}\')">Open</button></td>\n')
            f.write('</tr>\n')
        f.write('</tbody>\n</table>\n')
        f.write('<div class="flex-container">\n')
        f.write('<div class="entries-info" id="entriesInfo"></div>\n')
        f.write('<div id="pagination"></div>\n')
        f.write('</div>\n')
        f.write('</body>\n</html>\n')


def generate_csv_table(headers: List[str], rows: List[List[str]], output_file: str) -> None:
    with open(output_file, 'w', newline='') as csv_out:
        # noinspection PyTypeChecker
        writer = csv.writer(csv_out, delimiter='\t')
        writer.writerow(headers)
        for row in rows:
            writer.writerow(row)


def make_summary_report(args) -> None:
    if not is_valid_path(args.directory):
        print(f"\n-- Error: The directory '{args.directory}' does not exist. Input a valid directory.")
        return

    filter_dict = dict()
    filter_dict[pLDDT_index] = if_none_return_zero(args.filter_by_plddt)
    filter_dict[pTM_index] = if_none_return_zero(args.filter_by_ptm)
    filter_dict[ipTM_index] = if_none_return_zero(args.filter_by_iptm)
    filter_dict[ipTM_plus_pTM_index] = if_none_return_zero(args.filter_by_iptm_plus_ptm)
    filter_dict[actifpTM_index] = if_none_return_zero(args.filter_by_actifptm)
    filter_dict[actifpTM_plus_pTM_index] = if_none_return_zero(args.filter_by_actifptm_plus_ptm)
    filter_dict[min_pDockQ2_index] = if_none_return_zero(args.filter_by_min_pdockq2)
    filter_dict[mean_pDockQ2_index] = if_none_return_zero(args.filter_by_avg_pdockq2)

    summary_table_rows = get_summary_table_rows_from_result_path(args.directory)
    if len(summary_table_rows) == 0:
        print(f"\n-- Error: No results found in the provided path below:\n\t'{os.path.realpath(args.directory)}'\n")
        return

    clean_sum_tab_headers, clean_sum_tab_rows = remove_na_columns(summary_table_headers, summary_table_rows)
    output_html_file_path = os.path.join(os.path.realpath(args.output), 'summary.html')
    output_csv_file_path = os.path.join(os.path.realpath(args.output), 'summary.csv')
    generate_html_table(clean_sum_tab_headers, clean_sum_tab_rows, output_html_file_path)
    generate_csv_table(clean_sum_tab_headers, clean_sum_tab_rows, output_csv_file_path)
    print(f"-- HTML report has been generated: {output_html_file_path}")
    print(f"-- CSV report has been generated: {output_csv_file_path}")

    if all([filter_dict[i] == 0 for i in filter_dict]):
        return

    filtered_table_rows = [row for row in summary_table_rows if
                           all([make_float(row[i]) >= filter_dict[i] for i in filter_dict])]

    if len(filtered_table_rows) == 0:
        print(f"\n-- Warning: No results met the filtering criteria. \n")
        return

    clean_sum_tab_headers, clean_filtered_tab_rows = remove_na_columns(summary_table_headers, filtered_table_rows)

    filtered_output_html_file_path = os.path.join(os.path.realpath(args.output), 'summary_filtered.html')
    filtered_output_csv_file_path = os.path.join(os.path.realpath(args.output), 'summary_filtered.csv')
    generate_html_table(clean_sum_tab_headers, clean_filtered_tab_rows, filtered_output_html_file_path)
    print(f"-- Filtered HTML report has been generated: {filtered_output_html_file_path}")
    generate_csv_table(clean_sum_tab_headers, clean_filtered_tab_rows, filtered_output_csv_file_path)
    print(f"-- Filtered CSV report has been generated: {filtered_output_csv_file_path}")
    return
