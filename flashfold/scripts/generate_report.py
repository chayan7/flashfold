# Author: Chayan Kumar Saha

import os
import sys
import csv
import py3Dmol
import re
from typing import List, Dict, Tuple, Union, Optional, Literal
from flashfold.utils import is_valid_path, load_json_file

query_index = 0
length_index = 1
stoichiometry_index = 2
plddt_index = 3
ptm_index = 4
iptm_index = 5
iptm_plus_ptm_index = 6
actifptm_index = 7
actifptm_plus_ptm_index = 8
min_pdockq2_index = 9
mean_pdockq2_index = 10
ranking_index = 11
model_index = 12
relaxed_model_index = 13
result_path_index = 14

summary_table_headers = [''] * 15
summary_table_headers[query_index] = 'Query'
summary_table_headers[length_index] = 'Length'
summary_table_headers[stoichiometry_index] = 'Stoichiometry'
summary_table_headers[plddt_index] = 'pLDDT'
summary_table_headers[iptm_index] = 'ipTM'
summary_table_headers[ptm_index] = 'pTM'
summary_table_headers[iptm_plus_ptm_index] = 'ipTM+pTM'
summary_table_headers[actifptm_index] = 'actifpTM'
summary_table_headers[actifptm_plus_ptm_index] = 'actifpTM+pTM'
summary_table_headers[min_pdockq2_index] = 'min_pDockQ2'
summary_table_headers[mean_pdockq2_index] = 'mean_pDockQ2'
summary_table_headers[ranking_index] = 'Ranking_score'
summary_table_headers[model_index] = 'Predicted model'
summary_table_headers[relaxed_model_index] = 'Predicted model (relaxed)'
summary_table_headers[result_path_index] = 'Path to result'


row_per_page_options = [1, 5, 10, 25, 50]

current_file_dir = os.path.dirname(os.path.abspath(__file__))
static_file_dir = os.path.join(os.path.dirname(current_file_dir), "utils", "static")
scripts_js = os.path.join(static_file_dir, "scripts.js")
styles_css = os.path.join(static_file_dir, "style.css")


def make_float(input_item: Union[str, float]) -> float:
    return 0 if input_item == 'n/a' else float(input_item)


def round_if_float(input_item: Union[str, float]) -> str:
    try:
        return str(round(float(input_item), 3))
    except ValueError:
        return str(input_item)


def return_float_if_float(input_name: str, input_score: str) -> float:
    if not input_score.replace('.', '', 1).isdigit():
        print(f"\n-- Warning: Please provide a valid '{input_name}' score for filtering. "
              f"Provided score: '{input_score}'\n")
        sys.exit()
    else:
        return float(input_score)



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


def generate_3dmol_html(file_path: str, file_type: Literal['pdb', 'cif']) -> str:
    if not is_valid_path(file_path):
        return 'n/a'

    with open(file_path, 'r') as f:
        file_data = f.read()

    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(file_data, file_type)
    viewer.setStyle({'cartoon': {'color': 'spectrum'}})
    viewer.zoomTo()
    text_to_replace = f".{file_type}"
    html_path = file_path.replace(text_to_replace, '.html')
    with open(html_path, 'w') as html_file:
        # noinspection PyProtectedMember
        html_file.write(viewer._make_html())
    return html_path


def get_summary_table_rows_from_result_path(path_to_results: str, is_af3: bool) -> List[List[str]]:
    result_directory = os.path.realpath(path_to_results)
    row_of_rows = []
    for root, _, files in os.walk(result_directory):
        for file in files:
            if not is_af3:
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
                    row[model_index] = generate_3dmol_html(model_path, 'pdb')
                    row[relaxed_model_index] = generate_3dmol_html(relaxed_model_path, 'pdb')
                    row[plddt_index] = best_score_from_tsv.get('pLDDT', 'n/a')
                    row[iptm_index] = best_score_from_tsv.get('ipTM', 'n/a')
                    row[ptm_index] = best_score_from_tsv.get('pTM', 'n/a')
                    row[iptm_plus_ptm_index] = best_score_from_tsv.get('ipTM+pTM', 'n/a')
                    row[actifptm_index] = best_score_from_tsv.get('actifpTM', 'n/a')
                    row[actifptm_plus_ptm_index] = best_score_from_tsv.get('actifpTM+pTM', 'n/a')
                    row[min_pdockq2_index] = best_score_from_tsv.get('min_pDockQ2', 'n/a')
                    row[mean_pdockq2_index] = best_score_from_tsv.get('mean_pDockQ2', 'n/a')
                    row[ranking_index] = 'n/a'
                    row_of_rows.append(row)
            else:
                if file.endswith('_summary_confidences.json'):
                    row = [''] * len(summary_table_headers)
                    json_file_path = os.path.join(root, file)
                    query_id = file.split("_summary_confidences.json")[0]
                    cif_file = os.path.join(root, f"{query_id}_model.cif")
                    loaded_json_data = load_json_file(json_file_path)
                    row[query_index] = query_id
                    row[length_index] = loaded_json_data.get('length', 'n/a')
                    row[stoichiometry_index] = loaded_json_data.get('stoichiometry', 'n/a')
                    row[plddt_index] = loaded_json_data.get('pLDDT', 'n/a')
                    iptm_score = loaded_json_data.get('iptm')
                    ptm_score = loaded_json_data.get('ptm')
                    iptm_plus_ptm_score = 'n/a' if not iptm_score or not ptm_score else \
                        (0.8 * float(iptm_score) + 0.2 * float(ptm_score))
                    row[iptm_index] = 'n/a' if not iptm_score else round_if_float(iptm_score)
                    row[ptm_index] = 'n/a' if not ptm_score else round_if_float(ptm_score)
                    row[iptm_plus_ptm_index] = round_if_float(iptm_plus_ptm_score)
                    row[actifptm_index] = 'n/a'
                    row[actifptm_plus_ptm_index] = 'n/a'
                    row[min_pdockq2_index] = 'n/a'
                    row[mean_pdockq2_index] = 'n/a'
                    ranking_score = loaded_json_data.get('ranking_score')
                    row[ranking_index] = 'n/a' if not ranking_score else round_if_float(ranking_score)
                    row[model_index] = generate_3dmol_html(cif_file, 'cif')
                    row[relaxed_model_index] = "n/a"
                    row[result_path_index] = root
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
    if args.filter:
        for f_name, f_score in args.filter:
            index_name = f"{f_name.lower()}_index"
            if index_name not in globals():
                print(f"\n-- Warning: Skipping the invalid filter '{f_name}'. Please use one of the following:\n\t"
                      f"'plddt', 'ptm', 'iptm', 'iptm_plus_ptm', 'actifptm', 'actifptm_plus_ptm', "
                      f"'min_pdockq2', 'mean_pdockq2'")
                continue
            else:
                filter_dict[globals()[index_name]] = return_float_if_float(f_name, f_score)

    summary_table_rows = get_summary_table_rows_from_result_path(args.directory, args.alphafold3)
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
