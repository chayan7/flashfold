/**
 * Author: 	Chayan Kumar Saha
 */

/**
 * Sorts the table by the specified column.
 * @param n
 */
function sortTable(n) {
    let table, rows, switching, i, x, y, shouldSwitch, dir, switchcount = 0;
    table = document.getElementById("summaryTable");
    switching = true;
    dir = "asc";
    while (switching) {
        switching = false;
        rows = table.rows;
        for (i = 1; i < (rows.length - 1); i++) {
            shouldSwitch = false;
            x = rows[i].getElementsByTagName("TD")[n];
            y = rows[i + 1].getElementsByTagName("TD")[n];
            if (dir === "asc") {
                if (x.innerHTML.toLowerCase() > y.innerHTML.toLowerCase()) {
                    shouldSwitch = true;
                    break;
                }
            } else if (dir === "desc") {
                if (x.innerHTML.toLowerCase() < y.innerHTML.toLowerCase()) {
                    shouldSwitch = true;
                    break;
                }
            }
        }
        if (shouldSwitch) {
            rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
            switching = true;
            switchcount++;
        } else {
            if (switchcount == 0 && dir == "asc") {
                dir = "desc";
                switching = true;
            }
        }
    }
}

/**
 * Filters the table based on the search input.
 */
function searchTable() {
    let input, filter, table, tr, td, i, j, txtValue;
    input = document.getElementById("searchInput");
    filter = input.value.toLowerCase();
    table = document.getElementById("summaryTable");
    tr = table.getElementsByTagName("tr");
    for (i = 1; i < tr.length; i++) {
        tr[i].style.display = "none";
        td = tr[i].getElementsByTagName("td");
        for (j = 0; j < td.length; j++) {
            if (td[j]) {
                txtValue = td[j].textContent || td[j].innerText;
                if (txtValue.toLowerCase().indexOf(filter) > -1) {
                    tr[i].style.display = "";
                    break;
                }
            }
        }
    }
}

/**
 * Paginates the table.
 */
function paginateTable() {
    let table, rows, pagination, rowsPerPage, pageCount, currentPage;
    table = document.getElementById("summaryTable");
    rows = table.getElementsByTagName("tr");
    pagination = document.getElementById("pagination");
    rowsPerPage = parseInt(document.getElementById("rowsPerPage").value);
    pageCount = Math.ceil((rows.length - 1) / rowsPerPage);
    currentPage = 1;

    /**
     * Shows the specified page.
     * @param page
     */
    function showPage(page) {
        let start = (page - 1) * rowsPerPage + 1;
        let end = start + rowsPerPage;
        for (let i = 1; i < rows.length; i++) {
            if (i >= start && i < end) {
                rows[i].style.display = "";
            } else {
                rows[i].style.display = "none";
            }
        }
        updatePagination(page);
        updateEntriesInfo(start, Math.min(end - 1, rows.length - 1), rows.length - 1);
    }


    /**
     * Updates the pagination.
     * @param page
     */
    function updatePagination(page) {
        pagination.innerHTML = "";

        // Previous button
        let prevButton = document.createElement("a");
        prevButton.href = "#";
        prevButton.innerHTML = "Previous";
        prevButton.className = (page === 1) ? "disabled" : "";
        prevButton.onclick = function() {
            if (page > 1) showPage(page - 1);
        };
        pagination.appendChild(prevButton);

        // Page numbers
        let maxPagesToShow = 3;
        let startPage = Math.max(1, page - Math.floor(maxPagesToShow / 2));
        let endPage = Math.min(pageCount, startPage + maxPagesToShow - 1);

        if (startPage > 1) {
            let firstPageLink = document.createElement("a");
            firstPageLink.href = "#";
            firstPageLink.innerHTML = "1";
            firstPageLink.onclick = function() {
                showPage(1);
            };
            pagination.appendChild(firstPageLink);

            if (startPage > 2) {
                let ellipsis = document.createElement("span");
                ellipsis.innerHTML = "...";
                pagination.appendChild(ellipsis);
            }
        }

        for (let i = startPage; i <= endPage; i++) {
            let pageLink = document.createElement("a");
            pageLink.href = "#";
            pageLink.innerHTML = i;
            pageLink.className = (i === page) ? "active" : "";
            pageLink.onclick = (function(i) {
                return function() {
                    showPage(i);
                };
            })(i);
            pagination.appendChild(pageLink);
        }

        if (endPage < pageCount) {
            if (endPage < pageCount - 1) {
                let ellipsis = document.createElement("span");
                ellipsis.innerHTML = "...";
                pagination.appendChild(ellipsis);
            }

            let lastPageLink = document.createElement("a");
            lastPageLink.href = "#";
            lastPageLink.innerHTML = pageCount;
            lastPageLink.onclick = function() {
                showPage(pageCount);
            };
            pagination.appendChild(lastPageLink);
        }

        // Next button
        let nextButton = document.createElement("a");
        nextButton.href = "#";
        nextButton.innerHTML = "Next";
        nextButton.className = (page === pageCount) ? "disabled" : "";
        nextButton.onclick = function() {
            if (page < pageCount) showPage(page + 1);
        };
        pagination.appendChild(nextButton);
    }


    /**
     * Updates the entries info.
     * @param start
     * @param end
     * @param total
     */
    function updateEntriesInfo(start, end, total) {
        let entriesInfo = document.getElementById("entriesInfo");
        entriesInfo.innerHTML = `Showing ${start} to ${end} of ${total} prediction(s)`;
    }

    showPage(currentPage);
}


/**
 * Visualizes the specified PDB file.
 * @param path
 */
function visualizePDB(path) {
    const windowFeatures = "left=100,top=100,width=720,height=640";
    window.open(
        path,
        "top",
        windowFeatures,
    );
}


/**
 * Opens the specified folder in the file manager.
 * @param path
 */
function openInFolder(path) {
    const windowFeatures = "left=100,top=100,width=1060,height=640";
    window.open(
        'file://' + path,
        "top",
        windowFeatures,
    );
}


window.onload = function() {
    paginateTable();
};