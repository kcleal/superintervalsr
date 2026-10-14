#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Developer dependencies: Rcpp, roxygen2, testthat, knitr, rmarkdown, BiocCheck.
Rscript -e 'Rcpp::compileAttributes("."); roxygen2::roxygenise(".")'
Rscript -e 'testthat::test_local(".")'
R CMD build .
package_version=$(Rscript -e 'cat(read.dcf("DESCRIPTION")[1, "Version"])')
tarball="superintervalsr_${package_version}.tar.gz"
check_args=(--no-manual)
if [[ "${1:-}" == "--manual" ]]; then
    check_args=()
fi
R CMD check "${check_args[@]}" "$tarball"
Rscript -e '
if (!requireNamespace("BiocCheck", quietly = TRUE)) {
    stop("Install BiocCheck in your Bioconductor development library")
}
checks <- BiocCheck::BiocCheck(commandArgs(TRUE)[1])
if (checks$getNum("error") || checks$getNum("warning")) {
    stop("BiocCheck reported errors or warnings; inspect its output")
}
' "$tarball"
