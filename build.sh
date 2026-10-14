#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Developer dependencies: Rcpp, roxygen2, testthat, knitr, rmarkdown, BiocCheck.
Rscript -e 'Rcpp::compileAttributes("."); roxygen2::roxygenise(".")'
Rscript -e 'testthat::test_local(".")'
R CMD build .
package_version=$(Rscript -e 'cat(read.dcf("DESCRIPTION")[1, "Version"])')
tarball="superintervalsr_${package_version}.tar.gz"
R CMD check --no-manual "$tarball"
Rscript -e 'if (!requireNamespace("BiocCheck", quietly = TRUE)) stop("Install BiocCheck in your Bioconductor development library"); BiocCheck::BiocCheck(commandArgs(TRUE)[1])' "$tarball"
