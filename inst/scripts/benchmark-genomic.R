# GRanges interface benchmark: real airway annotations and synthetic subjects.
# Run: Rscript benchmark-genomic.R [output-directory]
# The package's vignettes explain biological context and coordinate semantics.
# Assisted by OpenAI Codex (October 2026).
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(superintervalsr)
})
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1] else "genomic-benchmark"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data("airway", package = "airway")
models <- SummarizedExperiment::rowRanges(airway)
spans <- range(models)
single_span <- lengths(spans) == 1L
airway_genes <- unlist(spans[single_span], use.names = FALSE)
names(airway_genes) <- names(models)[single_span]
airway_genes <- airway_genes[width(airway_genes) > 0L]

median_time <- function(fun, repetitions = 5L, calls = 10L) {
    median(vapply(seq_len(repetitions), function(i) {
        invisible(gc())
        system.time(for (j in seq_len(calls)) fun())[["elapsed"]] / calls
    }, numeric(1)))
}
make_synthetic <- function(n) {
    GenomicRanges::GRanges(
        sample(paste0("chr", seq_len(22L)), n, replace = TRUE),
        IRanges::IRanges(sample.int(10000000L, n, replace = TRUE),
                         width = sample.int(10000L, n, replace = TRUE)),
        strand = sample(c("+", "-", "*"), n, replace = TRUE),
        feature_id = seq_len(n))
}
pairs <- function(x) sort(paste(queryHits(x), subjectHits(x)))
set.seed(100)
datasets <- list(airway = airway_genes,
                 synthetic_10000 = make_synthetic(10000L),
                 synthetic_100000 = make_synthetic(100000L))
results <- build_results <- list()
for (dataset in names(datasets)) {
    subject <- datasets[[dataset]]
    cat("Benchmarking", dataset, "with", length(subject), "subjects\n")
    build_results[[length(build_results) + 1L]] <- data.frame(
        dataset = dataset, subjects = length(subject),
        SuperIntervals_seconds = median_time(
            function() GenomicIntervalMap(subject), calls = 1L),
        GNCList_seconds = median_time(function() GNCList(subject), calls = 1L))
    index <- GenomicIntervalMap(subject)
    reference <- GNCList(subject)
    for (n_query in c(100L, 1000L, 10000L)) {
        query <- promoters(subject[sample(seq_along(subject), n_query,
                                          replace = TRUE)],
                           upstream = 1000L, downstream = 100L)
        for (ignore in c(FALSE, TRUE)) {
            counts <- countOverlaps(query, index, ignore.strand = ignore)
            reference_counts <- countOverlaps(query, reference,
                                             ignore.strand = ignore)
            stopifnot(identical(counts, reference_counts))
            actual <- findOverlaps(query, index, ignore.strand = ignore)
            expected <- findOverlaps(query, reference, ignore.strand = ignore)
            stopifnot(identical(pairs(actual), pairs(expected)))
            for (operation in c("countOverlaps", "findOverlaps")) {
                fun <- get(operation, envir = asNamespace("GenomicRanges"))
                methods <- list(SuperIntervals = index, GNCList = reference,
                                GRanges = subject)
                timings <- vapply(methods, function(indexed_subject) {
                    median_time(function() fun(query, indexed_subject,
                                              ignore.strand = ignore))
                }, numeric(1))
                results[[length(results) + 1L]] <- data.frame(
                    dataset = dataset, subjects = length(subject),
                    queries = n_query, ignore_strand = ignore,
                    operation = operation, hits = length(actual),
                    SuperIntervals_seconds = timings[["SuperIntervals"]],
                    GNCList_seconds = timings[["GNCList"]],
                    GRanges_seconds = timings[["GRanges"]])
            }
        }
    }
    rm(index, reference)
    invisible(gc())
}
results <- do.call(rbind, results)
results$GNCList_over_SuperIntervals <-
    results$GNCList_seconds / results$SuperIntervals_seconds
build_results <- do.call(rbind, build_results)
utils::write.csv(results, file.path(output_dir, "queries.csv"), row.names = FALSE)
utils::write.csv(build_results, file.path(output_dir, "build.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
writeLines(capture.output(Sys.info()), file.path(output_dir, "system-info.txt"))
print(results, row.names = FALSE)
print(build_results, row.names = FALSE)

# Scientific figure: elapsed times are measurements on this machine, not
# universal speed claims. Index construction is excluded from query timings.
grDevices::png(file.path(output_dir, "query-times.png"),
               width = 1800, height = 1400, res = 160)
old <- graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(1, 0, 2, 0))
colors <- c("#1864AB", "#D9480F", "#495057")
for (operation in c("countOverlaps", "findOverlaps")) {
    for (dataset in names(datasets)) {
        rows <- results[results$operation == operation &
                        results$dataset == dataset & !results$ignore_strand, ]
        times <- as.matrix(rows[c("SuperIntervals_seconds", "GNCList_seconds",
                                  "GRanges_seconds")]) * 1000
        graphics::matplot(rows$queries, pmax(times, 0.1), type = "b", log = "xy",
            col = colors, lty = 1L, pch = 16:18, xlab = "Query ranges",
            ylab = "Median elapsed time (ms)",
            main = paste(dataset, operation, sep = "\n"))
        graphics::legend("bottomright", c("SuperIntervals", "GNCList", "GRanges"),
                         col = colors, lty = 1L, pch = 16:18, bty = "n", cex = 0.7)
    }
}
graphics::mtext("GRanges interface: strand-aware repeated queries",
                outer = TRUE, font = 2L)
graphics::par(old)
grDevices::dev.off()
