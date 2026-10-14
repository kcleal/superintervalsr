# R wrapper maintenance assisted by OpenAI Codex (October 2026).
test_that("ingestion rejects invalid coordinates before integer conversion", {
    for (bad in list(NA_real_, Inf, 1.5, 2147483648)) {
        expect_error(IntervalMap(bad, ends = 10), "whole numbers")
        expect_error(add(IntervalMap(), bad, 10), "whole numbers")
    }
    expect_error(IntervalMap(5, ends = 1), "start")
    expect_error(IntervalMap(c(1, 2), ends = 3), "length")
    expect_error(add(IntervalMap(), c(1, 2), 10), "whole numbers")
    expect_error(reserve(IntervalMap(), -1), "non-negative")
    expect_error(flank(IntervalMap(), -1, 1), "non-negative")
    expect_error(gaps(IntervalMap(), 10, 1), "start")
})

test_that("coverage counts inclusive widths, including point intersections", {
    im <- IntervalMap(c(1, 5), ends = c(10, 15))
    expect_equal(coverage(im, 1, 15), list(count = 2L, total_coverage = 21))
    expect_equal(coverage(im, 10, 10), list(count = 2L, total_coverage = 2))
    expect_equal(coverage(im, 16, 20), list(count = 0L, total_coverage = 0))
    expect_equal(coverage(IntervalMap(), 1, 10),
                 list(count = 0L, total_coverage = 0))
    im <- IntervalMap(-2147483647, ends = 2147483647)
    expect_equal(coverage(im, -2147483647, 2147483647)$total_coverage,
                 4294967295)
})

test_that("set operations handle the full supported coordinate range", {
    hi <- .Machine$integer.max
    a <- IntervalMap(hi - 2L, ends = hi, values = "A")
    b <- IntervalMap(hi, ends = hi, values = "B")
    expect_equal(at(difference(a, b), 1)$end, hi - 1L)
    expect_length(gaps(a, hi - 2L, hi), 0L)
    expect_length(symmetric_difference(b, b), 0L)
    expect_equal(at(expand(a, hi, hi, lo = -hi, hi = hi), 1)$start, -2L)
    expect_length(expand(a, -hi, 0), 0L)
    low <- IntervalMap(-hi, ends = -hi)
    expect_equal(at(expand(low, hi, hi, lo = -hi, hi = hi), 1)$start, -hi)
    expect_length(flank(b, 0, hi), 0L)
    expect_equal(at(flank(a, hi, 0, lo = -hi), 1)$start, -2L)
})

test_that("result metadata survives source collection and later mutation", {
    make_result <- function() {
        value <- list(id = paste0("feature_", seq_len(50)), score = runif(50))
        src <- IntervalMap(c(10, 20), ends = c(15, 25),
                           values = list(value, value))
        result <- flank(src, 2, 2)
        clear(src)
        result
    }
    result <- make_result()
    invisible(gc())
    expect_equal(length(result), 4L)
    expect_equal(data_at(result, 1)$id, paste0("feature_", seq_len(50)))
    add(result, 30, 35, list(id = "new"))
    build(result)
    invisible(gc())
    expect_equal(data_at(result, 5)$id, "new")
    clear(result)
    invisible(gc())
    add(result, 1, 2, "reuse")
    build(result)
    expect_identical(search_values(result, 1, 2), list("reuse"))
})

test_that("restored external pointers fail with a useful error", {
    im <- IntervalMap(1, ends = 5)
    restored <- unserialize(serialize(im, NULL))
    expect_error(size(restored), "recreate maps")
})

test_that("single-sequence GRanges and source row mapping are explicit", {
    skip_if_not_installed("GenomicRanges")
    gr <- GenomicRanges::GRanges(c("chr1", "chr2", "chr1"),
                                IRanges::IRanges(c(20, 1, 1), c(30, 5, 10)),
                                id = c("C", "D", "A"))
    expect_error(IntervalMap(gr), "one sequence")
    im <- IntervalMap(gr, seqname = "chr1", value_column = "id")
    expect_equal(search_values(im, 1, 5), list("A"))
    hits <- search_items(im, 1, 5)
    out <- as_granges(hits, seqname = "chr1", valuename = "id")
    expect_equal(as.character(GenomicRanges::seqnames(out)), "chr1")
    expect_equal(unlist(S4Vectors::mcols(out)$id), "A")
    expect_error(IntervalMap(GenomicRanges::GRanges("chr1",
        IRanges::IRanges(start = 1, width = 0))), "start")
})

test_that("randomized closed queries agree with IRanges", {
    skip_if_not_installed("IRanges")
    set.seed(42)
    starts <- sample(-100:100, 200, replace = TRUE)
    ends <- starts + sample(0:40, 200, replace = TRUE)
    # Include duplicate, singleton and tied-start intervals.
    starts <- c(starts, 0, 0, 0)
    ends <- c(ends, 0, 0, 20)
    qs <- sample(-120:120, 100, replace = TRUE)
    qe <- qs + sample(0:20, 100, replace = TRUE)
    im <- IntervalMap(starts, ends = ends, values = seq_along(starts))
    reference <- IRanges::countOverlaps(IRanges::IRanges(qs, qe),
                                       IRanges::IRanges(starts, ends))
    expect_identical(count.batch(im, qs, qe), reference)
    expect_identical(as.integer(lengths(search_idxs.batch(im, qs, qe))), reference)
    for (i in seq_along(qs)) {
        expected <- which(starts <= qe[i] & ends >= qs[i])
        expect_setequal(as.integer(unlist(search_values(im, qs[i], qe[i]))), expected)
        widths <- pmin(ends[expected], qe[i]) -
                  pmax(starts[expected], qs[i]) + 1
        expect_equal(coverage(im, qs[i], qe[i])$total_coverage, sum(widths))
    }
})

test_that("adjacency is preserved and set results are indexed", {
    a <- IntervalMap(c(1, 6), ends = c(5, 10), values = c("A", "B"))
    expect_equal(length(merge_overlaps(a)), 2L)
    b <- IntervalMap(c(3, 8), ends = c(4, 12))
    expect_equal(count(intersection(a, b), 1, 20), 2L)
    expect_equal(count(difference(a, b), 1, 20), 3L)
    expect_equal(count(unique_intervals(a), 1, 20), 2L)
})

test_that("metadata references are released when result maps are collected", {
    state <- new.env(parent = emptyenv())
    state$n <- 0L
    finalized <- function(e) state$n <- state$n + 1L
    exercise <- function() {
        value <- new.env(parent = emptyenv())
        reg.finalizer(value, finalized)
        src <- IntervalMap(10, ends = 20, values = list(value))
        result <- flank(src, 2, 2) # two occurrences of the same reference
        clear(src)
        add(result, 30, 40, value) # a third occurrence
        build(result)
        invisible(NULL)
    }
    for (i in seq_len(10)) exercise()
    for (i in seq_len(3)) invisible(gc())
    expect_identical(state$n, 10L)
})

test_that("R-owned metadata survives incremental growth, sorting and collection", {
    im <- IntervalMap()
    for (i in seq_len(20L)) {
        add(im, 21L - i, 30L - i, list(id = i, text = paste0("value_", i)))
        invisible(gc())
    }
    build(im)
    expect_identical(vapply(seq_len(20L), function(i) data_at(im, i)$id,
        integer(1)), 20:1)
    result <- difference(im, IntervalMap(5L, ends = 10L))
    clear(im)
    invisible(gc())
    expect_true(length(result) > 0L)
    expect_true(is.list(data_at(result, 1L)))
})
