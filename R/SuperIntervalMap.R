# Genomic adapter and tests assisted by OpenAI Codex (October 2026).

#' An indexed collection of genomic ranges
#'
#' Index a GRanges once for repeated chromosome- and strand-aware overlap
#' queries.
#' The original ranges, sequence information, names and metadata are retained.
#'
#' @param x,object An indexed genomic map.
#' @slot ranges The original GRanges, in input order.
#' @slot index Internal per-sequence and per-strand indexes. Do not edit slots;
#'   construct a new object when the subject ranges change.
#' @return The S4 class stores indexed GRanges. Its length method returns an
#'   integer feature count; show prints a summary and returns invisible NULL.
#' @importFrom methods show
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps countOverlaps
#' @exportClass SuperIntervalMap
#' @name SuperIntervalMap-class
#' @seealso [SuperIntervalMap()]
methods::setClass("SuperIntervalMap",
    slots = c(ranges = "GRanges", index = "list")
)

#' Create an indexed genomic subject
#'
#' @param x A GRanges containing subject features, possibly on multiple
#' sequences
#'   and strands. Coordinates, names, seqinfo and metadata are preserved.
#' @return A SuperIntervalMap S4 object for use as the subject of
#'   `findOverlaps()` and `countOverlaps()`. `as_granges()` retrieves the
#' original
#'   GRanges, and `length()` returns its number of features.
#' @details Hits refer to the original subject row order. Strand matching
#' follows
#'   GenomicRanges: `"*"` matches either strand, and `ignore.strand = TRUE`
#'   ignores strand entirely. Ordinary positive-width overlaps use
#' SuperIntervals.
#'   Queries involving zero-width ranges, a non-default gap or overlap
#' threshold,
#'   or a type other than `"any"` delegate to GenomicRanges to preserve its
#'   semantics. The object stores up to three indexes per sequence, trading
#'   memory for repeated-query speed. Query functions assume the indexes are
#'   current; do not alter the slots after construction.
#'
#'   Like IntervalMap, the indexes contain external pointers. Save
#'   `as_granges(x)` and reconstruct the index after loading serialized data.
#' @export
#' @examples
#' library(GenomicRanges)
#' genes <- GRanges(c("chr1", "chr2", "chr1"),
#'     IRanges(c(20, 1, 1), c(30, 10, 10)),
#'     strand = c("+", "-", "*"), gene_id = c("C", "B", "A"))
#' index <- SuperIntervalMap(genes)
#' query <- GRanges(c("chr1", "chr2"), IRanges(c(5, 5), c(25, 8)),
#'     strand = c("+", "+"))
#' findOverlaps(query, index)
#' countOverlaps(query, index, ignore.strand = TRUE)
#' as_granges(index)
SuperIntervalMap <- function(x) {
    if (!methods::is(x, "GRanges")) {
        stop("x must be a GRanges", call. = FALSE)
    }
    # Zero-width features are kept for the GenomicRanges fallback, but cannot
    # enter the closed-interval index.
    rows <- which(GenomicRanges::width(x) > 0L)
    by_sequence <- split(rows, as.character(GenomicRanges::seqnames(x))[rows])
    strands <- as.character(GenomicRanges::strand(x))
    starts <- GenomicRanges::start(x)
    ends <- GenomicRanges::end(x)
    make_index <- function(ids) {
        map <- IntervalMap(starts[ids], ends = ends[ids], values = ids)
        list(map = map, rows = as.integer(unlist(get_all_values(map))))
    }
    index <- lapply(by_sequence, function(ids) {
        list(all = make_index(ids),
                plus = make_index(ids[strands[ids] != "-"]),
                minus = make_index(ids[strands[ids] != "+"]))
    })
    attr(index, "zero.width") <- any(GenomicRanges::width(x) == 0L)
    methods::new("SuperIntervalMap", ranges = x, index = index)
}

# Internal batch dispatch. No per-interval index-state checks are performed.
.genomic_query <- function(query, subject, ignore.strand, count.only = FALSE) {
    if (!is.logical(ignore.strand) || length(ignore.strand) != 1L ||
        is.na(ignore.strand)) {
        stop("ignore.strand must be TRUE or FALSE", call. = FALSE)
    }
    by_sequence <- split(seq_along(query),
                            as.character(GenomicRanges::seqnames(query)))
    strands <- as.character(GenomicRanges::strand(query))
    starts <- GenomicRanges::start(query)
    ends <- GenomicRanges::end(query)
    counts <- integer(length(query))
    names(counts) <- names(query)
    query_hits <- subject_hits <- list()
    for (seqname in names(by_sequence)) {
        indexes <- subject@index[[seqname]]
        if (is.null(indexes)) next
        ids <- by_sequence[[seqname]]
        groups <- if (ignore.strand) list(all = ids) else {
            list(plus = ids[strands[ids] == "+"],
                    minus = ids[strands[ids] == "-"],
                    all = ids[strands[ids] == "*"])
        }
        for (group in names(groups)) {
            qi <- groups[[group]]
            if (!length(qi)) next
            index <- indexes[[group]]
            if (count.only) {
                counts[qi] <- count_batch(index$map, starts[qi], ends[qi])
            } else {
                found <- search_idxs_batch(index$map, starts[qi], ends[qi])
                query_hits[[length(query_hits) + 1L]] <- rep(qi, lengths(found))
                subject_hits[[length(subject_hits) + 1L]] <-
                    index$rows[unlist(found, use.names = FALSE)]
            }
        }
    }
    if (count.only) return(counts)
    S4Vectors::Hits(as.integer(unlist(query_hits, use.names = FALSE)),
                    as.integer(unlist(subject_hits, use.names = FALSE)),
                    length(query), length(subject@ranges), sort.by.query = TRUE)
}

.genomic_fast_path <- function(query, subject, maxgap, minoverlap, type) {
    identical(maxgap, -1L) && length(minoverlap) == 1L &&
        minoverlap %in% c(0L, 1L) && type == "any" &&
        all(GenomicRanges::width(query) > 0L) &&
        !isTRUE(attr(subject@index, "zero.width"))
}

#' Query an indexed genomic subject
#'
#' @param query A GRanges of query features.
#' @param subject A SuperIntervalMap constructed once from the subject
#' GRanges.
#' @param maxgap,minoverlap,type,select Overlap options, as in
#'   [GenomicRanges::findOverlaps()].
#' @param ignore.strand Whether to ignore strand (default FALSE).
#' @return `findOverlaps()` returns a Hits object with indices in the original
#'   input order, or a subject-index vector when `select` is not `"all"`.
#'   `countOverlaps()` returns one integer count per query feature.
#' @rdname genomic-overlaps
#' @aliases findOverlaps,GRanges,SuperIntervalMap-method
#' @export
#' @examples
#' library(GenomicRanges)
#' subject <- GRanges(c("chr1", "chr2"), IRanges(c(1, 10), c(5, 20)))
#' query <- GRanges(c("chr1", "chr2"), IRanges(c(3, 15), c(8, 18)))
#' index <- SuperIntervalMap(subject)
#' hits <- findOverlaps(query, index)
#' as_granges(index)[subjectHits(hits)]
#' countOverlaps(query, index)
methods::setMethod("findOverlaps",
    signature(query = "GRanges", subject = "SuperIntervalMap"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
                type = c("any", "start", "end", "within", "equal"),
                select = c("all", "first", "last", "arbitrary"),
                ignore.strand = FALSE) {
        type <- match.arg(type)
        select <- match.arg(select)
        if (!.genomic_fast_path(query, subject, maxgap, minoverlap, type)) {
            return(GenomicRanges::findOverlaps(query, subject@ranges,
                maxgap = maxgap, minoverlap = minoverlap, type = type,
                select = select, ignore.strand = ignore.strand))
        }
        hits <- .genomic_query(query, subject, ignore.strand)
        S4Vectors::selectHits(hits, select = select)
    }
)

#' @rdname genomic-overlaps
#' @aliases countOverlaps,GRanges,SuperIntervalMap-method
#' @export
methods::setMethod("countOverlaps",
    signature(query = "GRanges", subject = "SuperIntervalMap"),
    function(query, subject, maxgap = -1L, minoverlap = 0L,
                type = c("any", "start", "end", "within", "equal"),
                ignore.strand = FALSE) {
        type <- match.arg(type)
        if (!.genomic_fast_path(query, subject, maxgap, minoverlap, type)) {
            return(GenomicRanges::countOverlaps(query, subject@ranges,
                maxgap = maxgap, minoverlap = minoverlap, type = type,
                ignore.strand = ignore.strand))
        }
        .genomic_query(query, subject, ignore.strand, count.only = TRUE)
    }
)

#' @rdname SuperIntervalMap-class
#' @export
methods::setMethod("length", "SuperIntervalMap", function(x) {
    length(x@ranges)
})

#' @rdname SuperIntervalMap-class
#' @export
methods::setMethod("show", "SuperIntervalMap", function(object) {
    cat("SuperIntervalMap with", length(object), "ranges on",
        length(object@index), "indexed sequences\n")
})
