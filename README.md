# superintervalsr

`superintervalsr` is an R wrapper for the SuperIntervals C++ library, designed
for repeated overlap queries on an indexed collection of intervals. It is being
prepared for Bioconductor submission and is not yet available there.

## Installation

Install dependencies from CRAN and Bioconductor, then install the source package:

```r
install.packages(c("Rcpp", "remotes"))
install.packages("BiocManager")
BiocManager::install(c("S4Vectors", "IRanges", "GenomicRanges"))
remotes::install_github("kcleal/superintervalsr")
```

A C++17 compiler is required. See `BIOCONDUCTOR.md` for development checks and
submission steps.

## Quick start

```r
library(superintervalsr)
im <- IntervalMap(c(20, 1, 5), ends = c(30, 10, 15),
                  values = c("C", "A", "B"))
count.batch(im, c(7, 25), c(12, 35))
search_items(im, 7, 12)
search_point(im, 10)
coverage(im, 1, 15) # count = 2, total_coverage = 21 inclusive bases

# Incremental loading: build once after adding all intervals.
im <- IntervalMap()
reserve(im, 1000)
add(im, 1, 10, "A")
add(im, 5, 15, "B")
build(im)
search_values(im, 7, 12)
```

Intervals and queries are **closed**: both endpoints are included. Inputs must
be finite whole numbers in the R integer range (excluding `NA`). Constructors
validate coordinates and build automatically. Query functions retain a thin,
low-overhead interface: callers must provide valid query ranges and call
`build()` after every batch of additions. Index state is not checked per query.

`build()` sorts starts ascending and ends descending within tied starts. Search
indices address this internal order, not original input rows. Hit order and the
choice of metadata among identical coordinate ties are unspecified. For mapping
back to source rows, store row IDs as values.

## Genomic ranges

Each map represents one sequence and ignores strand. Filter explicitly when
input contains multiple chromosomes; construct and query separate maps for each.

```r
library(GenomicRanges)
gr <- GRanges(c("chr1", "chr1", "chr2"),
              IRanges(c(20, 1, 1), c(30, 10, 10)),
              gene_id = c("C", "A", "D"))
im <- IntervalMap(gr, value_column = "gene_id", seqname = "chr1")
hits <- superintervalsr::search_items(im, 5, 25)
superintervalsr::as_granges(hits, seqname = "chr1", valuename = "gene_id")
```

## Set operations

`merge_overlaps()`, `union_with()`, `intersection()`, `difference()`,
`symmetric_difference()`, `gaps()`, `expand()`, `flank()` and
`unique_intervals()` return independent, already indexed maps. `span()` returns
the enclosing start/end pair or `NULL` for an empty map.

Merge and union coalesce overlaps, but keep adjacent non-overlapping intervals
separate. Intersection returns one piece per overlapping pair, retaining the
left input's value; use `merge_overlaps()` for disjoint output. Difference
retains the left input's metadata, while expansion and flanks inherit each
source interval's value. These operations do not reconcile all metadata from
both inputs. Consult the help pages for exact value semantics.

`coverage()` returns the number of hits and the sum of inclusive intersection
widths. Shared bases are counted once per stored interval, so this is not union
width or the per-base coverage returned by `IRanges::coverage()`.

Maps mutate by reference: `other <- im` aliases the same map. External-pointer
maps cannot be restored with `readRDS()`; save coordinates and values and
reconstruct them. Qualify names such as `superintervalsr::coverage()` and
`superintervalsr::flank()` when working alongside Bioconductor packages.
