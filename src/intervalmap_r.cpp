// R wrapper maintenance assisted by OpenAI Codex (October 2026).
#include <Rcpp.h>
#include "../inst/include/superintervals.hpp"

using namespace Rcpp;
using namespace si;

// Each stored occurrence owns one preserved reference, including result maps.
inline void preserve_sexp(SEXP value) {
    if (value != R_NilValue) R_PreserveObject(value);
}
inline void release_sexp(SEXP value) {
    if (value != R_NilValue) R_ReleaseObject(value);
}
void cleanup_intervalmap(IntervalMap<int, SEXP>* ptr) {
    if (!ptr) return;
    for (SEXP value : ptr->data) release_sexp(value);
    delete ptr;
}
using MapPtr = XPtr<IntervalMap<int, SEXP>, PreserveStorage,
                    cleanup_intervalmap>;

MapPtr checked_map(SEXP container) {
    if (TYPEOF(container) != EXTPTRSXP || !R_ExternalPtrAddr(container)) {
        stop("Invalid IntervalMap pointer; recreate maps after loading saved R objects");
    }
    return MapPtr(container);
}
void preserve_values(IntervalMap<int, SEXP>* ptr) {
    for (SEXP value : ptr->data) preserve_sexp(value);
}
SEXP wrap_new_intervalmap(IntervalMap<int, SEXP>* ptr) {
    MapPtr main_ptr(ptr, true);
    XPtr<std::vector<SEXP>> value_cache_ptr(new std::vector<SEXP>(), true);
    main_ptr.attr("value_cache") = value_cache_ptr;
    XPtr<std::vector<size_t>> idx_cache_ptr(new std::vector<size_t>(), true);
    main_ptr.attr("idx_cache") = idx_cache_ptr;
    return main_ptr;
}
// Widen arithmetic for operations that use end + 1 or signed resizing.
// Query maps retain 32-bit coordinates and their SIMD counting implementation.
IntervalMap<int64_t, SEXP> wide_map(const IntervalMap<int, SEXP>& src, bool indexed = false) {
    IntervalMap<int64_t, SEXP> out;
    out.reserve(src.starts.size());
    for (size_t k = 0; k < src.starts.size(); ++k)
        out.add(src.starts[k], src.ends[k], src.data[k]);
    if (indexed) out.build();
    return out;
}
SEXP wrap_wide_map(const IntervalMap<int64_t, SEXP>& src) {
    auto out = new IntervalMap<int, SEXP>();
    out->reserve(src.starts.size());
    for (size_t k = 0; k < src.starts.size(); ++k)
        out->add(static_cast<int>(src.starts[k]),
                 static_cast<int>(src.ends[k]), src.data[k]);
    preserve_values(out);
    return wrap_new_intervalmap(out);
}

// [[Rcpp::export]]
SEXP create_intervalmap() {
    Shield<SEXP> out(wrap_new_intervalmap(new IntervalMap<int, SEXP>()));
    return out;
}

// [[Rcpp::export]]
SEXP create_intervalmap_from_vectors(IntegerVector starts, IntegerVector ends, List values) {
    if (starts.size() != ends.size()) {
        stop("starts and ends must have the same length");
    }

    if (values.size() != 0 && values.size() != starts.size()) {
        stop("values must have the same length as starts and ends");
    }

    IntervalMap<int, SEXP>* ptr = new IntervalMap<int, SEXP>();
    ptr->reserve(starts.size());
    for (int i = 0; i < starts.size(); ++i) {
        SEXP value = (values.size() > 0) ? values[i] : R_NilValue;
        preserve_sexp(value);
        ptr->add(starts[i], ends[i], value);
    }
    ptr->build();  // Auto-build for convenience

    Shield<SEXP> out(wrap_new_intervalmap(ptr));
    return out;
}

// [[Rcpp::export]]
void add_interval(SEXP container, int start, int end, SEXP value = R_NilValue) {
    auto si = checked_map(container);

    // Preserve the SEXP reference
    preserve_sexp(value);

    // Add to the interval map - this is now a simple operation
    si->add(start, end, value);
}

// [[Rcpp::export]]
void build_index(SEXP container) {
    auto si = checked_map(container);
    si->build();
}

// [[Rcpp::export]]
List get_interval_at(SEXP container, int r_index) {
    auto si = checked_map(container);

    // Convert R's 1-based to C++'s 0-based
    int cpp_index = r_index - 1;

    if (si->size() == 0 || cpp_index < 0 || cpp_index >= (int)si->size()) {
        stop("Index out of range");
    }

    int start = si->starts[cpp_index];
    int end = si->ends[cpp_index];
    SEXP value = si->data[cpp_index];

    return List::create(
        Named("start") = start,
        Named("end") = end,
        Named("value") = value
    );
}

// [[Rcpp::export]]
int get_start_at(SEXP container, int r_index) {
    auto si = checked_map(container);

    // Convert R's 1-based to C++'s 0-based
    int cpp_index = r_index - 1;

    if (si->size() == 0 || cpp_index < 0 || cpp_index >= (int)si->size()) {
        stop("Index out of range");
    }

    return si->starts[cpp_index];
}

// [[Rcpp::export]]
int get_end_at(SEXP container, int r_index) {
    auto si = checked_map(container);

    // Convert R's 1-based to C++'s 0-based
    int cpp_index = r_index - 1;

    if (si->size() == 0 || cpp_index < 0 || cpp_index >= (int)si->size()) {
        stop("Index out of range");
    }

    return si->ends[cpp_index];
}

// [[Rcpp::export]]
SEXP get_data_at(SEXP container, int r_index) {
    auto si = checked_map(container);

    // Convert R's 1-based to C++'s 0-based
    int cpp_index = r_index - 1;

    if (si->size() == 0 || cpp_index < 0 || cpp_index >= (int)si->size()) {
        stop("Index out of range");
    }

    return si->data[cpp_index];
}

// [[Rcpp::export]]
IntegerVector get_all_starts(SEXP container) {
    auto si = checked_map(container);
    IntegerVector result(si->size());

    for (size_t i = 0; i < si->size(); ++i) {
        result[i] = si->starts[i];
    }

    return result;
}

// [[Rcpp::export]]
IntegerVector get_all_ends(SEXP container) {
    auto si = checked_map(container);
    IntegerVector result(si->size());

    for (size_t i = 0; i < si->size(); ++i) {
        result[i] = si->ends[i];
    }

    return result;
}

// [[Rcpp::export]]
List get_all_values(SEXP container) {
    auto si = checked_map(container);
    List result(si->size());

    for (size_t i = 0; i < si->size(); ++i) {
        result[i] = si->data[i];
    }

    return result;
}

// [[Rcpp::export]]
void clear_intervals(SEXP container) {
    auto si = checked_map(container);

    // Release all SEXP references before clearing
    for (const SEXP& value : si->data) {
        release_sexp(value);
    }

    si->clear();
    si->start_sorted = si->end_sorted = true;
}

// [[Rcpp::export]]
void reserve_intervals(SEXP container, int n) {
    auto si = checked_map(container);
    si->reserve(n);
}

// [[Rcpp::export]]
int get_size(SEXP container) {
    auto si = checked_map(container);
    return si->size();
}

// [[Rcpp::export]]
bool cpp_has_overlaps(SEXP container, int start, int end) {
    auto si = checked_map(container);
    return si->has_overlaps(start, end);
}

// [[Rcpp::export]]
LogicalVector has_overlaps_batch(SEXP container, IntegerVector starts, IntegerVector ends) {
    if (starts.size() != ends.size()) {
        stop("starts and ends must have the same length");
    }

    auto si = checked_map(container);
    LogicalVector results(starts.size());

    for (int i = 0; i < starts.size(); ++i) {
        results[i] = si->has_overlaps(starts[i], ends[i]);
    }

    return results;
}

// [[Rcpp::export]]
int count_overlaps(SEXP container, int start, int end) {
    auto si = checked_map(container);
    return si->count(start, end);
}

// [[Rcpp::export]]
List cpp_search_values(SEXP container, int start, int end) {
    auto si = checked_map(container);
    XPtr<std::vector<SEXP>> cache_ptr(si.attr("value_cache"));
    std::vector<SEXP>& found_values = *cache_ptr;
    found_values.clear();

    si->search_values(start, end, found_values);

    List result(found_values.size());
    for (size_t i = 0; i < found_values.size(); ++i) {
        result[i] = found_values[i];
    }

    return result;
}

// [[Rcpp::export]]
IntegerVector search_indexes(SEXP container, int start, int end) {
    auto si = checked_map(container);
    XPtr<std::vector<size_t>> cache_ptr(si.attr("idx_cache"));
    std::vector<size_t>& found_indices = *cache_ptr;
    found_indices.clear();

    si->search_idxs(start, end, found_indices);

    IntegerVector result(found_indices.size());
    for (size_t i = 0; i < found_indices.size(); ++i) {
        result[i] = found_indices[i] + 1;  // Convert to 1-based indexing for R
    }

    return result;
}

// [[Rcpp::export]]
List cpp_search_keys(SEXP container, int start, int end) {
    auto si = checked_map(container);
    XPtr<std::vector<size_t>> cache_ptr(si.attr("idx_cache"));
    std::vector<size_t>& found_indices = *cache_ptr;
    found_indices.clear();

    si->search_idxs(start, end, found_indices);

    if (found_indices.empty()) {
        return List::create(
            Named("start") = IntegerVector(),
            Named("end") = IntegerVector()
        );
    }

    IntegerVector starts_vec(found_indices.size());
    IntegerVector ends_vec(found_indices.size());

    for (size_t i = 0; i < found_indices.size(); ++i) {
        size_t idx = found_indices[i];
        starts_vec[i] = si->starts[idx];
        ends_vec[i] = si->ends[idx];
    }

    return List::create(
        Named("start") = starts_vec,
        Named("end") = ends_vec
    );
}

// [[Rcpp::export]]
List cpp_search_items(SEXP container, int start, int end) {
    auto si = checked_map(container);
    XPtr<std::vector<size_t>> cache_ptr(si.attr("idx_cache"));
    std::vector<size_t>& found_indices = *cache_ptr;
    found_indices.clear();

    si->search_idxs(start, end, found_indices);

    if (found_indices.empty()) {
        return List::create(
            Named("start") = IntegerVector(),
            Named("end") = IntegerVector(),
            Named("value") = List()
        );
    }

    IntegerVector starts_vec(found_indices.size());
    IntegerVector ends_vec(found_indices.size());
    List values_list(found_indices.size());

    for (size_t i = 0; i < found_indices.size(); ++i) {
        size_t idx = found_indices[i];
        starts_vec[i] = si->starts[idx];
        ends_vec[i] = si->ends[idx];
        values_list[i] = si->data[idx];
    }

    return List::create(
        Named("start") = starts_vec,
        Named("end") = ends_vec,
        Named("value") = values_list
    );
}

// [[Rcpp::export]]
List get_coverage(SEXP container, int start, int end) {
    auto si = checked_map(container);

    int count = 0;
    double total = 0;
    for (size_t k : si->search_idxs(start, end)) {
        ++count;
        total += static_cast<double>(std::min(si->ends[k], end)) -
                 static_cast<double>(std::max(si->starts[k], start)) + 1;
    }
    return List::create(Named("count") = count,
                        Named("total_coverage") = total);
}

// [[Rcpp::export]]
List search_idxs_batch(SEXP container, IntegerVector starts, IntegerVector ends) {
    if (starts.size() != ends.size()) {
        stop("starts and ends must have the same length");
    }

    auto si = checked_map(container);
    XPtr<std::vector<size_t>> cache_ptr(si.attr("idx_cache"));
    std::vector<size_t>& found_indices = *cache_ptr;

    List results(starts.size());

    for (int i = 0; i < starts.size(); ++i) {
        found_indices.clear();
        si->search_idxs(starts[i], ends[i], found_indices);

        IntegerVector result(found_indices.size());
        for (size_t j = 0; j < found_indices.size(); ++j) {
            result[j] = found_indices[j] + 1;  // Convert to 1-based
        }
        results[i] = result;
    }

    return results;
}

// [[Rcpp::export]]
IntegerVector count_batch(SEXP container, IntegerVector starts, IntegerVector ends) {
    if (starts.size() != ends.size()) {
        stop("starts and ends must have the same length");
    }

    auto si = checked_map(container);
    IntegerVector results(starts.size());

    for (int i = 0; i < starts.size(); ++i) {
        results[i] = si->count(starts[i], ends[i]);
    }

    return results;
}

// [[Rcpp::export]]
List search_values_batch(SEXP container, IntegerVector starts, IntegerVector ends) {
    if (starts.size() != ends.size()) {
        stop("starts and ends must have the same length");
    }

    auto si = checked_map(container);
    XPtr<std::vector<SEXP>> cache_ptr(si.attr("value_cache"));
    std::vector<SEXP>& found_values = *cache_ptr;

    List results(starts.size());

    for (int i = 0; i < starts.size(); ++i) {
        found_values.clear();
        si->search_values(starts[i], ends[i], found_values);

        List result(found_values.size());
        for (size_t j = 0; j < found_values.size(); ++j) {
            result[j] = found_values[j];
        }
        results[i] = result;
    }

    return results;
}
// [[Rcpp::export]]
List cpp_search_point(SEXP container, int point) {
    auto si = checked_map(container);
    XPtr<std::vector<SEXP>> cache_ptr(si.attr("value_cache"));
    std::vector<SEXP>& found_values = *cache_ptr;
    found_values.clear();

    si->search_point(point, found_values);

    List result(found_values.size());
    for (size_t i = 0; i < found_values.size(); ++i) {
        result[i] = found_values[i];
    }

    return result;
}

// [[Rcpp::export]]
SEXP cpp_merge_overlaps(SEXP container) {
    auto si = checked_map(container);
    IntervalMap<int, SEXP>* out = new IntervalMap<int, SEXP>(si->merge_overlaps());
    preserve_values(out);
    return wrap_new_intervalmap(out);
}

// [[Rcpp::export]]
SEXP cpp_union_with(SEXP container_a, SEXP container_b) {
    auto si_a = checked_map(container_a);
    auto si_b = checked_map(container_b);
    IntervalMap<int, SEXP>* out = new IntervalMap<int, SEXP>(si_a->union_with(*si_b));
    preserve_values(out);
    return wrap_new_intervalmap(out);
}

// [[Rcpp::export]]
SEXP cpp_intersection(SEXP container_a, SEXP container_b) {
    auto si_a = checked_map(container_a);
    auto si_b = checked_map(container_b);
    IntervalMap<int, SEXP>* out = new IntervalMap<int, SEXP>(si_a->intersection(*si_b));
    preserve_values(out);
    return wrap_new_intervalmap(out);
}

// [[Rcpp::export]]
SEXP cpp_difference(SEXP container_a, SEXP container_b) {
    auto si_a = checked_map(container_a);
    auto si_b = checked_map(container_b);
    return wrap_wide_map(wide_map(*si_a).difference(wide_map(*si_b, true)));
}

// [[Rcpp::export]]
SEXP cpp_symmetric_difference(SEXP container_a, SEXP container_b) {
    auto si_a = checked_map(container_a);
    auto si_b = checked_map(container_b);
    return wrap_wide_map(wide_map(*si_a, true).symmetric_difference(wide_map(*si_b, true)));
}

// [[Rcpp::export]]
SEXP cpp_gaps(SEXP container, int lo, int hi, SEXP fill = R_NilValue) {
    auto si = checked_map(container);
    return wrap_wide_map(wide_map(*si).gaps(lo, hi, fill));
}

// [[Rcpp::export]]
SEXP cpp_span(SEXP container) {
    auto si = checked_map(container);

    std::pair<int, int> result;
    bool ok = si->span(result);

    if (!ok) {
        return R_NilValue;
    }

    return List::create(
        Named("start") = result.first,
        Named("end") = result.second
    );
}

// [[Rcpp::export]]
SEXP cpp_expand(SEXP container, int left, int right, int lo, int hi) {
    auto si = checked_map(container);
    return wrap_wide_map(wide_map(*si).expand(left, right, lo, hi));
}

// [[Rcpp::export]]
SEXP cpp_flank(SEXP container, int left, int right, int lo, int hi) {
    auto si = checked_map(container);
    return wrap_wide_map(wide_map(*si).flank(left, right, lo, hi));
}

// [[Rcpp::export]]
SEXP cpp_unique_intervals(SEXP container) {
    auto si = checked_map(container);
    IntervalMap<int, SEXP>* out = new IntervalMap<int, SEXP>(si->unique());
    preserve_values(out);
    return wrap_new_intervalmap(out);
}
