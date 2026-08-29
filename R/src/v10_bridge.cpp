#include <Rcpp.h>
#include "../../C++/straw_v10.h"

bool strawIsV10(const std::string &path) { return straw_v10::isV10(path); }

Rcpp::DataFrame strawReadV10(const std::string &norm, const std::string &path,
                             const std::string &chr1, const std::string &chr2,
                             const std::string &unit, int32_t resolution,
                             const std::string &matrix) {
    std::vector<int32_t> x, y;
    std::vector<double> counts;
    straw_v10::File(path).stream(matrix, norm, chr1, chr2, unit, resolution,
        [&](const contactRecord &record) {
            x.push_back(record.binX);
            y.push_back(record.binY);
            counts.push_back(record.counts);
        });
    return Rcpp::DataFrame::create(Rcpp::Named("x") = x, Rcpp::Named("y") = y,
                                   Rcpp::Named("counts") = counts);
}

Rcpp::DataFrame strawChromosomesV10(const std::string &path) {
    Rcpp::IntegerVector index;
    Rcpp::CharacterVector name;
    Rcpp::NumericVector length;
    for (const auto &c : straw_v10::File(path).chromosomes()) {
        index.push_back(c.index); name.push_back(c.name); length.push_back(c.length);
    }
    return Rcpp::DataFrame::create(Rcpp::Named("index") = index, Rcpp::Named("name") = name,
                                   Rcpp::Named("length") = length);
}

Rcpp::NumericVector strawResolutionsV10(const std::string &path) {
    return Rcpp::wrap(straw_v10::File(path).resolutions());
}

Rcpp::CharacterVector strawNormalizationsV10(const std::string &path) {
    auto values = straw_v10::File(path).normalizations();
    values.push_back("NONE");
    return Rcpp::wrap(values);
}
