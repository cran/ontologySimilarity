// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sample_null
NumericVector sample_null(std::string type, RObject term_sets_data, bool use_mean, int group_size, int samples);
RcppExport SEXP _ontologySimilarity_sample_null(SEXP typeSEXP, SEXP term_sets_dataSEXP, SEXP use_meanSEXP, SEXP group_sizeSEXP, SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< RObject >::type term_sets_data(term_sets_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type use_mean(use_meanSEXP);
    Rcpp::traits::input_parameter< int >::type group_size(group_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_null(type, term_sets_data, use_mean, group_size, samples));
    return rcpp_result_gen;
END_RCPP
}
// sim_p
double sim_p(std::string type, RObject term_sets_data, bool use_mean, IntegerVector group, int min_its, int max_its, double signif, double dismiss);
RcppExport SEXP _ontologySimilarity_sim_p(SEXP typeSEXP, SEXP term_sets_dataSEXP, SEXP use_meanSEXP, SEXP groupSEXP, SEXP min_itsSEXP, SEXP max_itsSEXP, SEXP signifSEXP, SEXP dismissSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< RObject >::type term_sets_data(term_sets_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type use_mean(use_meanSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type min_its(min_itsSEXP);
    Rcpp::traits::input_parameter< int >::type max_its(max_itsSEXP);
    Rcpp::traits::input_parameter< double >::type signif(signifSEXP);
    Rcpp::traits::input_parameter< double >::type dismiss(dismissSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_p(type, term_sets_data, use_mean, group, min_its, max_its, signif, dismiss));
    return rcpp_result_gen;
END_RCPP
}
// group_sim
double group_sim(std::string type, RObject term_sets_data, bool use_mean, IntegerVector group);
RcppExport SEXP _ontologySimilarity_group_sim(SEXP typeSEXP, SEXP term_sets_dataSEXP, SEXP use_meanSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< RObject >::type term_sets_data(term_sets_dataSEXP);
    Rcpp::traits::input_parameter< bool >::type use_mean(use_meanSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(group_sim(type, term_sets_data, use_mean, group));
    return rcpp_result_gen;
END_RCPP
}
// sim_grid
NumericMatrix sim_grid(IntegerVector term_ids1, IntegerVector case_ids1, int num_cases1, IntegerVector term_ids2, IntegerVector case_ids2, int num_cases2, NumericMatrix ttsm);
RcppExport SEXP _ontologySimilarity_sim_grid(SEXP term_ids1SEXP, SEXP case_ids1SEXP, SEXP num_cases1SEXP, SEXP term_ids2SEXP, SEXP case_ids2SEXP, SEXP num_cases2SEXP, SEXP ttsmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type term_ids1(term_ids1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type case_ids1(case_ids1SEXP);
    Rcpp::traits::input_parameter< int >::type num_cases1(num_cases1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type term_ids2(term_ids2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type case_ids2(case_ids2SEXP);
    Rcpp::traits::input_parameter< int >::type num_cases2(num_cases2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ttsm(ttsmSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_grid(term_ids1, case_ids1, num_cases1, term_ids2, case_ids2, num_cases2, ttsm));
    return rcpp_result_gen;
END_RCPP
}
// calc_term_sim_mat
NumericMatrix calc_term_sim_mat(IntegerVector anc_start, IntegerVector anc_stop, IntegerVector ancestors, NumericVector info, IntegerVector terms1, IntegerVector terms2);
RcppExport SEXP _ontologySimilarity_calc_term_sim_mat(SEXP anc_startSEXP, SEXP anc_stopSEXP, SEXP ancestorsSEXP, SEXP infoSEXP, SEXP terms1SEXP, SEXP terms2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type anc_start(anc_startSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type anc_stop(anc_stopSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ancestors(ancestorsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type info(infoSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type terms1(terms1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type terms2(terms2SEXP);
    rcpp_result_gen = Rcpp::wrap(calc_term_sim_mat(anc_start, anc_stop, ancestors, info, terms1, terms2));
    return rcpp_result_gen;
END_RCPP
}
// sim_grid_ic
NumericMatrix sim_grid_ic(bool lin, IntegerVector anc_start, IntegerVector anc_stop, IntegerVector ancestors, NumericVector info, IntegerVector term_ids1, IntegerVector case_ids1, int num_cases1, IntegerVector term_ids2, IntegerVector case_ids2, int num_cases2);
RcppExport SEXP _ontologySimilarity_sim_grid_ic(SEXP linSEXP, SEXP anc_startSEXP, SEXP anc_stopSEXP, SEXP ancestorsSEXP, SEXP infoSEXP, SEXP term_ids1SEXP, SEXP case_ids1SEXP, SEXP num_cases1SEXP, SEXP term_ids2SEXP, SEXP case_ids2SEXP, SEXP num_cases2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type lin(linSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type anc_start(anc_startSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type anc_stop(anc_stopSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ancestors(ancestorsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type info(infoSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type term_ids1(term_ids1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type case_ids1(case_ids1SEXP);
    Rcpp::traits::input_parameter< int >::type num_cases1(num_cases1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type term_ids2(term_ids2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type case_ids2(case_ids2SEXP);
    Rcpp::traits::input_parameter< int >::type num_cases2(num_cases2SEXP);
    rcpp_result_gen = Rcpp::wrap(sim_grid_ic(lin, anc_start, anc_stop, ancestors, info, term_ids1, case_ids1, num_cases1, term_ids2, case_ids2, num_cases2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ontologySimilarity_sample_null", (DL_FUNC) &_ontologySimilarity_sample_null, 5},
    {"_ontologySimilarity_sim_p", (DL_FUNC) &_ontologySimilarity_sim_p, 8},
    {"_ontologySimilarity_group_sim", (DL_FUNC) &_ontologySimilarity_group_sim, 4},
    {"_ontologySimilarity_sim_grid", (DL_FUNC) &_ontologySimilarity_sim_grid, 7},
    {"_ontologySimilarity_calc_term_sim_mat", (DL_FUNC) &_ontologySimilarity_calc_term_sim_mat, 6},
    {"_ontologySimilarity_sim_grid_ic", (DL_FUNC) &_ontologySimilarity_sim_grid_ic, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_ontologySimilarity(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
