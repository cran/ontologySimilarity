#include "RcppExports.h"

using namespace Rcpp;
using namespace std;

RcppExport SEXP sim_p(
	SEXP R_sim_matrix,
	SEXP R_group,
	SEXP R_min_its,
	SEXP R_max_its,
	SEXP R_signif,
	SEXP R_log_dismiss
) {
BEGIN_RCPP
	Rcpp::NumericMatrix sim_matrix(R_sim_matrix);
	Rcpp::IntegerVector group(R_group);

	int min_its = Rcpp::as<int>(R_min_its);
	int max_its = Rcpp::as<int>(R_max_its);
	double signif = Rcpp::as<double>(R_signif);
	double log_dismiss = Rcpp::as<double>(R_log_dismiss);

	simple_measure sim(mean_sim);
	simple_sampler sampler(sim_matrix.nrow(), group.length());

	double group_sim = sim.calc_sim(sim_matrix, group);

	return Rcpp::NumericVector::create(
		get_sim_p_val<simple_measure, simple_sampler>(
			sim,
			sampler,
			sim_matrix,
			group_sim,
			min_its,
			max_its,
			signif,
			log_dismiss
		)
	);
END_RCPP
}

RcppExport SEXP stratified_sim_p(
	SEXP R_sim_matrix,
	SEXP R_group,
	SEXP R_strata_sizes,
	SEXP R_strata_sample_sizes,
	SEXP R_min_its,
	SEXP R_max_its,
	SEXP R_signif,
	SEXP R_log_dismiss
) {
BEGIN_RCPP
	Rcpp::NumericMatrix sim_matrix(R_sim_matrix);
	Rcpp::IntegerVector group(R_group);
	Rcpp::IntegerVector strata_sizes(R_strata_sizes);
	Rcpp::IntegerVector strata_sample_sizes(R_strata_sample_sizes);

	int min_its = Rcpp::as<int>(R_min_its);
	int max_its = Rcpp::as<int>(R_max_its);
	double signif = Rcpp::as<double>(R_signif);
	double log_dismiss = Rcpp::as<double>(R_log_dismiss);

	simple_measure sim(mean_sim);
	stratified_sampler sampler(strata_sizes, strata_sample_sizes);

	double group_sim = sim.calc_sim(sim_matrix, group);

	return Rcpp::NumericVector::create(
		get_sim_p_val<simple_measure, stratified_sampler>(
			sim,
			sampler,
			sim_matrix,
			group_sim,
			min_its,
			max_its,
			signif,
			log_dismiss
		)
	);
END_RCPP
}

RcppExport SEXP stratified_best_subgroup_p(
	SEXP R_subgroup_size,
	SEXP R_sim_matrix,
	SEXP R_group,
	SEXP R_strata_sizes,
	SEXP R_strata_sample_sizes,
	SEXP R_min_its,
	SEXP R_max_its,
	SEXP R_signif,
	SEXP R_log_dismiss
) {
BEGIN_RCPP
	int subgroup_size = Rcpp::as<int>(R_subgroup_size);
	Rcpp::NumericMatrix sim_matrix(R_sim_matrix);
	Rcpp::IntegerVector group(R_group);
	Rcpp::IntegerVector strata_sizes(R_strata_sizes);
	Rcpp::IntegerVector strata_sample_sizes(R_strata_sample_sizes);

	int min_its = Rcpp::as<int>(R_min_its);
	int max_its = Rcpp::as<int>(R_max_its);
	double signif = Rcpp::as<double>(R_signif);
	double log_dismiss = Rcpp::as<double>(R_log_dismiss);

	best_subgroup_measure sim(mean_sim, subgroup_size);
	stratified_sampler sampler(strata_sizes, strata_sample_sizes);

	double group_sim = sim.calc_sim(sim_matrix, group);

	return Rcpp::NumericVector::create(
		get_sim_p_val<best_subgroup_measure, stratified_sampler>(
			sim,
			sampler,
			sim_matrix,
			group_sim,
			min_its,
			max_its,
			signif,
			log_dismiss
		)
	);
END_RCPP
}

RcppExport SEXP R_get_sim_matrix(
	SEXP R_term_ids,
	SEXP R_case_ids,
	SEXP R_num_cases,
	SEXP R_ttsm
) {
BEGIN_RCPP
	NumericMatrix ttsm(R_ttsm);

	IntegerVector term_ids(R_term_ids);
	IntegerVector case_ids(R_case_ids);

	term_list h(term_ids, case_ids, as<int>(R_num_cases));

	return get_sim_matrix(ttsm, h, h);
END_RCPP
}

RcppExport SEXP R_get_sim_matrix_2_sets(
	SEXP R_term_ids1,
	SEXP R_case_ids1,
	SEXP R_num_cases1,
	SEXP R_term_ids2,
	SEXP R_case_ids2,
	SEXP R_num_cases2,
	SEXP R_ttsm
) {
BEGIN_RCPP
	NumericMatrix ttsm(R_ttsm);

	IntegerVector term_ids1(R_term_ids1);
	IntegerVector case_ids1(R_case_ids1);

	term_list h1(term_ids1, case_ids1, as<int>(R_num_cases1));

	IntegerVector term_ids2(R_term_ids2);
	IntegerVector case_ids2(R_case_ids2);

	term_list h2(term_ids2, case_ids2, as<int>(R_num_cases2));

	return get_sim_matrix(ttsm, h1, h2);
END_RCPP
}

RcppExport SEXP R_get_sim_sample(
	SEXP R_sim_mat,
	SEXP R_group_size,
	SEXP R_sample_size
) {
BEGIN_RCPP
	int group_size = as<int>(R_group_size);
	int sample_size = as<int>(R_sample_size);

	NumericMatrix sim_mat(R_sim_mat);
	simple_measure sim(mean_sim);
	simple_sampler sampler(sim_mat.nrow(), group_size);

	return get_sim_sample<simple_measure, simple_sampler>(sim, sampler, sim_mat, sample_size);
END_RCPP
}

RcppExport SEXP R_get_sim_to_profile_p(
	SEXP R_ttsm,
	SEXP R_sample_from,
	SEXP R_profile,
	SEXP R_pheno,
	SEXP R_min_its,
	SEXP R_max_its,
	SEXP R_signif,
	SEXP R_log_dismiss
) {
BEGIN_RCPP
	NumericMatrix ttsm = as<NumericMatrix>(R_ttsm);
	IntegerVector sample_from = as<IntegerVector>(R_sample_from);
	IntegerVector profile = as<IntegerVector>(R_profile);
	IntegerVector pheno = as<IntegerVector>(R_pheno);
	double bma = bma_across_target_pheno(
		ttsm,
		profile,
		pheno
	);

	return wrap<double>(get_sim_to_profile_p(
		ttsm,
		sample_from,
		profile,
		pheno.length(),
		bma,
		as<int>(R_min_its),
		as<int>(R_max_its),
		as<double>(R_signif),
		as<double>(R_log_dismiss)
	));
END_RCPP
}

