#include <Rcpp.h>
#include <cmath>
#include "SimMatrix.h"
#include "Sampling.h"
#include "PhenoSim.h"
#include "PVal.h"

#ifndef INCLUDED_RCPP_EXPORTS_H
#define INCLUDED_RCPP_EXPORTS_H

using namespace Rcpp;
using namespace std;

RcppExport SEXP sim_p(
	SEXP R_sim_matrix,
	SEXP R_group,
	SEXP R_min_its,
	SEXP R_max_its,
	SEXP R_signif,
	SEXP R_log_dismiss
);

RcppExport SEXP stratified_sim_p(
	SEXP R_sim_matrix,
	SEXP R_group,
	SEXP R_strata_sizes,
	SEXP R_strata_sample_sizes,
	SEXP R_min_its,
	SEXP R_max_its,
	SEXP R_signif,
	SEXP R_log_dismiss
);

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
);

RcppExport SEXP R_get_sim_matrix(
	SEXP R_term_ids,
	SEXP R_case_ids,
	SEXP R_num_cases,
	SEXP R_ttsm
);

RcppExport SEXP R_get_sim_sample(
	SEXP R_sim_mat,
	SEXP R_group_size,
	SEXP R_sample_size
);

RcppExport SEXP R_get_sim_to_profile_p(
	SEXP R_ttsm,
	SEXP R_sample_from,
	SEXP R_profile,
	SEXP R_pheno,
	SEXP R_min_its,
	SEXP R_max_its,
	SEXP R_signif,
	SEXP R_log_dismiss
);

#endif
