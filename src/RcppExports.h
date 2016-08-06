#include <Rcpp.h>
#include <cmath>
#include "SimMatrix.h"
#include "Sampling.h"
#include "TermSetSimData.h"
#include "PVal.h"

#ifndef INCLUDED_RCPP_EXPORTS_H
#define INCLUDED_RCPP_EXPORTS_H

using namespace Rcpp;
using namespace std;

GroupSim* sim_matrix_from_data(int type, ReduceSim r, SEXP term_sets_data);

RcppExport SEXP sample_null(
	SEXP data_type,
	SEXP term_sets_data,
	SEXP R_mean_not_min,
	SEXP group_size,
	SEXP samples
);

RcppExport SEXP sim_p(
	SEXP data_type,
	SEXP term_sets_data,
	SEXP R_mean_not_min,
	SEXP R_group,
	SEXP min_its,
	SEXP max_its,
	SEXP signif,
	SEXP dismiss
);

RcppExport SEXP group_sim(
	SEXP data_type,
	SEXP term_sets_data,
	SEXP R_mean_not_min,
	SEXP R_group
);

RcppExport SEXP R_get_sim_grid(
	SEXP R_term_ids1,
	SEXP R_case_ids1,
	SEXP R_num_cases1,
	SEXP R_term_ids2,
	SEXP R_case_ids2,
	SEXP R_num_cases2,
	SEXP R_ttsm
);

RcppExport SEXP calc_term_sim_mat(
	SEXP R_anc_start,
	SEXP R_anc_stop,
	SEXP R_ancestors,
	SEXP R_info,
	SEXP R_terms1,
	SEXP R_terms2
);

RcppExport SEXP R_get_sim_grid_ic(
	SEXP lin,
	SEXP anc_start,
	SEXP anc_stop,
	SEXP ancestors,
	SEXP info,
	SEXP R_term_ids1,
	SEXP R_case_ids1,
	SEXP R_num_cases1,
	SEXP R_term_ids2,
	SEXP R_case_ids2,
	SEXP R_num_cases2
);

#endif
