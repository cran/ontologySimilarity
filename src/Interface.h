#include <Rcpp.h>
#include <cmath>
#include "SimMatrix.h"
#include "Sampling.h"
#include "TermSetSimData.h"
#include "PVal.h"

#ifndef INCLUDED_INTERFACE_H
#define INCLUDED_INTERFACE_H

using namespace Rcpp;
using namespace std;

GroupSim* sim_matrix_from_data(int type, ReduceSim r, RObject term_sets_data);

// [[Rcpp::export]]
NumericVector sample_null(
	std::string type,
	RObject term_sets_data,
	bool use_mean,
	int group_size,
	int samples
);

// [[Rcpp::export]]
double sim_p(
	std::string type,
	RObject term_sets_data,
	bool use_mean,
	IntegerVector group,
	int min_its,
	int max_its,
	double signif,
	double dismiss
);

// [[Rcpp::export]]
double group_sim(
	std::string type,
	RObject term_sets_data,
	bool use_mean,
	IntegerVector group
);

// [[Rcpp::export]]
NumericMatrix sim_grid(
	IntegerVector term_ids1,
	IntegerVector case_ids1,
	int num_cases1,
	IntegerVector term_ids2,
	IntegerVector case_ids2,
	int num_cases2,
	NumericMatrix ttsm
);

// [[Rcpp::export]]
NumericMatrix calc_term_sim_mat(
	IntegerVector anc_start,
	IntegerVector anc_stop,
	IntegerVector ancestors,
	NumericVector info,
	IntegerVector terms1,
	IntegerVector terms2
);

// [[Rcpp::export]]
NumericMatrix sim_grid_ic(
	bool lin,
	IntegerVector anc_start,
	IntegerVector anc_stop,
	IntegerVector ancestors,
	NumericVector info,
	IntegerVector term_ids1,
	IntegerVector case_ids1,
	int num_cases1,
	IntegerVector term_ids2,
	IntegerVector case_ids2,
	int num_cases2
);

#endif
