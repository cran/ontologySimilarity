#include <Rcpp.h>
#include <cmath>
#include "TermList.h"

#ifndef SIMILARITY_H
#define SIMILARITY_H

using namespace Rcpp;
using namespace std;

NumericVector average_across_h(
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms
);

NumericVector average_across_phi(
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms
);

double get_phi_lik(
	bool gamma,
	IntegerVector pseudo_phi_marginal_prior_each,
	IntegerVector pseudo_phi_marginal_prior,
	LogicalMatrix row_is_column_anc,
	NumericVector lit_similarities,
	double rate,
	IntegerVector phi
);

NumericMatrix get_sim_matrix(
	NumericMatrix term_term_sim_mat,
	term_list terms1,
	term_list terms2
);

#endif
