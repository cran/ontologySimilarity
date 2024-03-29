#include <Rcpp.h>
#include <cmath>
#include "TermList.h"
#include "TermSetSimData.h"

#ifndef SIMILARITY_H
#define SIMILARITY_H

using namespace Rcpp;
using namespace std;

NumericMatrix get_sim_matrix(
	NumericMatrix term_term_sim_mat,
	term_list terms1,
	term_list terms2
);

NumericMatrix get_sim_grid_ic(
	bool lin,
	IntegerVector anc_start,
	IntegerVector anc_stop,
	IntegerVector ancestors,
	NumericVector info,
	term_list terms1,
	term_list terms2
);

#endif
