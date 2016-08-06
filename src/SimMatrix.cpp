#include "SimMatrix.h"

using namespace Rcpp;
using namespace std;

NumericMatrix get_sim_matrix(
	NumericMatrix term_term_sim_mat,
	term_list terms1,
	term_list terms2
) {
	NumericMatrix result(terms1.num_cases, terms2.num_cases);
	for (int i1 = 0; i1 < terms1.num_cases; i1++) {
		for (int i2 = 0; i2 < terms2.num_cases; i2++) {
			double total = 0.0;
			for (int t1_ind = terms1.case_from[i1]; t1_ind < terms1.case_to[i1]; t1_ind++) {
				int t1 = terms1.term_ids[t1_ind];
			
				double best_term = 0.0;
				for (int t2_ind = terms2.case_from[i2]; t2_ind < terms2.case_to[i2]; t2_ind++) {
					int t2 = terms2.term_ids[t2_ind];
					best_term = max(term_term_sim_mat(t1, t2), best_term);
				}
				total += best_term;
			}
			result(i1, i2) = (terms1.n_terms[i1] == 0) ? 0.0 : (total / (double)terms1.n_terms[i1]);
		}
	}
	return result;
}

NumericMatrix get_sim_grid_ic(
	bool lin,
	IntegerVector anc_start,
	IntegerVector anc_stop,
	IntegerVector ancestors,
	NumericVector info,
	term_list terms1,
	term_list terms2
) {
	NumericMatrix result(terms1.num_cases, terms2.num_cases);
	for (int i1 = 0; i1 < terms1.num_cases; i1++) {
		for (int i2 = 0; i2 < terms2.num_cases; i2++) {
			result(i1, i2) = sim(lin, anc_start, anc_stop, ancestors, info, terms1, terms2, i1, i2);
		}
	}
	return result;
}


