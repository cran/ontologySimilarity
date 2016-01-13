#include "SimMatrix.h"

using namespace Rcpp;
using namespace std;

NumericVector average_across_h(
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms
) {
	NumericVector result(terms.num_cases);
	IntegerVector terms_per_case(terms.num_cases);

	int phi_terms = phi.length();
	for (int i = 0; i < terms.num_cases; i++) {
		terms_per_case[i] = 0;
		result[i] = 0.0;
	}

	for (int i = 0; i < terms.case_ids.length(); i++) {
		double best_phi_match = 0.0;
		for (int j = 0; j < phi_terms; j++)
			best_phi_match = (best_phi_match < term_term_sim_mat(terms.term_ids[i], phi[j])) ? term_term_sim_mat(terms.term_ids[i], phi[j]) : best_phi_match;

		result[terms.case_ids[i]] += best_phi_match;
		terms_per_case[terms.case_ids[i]]++;
	}

	for (int i = 0; i < terms.num_cases; i++)
		result[i] = result[i] / (double)terms_per_case[i];
	
	return result;
}

NumericVector average_across_phi(
	NumericMatrix term_term_sim_mat,
	IntegerVector phi,
	term_list terms
) {
	NumericVector result(terms.num_cases);
	for (int case_index = 0; case_index < terms.num_cases; case_index++)
		result[case_index] = 0.0;

	int phi_terms = phi.length();
	for (int phi_term_index = 0; phi_term_index < phi_terms; phi_term_index++) {
		NumericVector best_match_in_cases(terms.num_cases);
		for (int case_index = 0; case_index < terms.num_cases; case_index++) 
			best_match_in_cases[case_index] = 0.0;
		for (int term = 0; term < terms.case_ids.length(); term++)
			best_match_in_cases[terms.case_ids[term]] = (best_match_in_cases[terms.case_ids[term]] < term_term_sim_mat(terms.term_ids[term], phi[phi_term_index])) ? term_term_sim_mat(terms.term_ids[term], phi[phi_term_index]) : best_match_in_cases[terms.case_ids[term]];
		for (int case_index = 0; case_index < terms.num_cases; case_index++)
			result[case_index] += best_match_in_cases[case_index] / (double)phi_terms;
	}
	return result;
}

NumericMatrix get_sim_matrix(
	NumericMatrix term_term_sim_mat,
	term_list terms1,
	term_list terms2
) {
	NumericMatrix result(terms1.num_cases, terms2.num_cases);
	//not the most efficient way of selecting the terms for a case, but won't be the bottleneck in most cases...
	for (int case_index = 0; case_index < terms1.num_cases; case_index++) {
		int number_of_terms = 0;
		for (int case_term_ref = 0; case_term_ref < terms1.case_ids.length(); case_term_ref++) if (terms1.case_ids[case_term_ref] == case_index) number_of_terms++;
		int cursor = 0;
		IntegerVector case_terms(number_of_terms);
		for (int case_term_ref = 0; case_term_ref < terms1.case_ids.length(); case_term_ref++) if (terms1.case_ids[case_term_ref] == case_index) {
			case_terms[cursor] = terms1.term_ids[case_term_ref];
			cursor++;
		}
		NumericVector sims_to_case = average_across_phi(term_term_sim_mat, case_terms, terms2);	
		for (int i = 0; i < terms2.num_cases; i++)
			result(case_index, i) = sims_to_case[i];
	}

	return result;
}
