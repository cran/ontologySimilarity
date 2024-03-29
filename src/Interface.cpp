#include "Interface.h"

using namespace Rcpp;
using namespace std;

GroupSim* sim_matrix_from_data(std::string type, ReduceSim r, RObject term_sets_data) {
	GroupSim* result = NULL;
	if (type.compare("matrix") == 0) {
		result = new pre_computed_matrix(as<NumericMatrix>(term_sets_data), r);
	} 
	else if (type.compare("sim_index") == 0) {
		List l(term_sets_data);
		result = new sim_term_set_list(
			as<bool>(l["lin"]),
			as<IntegerVector>(l["start"]),
			as<IntegerVector>(l["stop"]),
			as<IntegerVector>(l["ancs"]),
			as<NumericVector>(l["info"]),
			as<IntegerVector>(l["t"]),
			as<IntegerVector>(l["c"]),
			as<int>(l["n"]),
			r,
			as<bool>(l["average_not_product"]) ? &meansim : &multiply
		);
	}
	else if (type.compare("numeric") == 0) {
		result = new VectorSim(as<NumericVector>(term_sets_data), r);
	}
	return result;
}

NumericVector sample_null(
	std::string type,
	RObject term_sets_data,
	bool use_mean,
	int group_size,
	int samples
) {
	ReduceSim r(use_mean ? &add : &worst, use_mean ? &by_size : &identity, use_mean ? 0.0 : INFINITY);

	GroupSim* sim_mat = sim_matrix_from_data(type, r, term_sets_data);
	simple_sampler simple = simple_sampler(sim_mat->population_size(), group_size);
	Sampler* sampler = &simple;
	NumericVector result = null(sampler, sim_mat, samples);
	delete sim_mat;
	return result;
}

double sim_p(
	std::string type,
	RObject term_sets_data,
	bool use_mean,
	IntegerVector group,
	int min_its,
	int max_its,
	double signif,
	double dismiss
) {
	ReduceSim r(use_mean ? &add : &worst, use_mean ? &by_size : &identity, use_mean ? 0.0 : INFINITY);

	GroupSim* sim_mat = sim_matrix_from_data(type, r, term_sets_data);

	double sim = sim_mat->groupsim(group);

	simple_sampler simple = simple_sampler(sim_mat->population_size(), group.length());

	Sampler* sampler = &simple;

	double p_val = p(
		sampler,
		sim_mat,
		sim,
		min_its,
		max_its,
		signif,
		dismiss
	);

	delete sim_mat;
	return p_val;
}

double group_sim(
	std::string type,
	RObject term_sets_data,
	bool use_mean,
	IntegerVector group
) {
	ReduceSim r(use_mean ? &add : &worst, use_mean ? &by_size : &identity, use_mean ? 0.0 : INFINITY);

	GroupSim* sim_mat = sim_matrix_from_data(type, r, term_sets_data);

	double sim = sim_mat->groupsim(group);

	delete sim_mat;
	return sim;
}

NumericMatrix sim_grid(
	IntegerVector term_ids1,
	IntegerVector case_ids1,
	int num_cases1,
	IntegerVector term_ids2,
	IntegerVector case_ids2,
	int num_cases2,
	NumericMatrix ttsm
) {
	term_list h1(term_ids1, case_ids1, num_cases1);
	term_list h2(term_ids2, case_ids2, num_cases2);
	return get_sim_matrix(ttsm, h1, h2);
}

NumericMatrix calc_term_sim_mat(
	IntegerVector anc_start,
	IntegerVector anc_stop,
	IntegerVector ancestors,
	NumericVector info,
	IntegerVector terms1,
	IntegerVector terms2
) {
	NumericMatrix result(terms1.length(), terms2.length());

	for (int i1 = 0; i1 < terms1.length(); i1++) {
		for (int i2 = 0; i2 < terms2.length(); i2++) {
			result(i1, i2) = 0.0;
			int t1 = terms1[i1];
			int t2 = terms2[i2];
			int cur_t2_anc_ind = anc_start[t2];
			for (int a1_ind = anc_start[t1]; a1_ind < anc_stop[t1]; a1_ind++) {
				int a1 = ancestors[a1_ind];
				while ((cur_t2_anc_ind < (anc_stop[t2]-1)) && (ancestors[cur_t2_anc_ind] < a1)) {
					cur_t2_anc_ind++;
				}
				if (ancestors[cur_t2_anc_ind] == a1) {
					result(i1, i2) = info[a1];
					break;
				}
			}
		}
	}
	return result;
}

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
) {
	term_list h1(term_ids1, case_ids1, num_cases1);
	term_list h2(term_ids2, case_ids2, num_cases2);

	return get_sim_grid_ic(
		lin,
		anc_start,
		anc_stop,
		ancestors,
		info,
		h1,
		h2
	);
}
