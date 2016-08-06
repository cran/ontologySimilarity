#include "RcppExports.h"

using namespace Rcpp;
using namespace std;

GroupSim* sim_matrix_from_data(string type, ReduceSim r, SEXP term_sets_data) {
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

RcppExport SEXP sample_null(
	SEXP data_type,
	SEXP term_sets_data,
	SEXP R_mean_not_min,
	SEXP group_size,
	SEXP samples
) {
BEGIN_RCPP
	RNGScope scope;
	string type = as<string>(data_type);

	bool use_mean = as<bool>(R_mean_not_min);

	ReduceSim r(use_mean ? &add : &worst, use_mean ? &by_size : &identity, use_mean ? 0.0 : INFINITY);

	GroupSim* sim_mat = sim_matrix_from_data(type, r, term_sets_data);
	simple_sampler simple = simple_sampler(sim_mat->population_size(), as<int>(group_size));
	Sampler* sampler = &simple;
	NumericVector result = null(sampler, sim_mat, as<int>(samples));
	delete sim_mat;
	return result;
END_RCPP
}

RcppExport SEXP sim_p(
	SEXP data_type,
	SEXP term_sets_data,
	SEXP R_mean_not_min,
	SEXP R_group,
	SEXP min_its,
	SEXP max_its,
	SEXP signif,
	SEXP dismiss
) {
BEGIN_RCPP
	RNGScope scope;
	Rcpp::RObject __result;

	string type = as<string>(data_type);

	bool use_mean = as<bool>(R_mean_not_min);

	ReduceSim r(use_mean ? &add : &worst, use_mean ? &by_size : &identity, use_mean ? 0.0 : INFINITY);

	GroupSim* sim_mat = sim_matrix_from_data(type, r, term_sets_data);

	IntegerVector group = as<IntegerVector>(R_group);

	double sim = sim_mat->groupsim(group);

	simple_sampler simple = simple_sampler(sim_mat->population_size(), group.length());

	Sampler* sampler = &simple;

	double p_val = p(
		sampler,
		sim_mat,
		sim,
		as<int>(min_its),
		as<int>(max_its),
		as<double>(signif),
		as<double>(dismiss)
	);

	delete sim_mat;

	__result = wrap<double>(p_val);
	return __result;
END_RCPP
}

RcppExport SEXP group_sim(
	SEXP data_type,
	SEXP term_sets_data,
	SEXP R_mean_not_min,
	SEXP R_group
) {
BEGIN_RCPP
	Rcpp::RObject __result;
	string type = as<string>(data_type);

	bool use_mean = as<bool>(R_mean_not_min);

	ReduceSim r(use_mean ? &add : &worst, use_mean ? &by_size : &identity, use_mean ? 0.0 : INFINITY);

	GroupSim* sim_mat = sim_matrix_from_data(type, r, term_sets_data);

	IntegerVector group = as<IntegerVector>(R_group);

	double sim = sim_mat->groupsim(group);

	delete sim_mat;

	__result = wrap<double>(sim);
	return __result;
END_RCPP
}

RcppExport SEXP R_get_sim_grid(
	SEXP R_term_ids1,
	SEXP R_case_ids1,
	SEXP R_num_cases1,
	SEXP R_term_ids2,
	SEXP R_case_ids2,
	SEXP R_num_cases2,
	SEXP R_ttsm
) {
BEGIN_RCPP
	Rcpp::RObject __result;
	NumericMatrix ttsm(R_ttsm);

	IntegerVector term_ids1(R_term_ids1);
	IntegerVector case_ids1(R_case_ids1);

	term_list h1(term_ids1, case_ids1, as<int>(R_num_cases1));

	IntegerVector term_ids2(R_term_ids2);
	IntegerVector case_ids2(R_case_ids2);

	term_list h2(term_ids2, case_ids2, as<int>(R_num_cases2));

	__result = get_sim_matrix(ttsm, h1, h2);
	return __result;
END_RCPP
}

RcppExport SEXP calc_term_sim_mat(
	SEXP R_anc_start,
	SEXP R_anc_stop,
	SEXP R_ancestors,
	SEXP R_info,
	SEXP R_terms1,
	SEXP R_terms2
) {
	IntegerVector anc_start = as<IntegerVector>(R_anc_start);
	IntegerVector anc_stop = as<IntegerVector>(R_anc_stop);
	IntegerVector ancestors = as<IntegerVector>(R_ancestors);
	NumericVector info = as<NumericVector>(R_info);
	IntegerVector terms1 = as<IntegerVector>(R_terms1);
	IntegerVector terms2 = as<IntegerVector>(R_terms2);

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

RcppExport SEXP R_get_sim_grid_ic(
	SEXP lin,
	SEXP anc_start,
	SEXP anc_stop,
	SEXP ancestors,
	SEXP info,
	SEXP term_ids1,
	SEXP case_ids1,
	SEXP num_cases1,
	SEXP term_ids2,
	SEXP case_ids2,
	SEXP num_cases2
) {
	Rcpp::RObject __result;
	term_list h1(as<IntegerVector>(term_ids1), as<IntegerVector>(case_ids1), as<int>(num_cases1));
	term_list h2(as<IntegerVector>(term_ids2), as<IntegerVector>(case_ids2), as<int>(num_cases2));

	__result = get_sim_grid_ic(
		as<bool>(lin),
		as<IntegerVector>(anc_start),
		as<IntegerVector>(anc_stop),
		as<IntegerVector>(ancestors),
		as<NumericVector>(info),
		h1,
		h2
	);
	return __result;
}
