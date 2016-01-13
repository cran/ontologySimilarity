#include <Rcpp.h>
#include <cmath>

#include "Sampling.h"
#include "PhenoSim.h"

#ifndef PVAL_H
#define PVAL_H

using namespace Rcpp;
using namespace std;

template <typename similarity_policy_type, typename sampling_policy_type>
double get_sim_p_val(
	similarity_policy_type similarity_policy,
	sampling_policy_type sampling_policy,
	NumericMatrix sim_matrix,
	double group_sim,
	int min_its = 1E3,
	int max_its = 1E6,
	double signif = 0.05,
	double log_dismiss = -15.0
) {
	RNGScope scope;

	int as_sim = 0;

	for (int i = 1; i < max_its; i++) {
		IntegerVector sample = sampling_policy.new_sample();
		
		if (group_sim <= similarity_policy.calc_sim(sim_matrix, sample)) {
			as_sim++;
			if ((i + 1) >= min_its)
				//use norm instead of binom as presumably lower complexity
				if (log_dismiss > p_norm(
					(double)i * signif,
					(double)i * signif * (1.0 - signif),
					(double)as_sim
				)) return (double)(as_sim + 1) / (double)(i + 1);
		}
	}

	return (double)(as_sim + 1) / (double)max_its;
}

double bma_across_target_pheno(
	NumericMatrix ttsm,
	IntegerVector target_pheno,
	IntegerVector pheno
);
	
double get_sim_to_profile_p(
	NumericMatrix ttsm,
	IntegerVector sample_from,
	IntegerVector target_pheno,
	int pheno_size,
	double similarity,
	int min_its = 1E3,
	int max_its = 1E6,
	double signif = 0.05,
	double log_dismiss = -15.0
);

template <typename similarity_policy_type, typename sampling_policy_type>
NumericVector get_sim_sample(
	similarity_policy_type similarity_policy,
	sampling_policy_type sampling_policy,
	NumericMatrix sim_matrix,
	int sample_size
) {
	RNGScope scope;

	NumericVector result(sample_size);

	for (int i = 0; i < sample_size; i++) {
		IntegerVector sample = sampling_policy.new_sample();
		
		result[i] = similarity_policy.calc_sim(sim_matrix, sample);
	}

	return result;
}

#endif
