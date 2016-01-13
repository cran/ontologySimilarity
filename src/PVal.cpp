#include "PVal.h"

using namespace Rcpp;
using namespace std;

double bma_across_target_pheno(
	NumericMatrix ttsm,
	IntegerVector target_pheno,
	IntegerVector pheno
) {
	double total = 0.0;
	for (int i = 0; i < target_pheno.length(); i++) {
		double best = 0.0;
		for (int j = 0; j < pheno.length(); j++) {
			if (best < ttsm(target_pheno[i], pheno[j]))
				best = ttsm(target_pheno[i], pheno[j]);
		}
		total += best;
	}
	return total / (double)target_pheno.length();
}

double get_sim_to_profile_p(
	NumericMatrix ttsm,
	IntegerVector sample_from,
	IntegerVector target_pheno,
	int pheno_size,
	double similarity,
	int min_its,
	int max_its,
	double signif,
	double log_dismiss
) {
	RNGScope scope;

	int sample_from_size = sample_from.length();
	int as_sim = 0;

	for (int i = 1; i < max_its; i++) {
		IntegerVector sample(pheno_size);
		for (int j = 0; j <= pheno_size; j++)
			sample[j] = sample_from[random_integer(sample_from_size)];
		
		//Get bma... 
		if (similarity <= bma_across_target_pheno(ttsm, target_pheno, sample)) {
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
