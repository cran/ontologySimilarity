#include "Sampling.h"

using namespace Rcpp;
using namespace std;

Rcpp::IntegerVector sample_int(int n, int r) {
	Rcpp::IntegerVector result(r);
	Rcpp::LogicalVector still_in(n);
	for (int i = 0; i < n; i++)
		still_in[i] = true;

	for (int i = 0; i < r; i++) {
		do {
			result[i] = random_integer(n);
		}
		while (!still_in[result[i]]);
		still_in[result[i]] = false;
	}
	return result;
}

void set_sample(
	Rcpp::IntegerVector sample,
	int set_from,
	int exc_set_to,
	int min_val,
	int exc_max_val
) {
	Rcpp::LogicalVector still_in(exc_max_val - min_val);
	for (int i = 0; i < (exc_max_val - min_val); i++)
		still_in[i] = true;

	for (int i = set_from; i < exc_set_to; i++) {
		do {
			//can't use c++ rand() if we want to share RNG state with R...
			//sample[i] = min_val + rand() % (exc_max_val - min_val);
			sample[i] = min_val + random_integer(exc_max_val - min_val);
		}
		while (!still_in[sample[i] - min_val]);
		still_in[sample[i] - min_val] = false;
	}
}

Rcpp::IntegerVector stratified_sample_int(
	Rcpp::IntegerVector strata_sizes,
	Rcpp::IntegerVector strata_sample_sizes
) {
	int num_strata = strata_sizes.length();
	int total_sample_size = 0;
	int total_items = 0;
	for (int i = 0; i < num_strata; i++) {
		total_sample_size += strata_sample_sizes[i];
		total_items += strata_sizes[i];
	}

	Rcpp::IntegerVector result(total_sample_size);
	int set_from = 0;
	int min_val = 0;

	for (int i = 0; i < num_strata; i++) {
		set_sample(
			result,
			set_from,
			set_from + strata_sample_sizes[i],
			min_val,
			min_val + strata_sizes[i]
		);

		set_from += strata_sample_sizes[i];
		min_val += strata_sizes[i];
	}

	return result;
}

simple_sampler::simple_sampler(int in_n, int in_r) {
	n = in_n;
	r = in_r;
}

Rcpp::IntegerVector simple_sampler::new_sample(void) {
	return sample_int(n, r);
}

stratified_sampler::stratified_sampler(
	Rcpp::IntegerVector in_strata_sizes,
	Rcpp::IntegerVector in_strata_sample_sizes
) {
	strata_sizes = in_strata_sizes;
	strata_sample_sizes = in_strata_sample_sizes;
}

Rcpp::IntegerVector stratified_sampler::new_sample(void) {
	return stratified_sample_int(strata_sizes, strata_sample_sizes);
}
