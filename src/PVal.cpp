#include "PVal.h"

using namespace Rcpp;
using namespace std;

NumericVector null(
	Sampler* sampler,
	GroupSim* data,
	int N
) {
	NumericVector result(N);
	for (int i = 0; i < N; i++) {
		IntegerVector sample = sampler->new_sample();
		result[i] = data->groupsim(sample);
	}
	return result;
}

double p(
	Sampler* sampler,
	GroupSim* data,
	double sim,
	int min_its,
	int max_its,
	double signif,
	double log_dismiss
) {
	int as_sim = 0;
	int samples = 0;
	do {
		IntegerVector sample = sampler->new_sample();
		samples++;
		as_sim += (int)(sim <= data->groupsim(sample));
	}
	while (samples < min_its || ((R::pnorm((double)as_sim, (double)samples*signif, sqrt((double)samples*signif*(1.0-signif)),false,true) > log_dismiss) && (samples < max_its)));
	return (double)(as_sim + 1) / (double)(samples+1);
}
