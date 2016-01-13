#include "PhenoSim.h"

using namespace Rcpp;
using namespace std;

double mean_sim(
	Rcpp::NumericMatrix sim_matrix,
	Rcpp::IntegerVector group
) {
	double result = 0.0;
	int n = group.length();
	for (int row = 1; row < n; row++)
		for (int col = 0; col < row; col++)
			result += sim_matrix(group[row], group[col]);
	return result / ((double)(n * (n - 1)) * 0.5);
}

double min_sim(
	Rcpp::NumericMatrix sim_matrix,
	Rcpp::IntegerVector group
) {
	double result = 0.0;
	int n = group.length();
	for (int row = 1; row < n; row++)
		for (int col = 0; col < row; col++) {
			double next_contestant = sim_matrix(group[row], group[col]);
			if (next_contestant < result)
				result = next_contestant;
		}

	return result;
}

void first_combination(Rcpp::IntegerVector item, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        item[i] = i;
    }
}

bool next_combination(Rcpp::IntegerVector item, size_t n, size_t N)
{
    for (size_t i = 1; i <= n; ++i) {
        if (item[n-i] < N-i) {
            ++item[n-i];
            for (size_t j = n-i+1; j < n; ++j) {
                item[j] = item[j-1] + 1;
            }
            return true;
        }
    }
    return false;
}

double best_n(
	similarity_measure sim_measure,
	int n,
	Rcpp::NumericMatrix sim_matrix,
	Rcpp::IntegerVector group
) {
	double result = 0.0;
	int group_size = group.length();
	Rcpp::IntegerVector test_group(n);
	first_combination(test_group, n);
	do {
		double test_sim = sim_measure(sim_matrix, group[test_group]);
		if (test_sim > result)
			result = test_sim;
	} 
	while (next_combination(test_group, n, group_size));

	return result;
}

double p_norm(
	double mean,
	double variance,
	double value
) {
	return R::pnorm(
		value,
		mean,
		sqrt(variance),
		0,
		1	
	);
}

simple_measure::simple_measure(similarity_measure in_measure) {
	measure = in_measure;
}

double simple_measure::calc_sim(Rcpp::NumericMatrix m, Rcpp::IntegerVector g) {
	return measure(m, g);
}

best_subgroup_measure::best_subgroup_measure(
	similarity_measure in_measure,
	int in_subgroup_size
) {
	subgroup_size = in_subgroup_size;
	measure = in_measure;
}

double best_subgroup_measure::calc_sim(Rcpp::NumericMatrix m, Rcpp::IntegerVector g) {
	return best_n(measure, subgroup_size, m, g);
}
