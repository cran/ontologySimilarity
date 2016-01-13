#include <Rcpp.h>
#include <cmath>
#include "Sampling.h"

#ifndef INCLUDED_PHENO_SIM_H
#define INCLUDED_PHENO_SIM_H

typedef double (*similarity_measure)(Rcpp::NumericMatrix,Rcpp::IntegerVector);

using namespace Rcpp;
using namespace std;

double mean_sim(
	Rcpp::NumericMatrix sim_matrix,
	Rcpp::IntegerVector group
);

double min_sim(
	Rcpp::NumericMatrix sim_matrix,
	Rcpp::IntegerVector group
);

void first_combination(Rcpp::IntegerVector item, size_t n);
bool next_combination(Rcpp::IntegerVector item, size_t n, size_t N);

double best_n(
	similarity_measure sim_measure,
	int n,
	Rcpp::NumericMatrix sim_matrix,
	Rcpp::IntegerVector group
);

double p_norm(
	double mean,
	double variance,
	double value
);

class simple_measure
{
private:
	similarity_measure measure;
public:
	simple_measure(similarity_measure in_measure);
	double calc_sim(Rcpp::NumericMatrix, Rcpp::IntegerVector);
};

class best_subgroup_measure
{
private:
	similarity_measure measure;
	int subgroup_size;
public:
	best_subgroup_measure(similarity_measure, int);
	double calc_sim(Rcpp::NumericMatrix, Rcpp::IntegerVector);
};

#endif

