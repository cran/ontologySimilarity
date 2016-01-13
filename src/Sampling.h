#include <Rcpp.h>
#include <cmath>

#ifndef SAMPLING_H
#define SAMPLING_H

using namespace Rcpp;
using namespace std;

inline int random_integer(int exc_max)
{
	return (int)(unif_rand() * (double)exc_max) % exc_max;
}

Rcpp::IntegerVector sample_int(int n, int r);
void set_sample(
	Rcpp::IntegerVector sample,
	int set_from,
	int set_to,
	int min_val,
	int exc_max_val
);

void set_sample(
	Rcpp::IntegerVector sample,
	int set_from,
	int exc_set_to,
	int min_val,
	int exc_max_val
);

Rcpp::IntegerVector stratified_sample_int(
	Rcpp::IntegerVector strata_sizes,
	Rcpp::IntegerVector strata_sample_sizes
);

class simple_sampler
{
private:
	int n;
	int r;
public:
	simple_sampler(int, int);
	Rcpp::IntegerVector new_sample(void);
};

class stratified_sampler
{
private:
	Rcpp::IntegerVector strata_sizes;
	Rcpp::IntegerVector strata_sample_sizes;
public:
	stratified_sampler(Rcpp::IntegerVector, Rcpp::IntegerVector);
	Rcpp::IntegerVector new_sample(void);
};

#endif
