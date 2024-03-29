#include <Rcpp.h>
#include <cmath>

#ifndef SAMPLING_H
#define SAMPLING_H

using namespace Rcpp;
using namespace std;

class Sampler
{
	public:
		virtual IntegerVector new_sample(void) = 0;
};

inline int random_integer(int exc_max)
{
	return (int)(unif_rand() * (double)exc_max) % exc_max;
}

void first_combination(IntegerVector item, size_t n);
bool next_combination(IntegerVector item, size_t n, int N);

IntegerVector sample_int(int n, int r);
void set_sample(
	IntegerVector sample,
	int set_from,
	int set_to,
	int min_val,
	int exc_max_val
);

void set_sample(
	IntegerVector sample,
	int set_from,
	int exc_set_to,
	int min_val,
	int exc_max_val
);

IntegerVector stratified_sample_int(
	IntegerVector strata_sizes,
	IntegerVector strata_sample_sizes
);

class simple_sampler : public Sampler
{
private:
	int n;
	int r;
public:
	simple_sampler(int, int);
	IntegerVector new_sample(void);
};

class stratified_sampler : public Sampler
{
private:
	IntegerVector strata_sizes;
	IntegerVector strata_sample_sizes;
public:
	stratified_sampler(IntegerVector, IntegerVector);
	IntegerVector new_sample(void);
};

#endif
