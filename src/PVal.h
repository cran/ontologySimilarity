#include <Rcpp.h>
#include <cmath>

#include "Sampling.h"
#include "TermSetSimData.h"

#ifndef PVAL_H
#define PVAL_H

using namespace Rcpp;
using namespace std;

NumericVector null(
	Sampler*,
	GroupSim*,
	int
);

double p(
	Sampler*,
	GroupSim*,
	double,
	int,
	int,
	double,
	double
);

#endif
