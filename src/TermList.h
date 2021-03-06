#include <Rcpp.h>
#include <cmath>

#ifndef TERM_LIST_H
#define TERM_LIST_H

using namespace Rcpp;
using namespace std;

struct term_list
{
	term_list(Rcpp::IntegerVector, Rcpp::IntegerVector, int);
	IntegerVector n_terms;
	IntegerVector case_from;
	IntegerVector case_to;
	IntegerVector case_ids;
	IntegerVector term_ids;
	int num_cases;
};

#endif

