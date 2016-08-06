#include <Rcpp.h>
#include <cmath>

#include "TermList.h"

#ifndef TERMSETSIMDATA_H
#define TERMSETSIMDATA_H

using namespace Rcpp;
using namespace std;

typedef double (*combine_sim)(double, double);
typedef double (*normalise)(double, int);

double add(double a, double b);
double meansim(double a, double b);
double multiply(double a, double b);
double a(double a, double b);
double worst(double a, double b);
double by_size(double total, int n);
double identity(double a, int b);

class ReduceSim
{
public:
	combine_sim reduce;
	normalise norm;
	double sim0;
	ReduceSim(ReduceSim& r);
	ReduceSim(combine_sim, normalise, double);
};

double sim(
	bool& lin,
	IntegerVector& anc_start,
	IntegerVector& anc_stop,
	IntegerVector& ancestors,
	NumericVector& info,
	term_list& terms1,
	term_list& terms2,
	int i1,
	int i2
);

struct GroupSim {
public:
	ReduceSim reducer;
	virtual int population_size() = 0;
	virtual double groupsim(IntegerVector) = 0;
	virtual ~GroupSim();
	GroupSim(ReduceSim);
};

struct VectorSim : public GroupSim
{
public:
	NumericVector vec;
	VectorSim(NumericVector, ReduceSim);
	double groupsim(IntegerVector);
	int population_size();
};

struct SimMatrix : public GroupSim {
  public:
    virtual double pairsim(int, int) = 0;
	double groupsim(IntegerVector);
	virtual ~SimMatrix();
	SimMatrix(ReduceSim);
};

struct sim_term_set_list : public SimMatrix
{
	term_list A;
	bool lin;
	IntegerVector anc_start;
	IntegerVector anc_stop;
	IntegerVector ancestors;
	NumericVector info;
	combine_sim combine;

	sim_term_set_list(
		bool,
		IntegerVector,
		IntegerVector,
		IntegerVector,
		NumericVector,
		IntegerVector,
		IntegerVector,
		int,
		ReduceSim,
		combine_sim
	);

	double asym_sim(int, int);
	double pairsim(int, int);
	int population_size();
};

struct pre_computed_matrix : public SimMatrix
{
	NumericMatrix mat;

	pre_computed_matrix(NumericMatrix, ReduceSim);
	int population_size();
	double pairsim(int, int);
};

#endif
