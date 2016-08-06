#include "TermSetSimData.h"

using namespace Rcpp;
using namespace std;

double add(double a, double b) { return a + b; }
double meansim(double a, double b) { return (a + b) / 2.0; }
double multiply(double a, double b) { return a * b; }
double a(double a, double b) { return a; }
double worst(double a, double b) { return min(a, b); }
double by_size(double total, int n) { return total / (double)n; }
double identity(double a, int b) { return a; }

GroupSim::~GroupSim() {}
SimMatrix::~SimMatrix() {}

ReduceSim::ReduceSim(ReduceSim& r) {
	norm = r.norm;
	reduce = r.reduce;
	sim0 = r.sim0;
}

ReduceSim::ReduceSim(combine_sim r, normalise n, double ini) {
	sim0 = ini;
	reduce = r;
	norm = n;
}

GroupSim::GroupSim(ReduceSim r) : reducer(r) { }

VectorSim::VectorSim(NumericVector v, ReduceSim r) : GroupSim(r) {
	vec = v;
}

int VectorSim::population_size() { return vec.length(); }

double VectorSim::groupsim(IntegerVector group) {
	double agg = reducer.sim0;
	int n = group.length();
	for (int i = 0; i < n; i++) {
		double sim = vec[group[i]];
		agg = reducer.reduce(agg, sim);
	}

	return reducer.norm(agg, n);
}

SimMatrix::SimMatrix(ReduceSim r) : GroupSim(r) {}

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
) {
	double total = 0.0;
	for (int t1_ind = terms1.case_from[i1]; t1_ind < terms1.case_to[i1]; t1_ind++) {
		int t1 = terms1.term_ids[t1_ind];
		double best_term = 0.0;
		for (int t2_ind = terms2.case_from[i2]; t2_ind < terms2.case_to[i2]; t2_ind++) {
			double best_anc = 0.0;
			int t2 = terms2.term_ids[t2_ind];
			int cur_t2_anc_ind = anc_start[t2];
			for (int a1_ind = anc_start[t1]; a1_ind < anc_stop[t1]; a1_ind++) {
				int a1 = ancestors[a1_ind];
				if (a1 < ancestors[cur_t2_anc_ind]) continue;
				while ((cur_t2_anc_ind < (anc_stop[t2]-1)) && (ancestors[cur_t2_anc_ind] < a1)) {
					cur_t2_anc_ind++;
				}
				if (ancestors[cur_t2_anc_ind] == a1) {
					best_anc = info[a1];
					break;
				}
			}

			double score;
			if (lin) {
				if (best_anc > 0.0)
					score = (2.0 * best_anc / (info[t1] + info[t2]));
				else 
					score = 0.0;
			}
			else {
				score = best_anc;
			}

			best_term = max(score, best_term);
		}
		total += best_term;
	}
	return (terms1.n_terms[i1] == 0) ? 0.0 : (total / (double)terms1.n_terms[i1]);
}


double SimMatrix::groupsim(
	IntegerVector group
) {
	double agg = reducer.sim0;
	int n = group.length();
	for (int row = 1; row < n; row++)
		for (int col = 0; col < row; col++) {
			double sim = pairsim(group[row], group[col]);
			agg = reducer.reduce(agg, sim);
		}

	return reducer.norm(agg, ((double)(n * (n - 1)) * 0.5));
}

sim_term_set_list::sim_term_set_list(
	bool in_lin,
	IntegerVector in_anc_start,
	IntegerVector in_anc_stop,
	IntegerVector in_ancestors,
	NumericVector in_info,
	IntegerVector in_t1,
	IntegerVector in_c1,
	int in_n1,
	ReduceSim r,
	combine_sim comb
) : SimMatrix(r), A(in_t1, in_c1, in_n1) {
	lin = in_lin;
	anc_start = in_anc_start;
	anc_stop = in_anc_stop;
	ancestors = in_ancestors;
	info = in_info;
	combine = comb;
}

pre_computed_matrix::pre_computed_matrix(NumericMatrix m, ReduceSim r) : SimMatrix(r) {
	mat = m;
}

double pre_computed_matrix::pairsim(int i1, int i2) { return mat(i1, i2); }
int pre_computed_matrix::population_size() {
	return mat.nrow();
}

double sim_term_set_list::asym_sim(int i1, int i2) {
	return sim(lin, anc_start, anc_stop, ancestors, info, A, A, i1, i2);
}

double sim_term_set_list::pairsim(int i1, int i2) {
	return combine(asym_sim(i1, i2),asym_sim(i2, i1));
}

int sim_term_set_list::population_size() {
	return A.num_cases;
}
