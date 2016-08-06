#include "TermList.h"

using namespace Rcpp;
using namespace std;

term_list::term_list(
	IntegerVector in_term_ids, 
	IntegerVector in_case_ids, 
	int in_num_cases
) {
	term_ids = in_term_ids;
	case_ids = in_case_ids;
	num_cases = in_num_cases;

	n_terms = IntegerVector(num_cases);
	case_from = IntegerVector(num_cases);
	case_to = IntegerVector(num_cases);
	
	for (int i = 0; i < case_ids.length(); i++)
		n_terms[case_ids[i]]++;
	
	case_from[0] = 0;
	case_to[num_cases-1] = case_ids.length();
	for (int i = 0; i < (num_cases-1); i++) {
		case_from[i+1] = case_from[i] + n_terms[i];
		case_to[i] = case_from[i+1];
	}
}


