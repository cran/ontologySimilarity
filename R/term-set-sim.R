#' Get matrix of similarity rank from similarity matrix
#'
#' Given a lower triangular similarity matrix, construct a distance matrix where the rows are the ranks of the column cases with respect to similarity to the row case. If relative similarity is of interest, this rank-transformation may reduce bias in favour of high similarity scores in downstream analysis.
#'
#' @param similarity_matrix Lower triangular numeric matrix of similarities, where the rownames and colnames are identical to the case IDs.
#' @param symmetric Logical value determining whether to `symmetrify' resultant matrix by averaging rank similarity of A -> B and B -> A.
#' @return Matrix of rank similarities.
#' @export
get_similarity_rank_matrix <- function(similarity_matrix, symmetric=TRUE) {
	similarity_matrix[upper.tri(similarity_matrix, diag=FALSE)] <- 0
	similarity_matrix <- similarity_matrix + t(similarity_matrix)
	diag(similarity_matrix) <- NA
	ranked.sim.matrix <- apply(similarity_matrix, 1, function(x) rank(-x, ties.method="max", na.last="keep"))
	if (symmetric)
		(t(ranked.sim.matrix) + ranked.sim.matrix)/2
	else
		t(ranked.sim.matrix)
}

#' Get `term sets to term' similarity matrix
#'
#' Create a numeric matrix of similarities between term sets and individual terms.
#'
#' @template terms
#' @template term_sets
#' @param ... Other arguments to be passed to \code{\link{get_sim_grid}}.
#' @return Numeric matrix of term set-to-term similarities
#' @seealso \code{\link{get_sim_grid}}
#' @export
get_term_set_to_term_sims <- function(term_sets, terms, ...) {
	get_sim_grid(
		...,
		term_sets=term_sets,
		term_sets2=setNames(nm=terms, as.list(terms)),
		combine=function(x, y) x
	)
}

#' Get asymmetrical similarity matrix
#'
#' Create a numeric matrix of similarities between two lists of term sets, but only averaging over the terms in sets from \code{A} the similarities of the best matches in sets from \code{B}.
#'
#' @param A List of term sets.
#' @param B List of term sets.
#' @param ... Other arguments to be passed to \code{\link{get_sim_grid}}.
#' @return Numeric matrix of similarities
#' @seealso \code{\link{get_sim_grid}} \code{\link{get_profile_sims}}
#' @export
get_asym_sim_grid <- function(A, B, ...) {
	get_sim_grid(
		...,
		term_sets=A,
		term_sets2=B,
		combine=function(x, y) x
	)
}

#' Get similarities of term sets to profile
#'
#' Get numeric vector of similarities between each item in a list of term sets and another `ontological profile', i.e. a single term set. Similarity averaging over terms in \code{term_sets}.
#'
#' @param profile Character vector of term IDs.
#' @template term_sets
#' @param ... Other arguments to pass to \code{\link{get_sim_grid}}.
#' @return Numeric vector of profile similarities.
#' @seealso \code{\link{get_asym_sim_grid}} \code{\link{get_sim_grid}}
#' @export
get_profile_sims <- function(profile, term_sets, ...) {
	sims <- get_asym_sim_grid(..., A=term_sets, B=list(profile))
	setNames(as.numeric(sims), nm=rownames(sims))
}

#' Get similarity matrix of pairwise similarities of term sets.
#'
#' Using either an \code{ontology_index} object and numeric vector of information content per term - or a matrix of between-term similarities (e.g. the output of \code{\link{get_term_sim_mat}}), create a numeric matrix of `between-term set' similarities. Either the `best-match-average' or `best-match-product' approach (i.e. where the 2 scores obtained by applying the asymmetric `best-match' similarity function to two term sets in each order are combined by taking the average or the product respectively). Either Lin's (default) or Resnik's definition of term similarity can be used. If \code{information_content} is not specified, a default value from \code{\link{descendants_IC}} is generated.
#'
#' Note that if any term set within \code{term_sets} has 0 terms associated with it, it will get a similarity of 0 to any other set. If you do not want to compare term sets with no annotation, take care to filter out empty sets first, e.g. by `term_sets=term_sets[sapply(term_sets, length) > 0]`.
#'
#' @template ontology
#' @template information_content
#' @template term_sim_method
#' @template term_sim_mat
#' @template term_sets
#' @param term_sets2 Second set of term sets.
#' @template combine
#' @return Numeric matrix of pairwise term set similarities.
#' @seealso \code{\link{get_term_sim_mat}} \code{\link{get_sim_p}} \code{\link{get_asym_sim_grid}}
#' @export
#' @examples
#' library(ontologyIndex)
#' data(hpo)
#' get_sim_grid(ontology=hpo, term_sets=list(
#'   `case 1`=c("HP:0001873","HP:0011877"),
#'   `case 2`=c("HP:0001892","HP:0001873"),
#'   `case 3`=c("HP:0001872","HP:0000707")))
get_sim_grid <- function(ontology, information_content, term_sim_method, term_sim_mat, term_sets, term_sets2=term_sets, combine="average") {
	stopifnot(class(term_sets) == "list")
	stopifnot(class(term_sets2) == "list")
	if (!xor(missing(ontology), missing(term_sim_mat))) stop("Must pass either a 'term_sim_mat' or 'ontology' argument")
	if (is.character(combine)) stopifnot(combine %in% c("average","product"))
	if (!is.character(combine)) stopifnot(is.function(combine))

	if (!missing(term_sim_mat)) {
		get_sim_grid_from_tsm(term_sim_mat=term_sim_mat, term_sets=term_sets, term_sets2=term_sets2, combine=combine)
	} else {
		if (missing(information_content))
			information_content <- descendants_IC(ontology)

		if (missing(term_sim_method))
			term_sim_method <- "lin"

		get_sim_grid_from_ic(ontology, information_content, term_sets, term_sets2, combine=combine, term_sim_method=term_sim_method)
	}
}

#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_sim_grid_from_ic <- function(
	ontology,
	information_content,
	term_sets,
	term_sets2,
	combine,
	term_sim_method
) {
	stopifnot(class(ontology) == "ontology_index")
	stopifnot(is.numeric(information_content))
	stopifnot(term_sim_method %in% c("lin", "resnik"))

	args <- components_for_calc_sim_from_ic(ontology, information_content, list(term_sets, term_sets2))

	result <- (if (is.character(combine)) (if (combine == "average") function(x, y) { (x + y) / 2 } else match.fun("*")) else combine)(
		do.call(what=.Call, c(list("R_get_sim_grid_ic", PACKAGE="ontologySimilarity", term_sim_method == "lin"), args[c("start", "stop", "ancs", "info")], args[["term_sets"]][[1]], args[["term_sets"]][[2]])),
		t(do.call(what=.Call, c(list("R_get_sim_grid_ic", PACKAGE="ontologySimilarity", term_sim_method == "lin"), args[c("start", "stop", "ancs", "info")], args[["term_sets"]][[2]], args[["term_sets"]][[1]])))
	)

	rownames(result) <- names(term_sets)
	colnames(result) <- names(term_sets2)
	result
}

#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_sim_grid_from_tsm <- function(
	term_sim_mat,
	term_sets,
	term_sets2=term_sets,
	combine
) {
	if (!all(unlist(use.names=FALSE, term_sets) %in% rownames(term_sim_mat)))
		stop("Terms in 'term_sets' missing from rownames of 'term_sim_mat'")
	if (!all(unlist(use.names=FALSE, term_sets2) %in% colnames(term_sim_mat)))
		stop("Terms in 'term_sets2' missing from colnames of 'term_sim_mat'")

	t1 <- as.integer(match(unlist(use.names=FALSE, term_sets), rownames(term_sim_mat)))-1
	c1 <- unlist(use.names=FALSE, mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(term_sets)-1), sapply(term_sets, length)))
	n1 <- length(term_sets)

	t2 <- as.integer(match(unlist(use.names=FALSE, term_sets2), colnames(term_sim_mat)))-1
	c2 <- unlist(use.names=FALSE, mapply(SIMPLIFY=FALSE, FUN=rep, 0:(length(term_sets2)-1), sapply(term_sets2, length)))
	n2 <- length(term_sets2)
	result <- (if (is.character(combine)) (if (combine == "average") function(x, y) { (x + y) / 2 } else match.fun("*")) else combine)(
		.Call("R_get_sim_grid", PACKAGE="ontologySimilarity", t1, c1, n1, t2, c2, n2, term_sim_mat),
		t(.Call("R_get_sim_grid", PACKAGE="ontologySimilarity", t2, c2, n2, t1, c1, n1, term_sim_mat))
	)

	rownames(result) <- names(term_sets)
	colnames(result) <- names(term_sets2)
	result
}

resnik_asym <- function(
	ontology,
	information_content,
	term_set_1,
	term_set_2
) {
	mean(sapply(
		term_set_1,	
		function(x) max(c(0, information_content[
			intersect(ontology$ancestors[[x]], get_ancestors(ontology, term_set_2))
		]))
	))
}

#' Calculate Resnik similarity score of two term sets
#'
#' Warning! This function is slow - performing large numbers of `between term-set' similarity calculations should be done using \code{\link{get_sim_grid}}.
#'
#' @template ontology
#' @template information_content
#' @param term_set_1 Character vector of terms.
#' @param term_set_2 Character vector of terms.
#' @return Numeric value.
#' @seealso \code{\link{lin}}, \code{\link{get_term_sim_mat}}
#' @references Resnik, P. (1995). `Using information content to evaluate semantic similarity in a taxonomy'. Proceedings of the 14th IJCAI 1, 448-453.
#' @export
resnik <- function(
	ontology,
	information_content,
	term_set_1,
	term_set_2
) {
	if (min(c(length(term_set_1), length(term_set_2))) == 0)
		return(0)
	(resnik_asym(ontology, information_content, term_set_1, term_set_2) + resnik_asym(ontology, information_content, term_set_2, term_set_1))/2	
}

lin_asym <- function(
	ontology,
	information_content,
	term_set_1,
	term_set_2
) {
	mean(sapply(
		term_set_1,	
		function(t1) {
			max(sapply(term_set_2, function(t2) 2*max(c(0, information_content[intersect(ontology$ancestors[[t2]], ontology$ancestors[[t1]])]))/(information_content[t1]+information_content[t2])))
		}
	))
}

#' Calculate Lin similarity score of two term sets
#'
#' Warning! This function is slow - performing large numbers of `between term-set' similarity calculations should be done using \code{\link{get_sim_grid}}.
#'
#' @template ontology
#' @template information_content
#' @param term_set_1 Character vector of terms.
#' @param term_set_2 Character vector of terms.
#' @return Numeric value.
#' @seealso \code{\link{resnik}}, \code{\link{get_term_sim_mat}}
#' @references Lin D (1998). `An Information-Theoretic Definition of Similarity.' In Shavlik JW (ed.), _Proceedings of the Fifteenth International Conference on Machine Learning (ICML 1998), Madison, Wisconsin, USA, July 24-27, 1998_, pp. 296-304.
#' @export
lin <- function(
	ontology,
	information_content,
	term_set_1,
	term_set_2
) {
	if (min(c(length(term_set_1), length(term_set_2))) == 0)
		return(0)
	(lin_asym(ontology, information_content, term_set_1, term_set_2) + lin_asym(ontology, information_content, term_set_2, term_set_1))/2	
}

