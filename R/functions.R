#' @title Get term-term similarity matrix 
#' @description Get matrix of pairwise similarity of individual terms based on Lin's (default) or Resnik's information content-based expression.
#' @template ontology
#' @template information_content 
#' @param method Character value equalling either "lin" or "resnik" to use Lin or Resnik's expression for similarity of terms respectively.
#' @param row_terms Character vector of term IDs to appear as rows of result matrix.
#' @param col_terms Character vector of term IDs to appear as cols of result matrix.
#' @return Numeric matrix of pairwise term similarities.
#' @seealso \code{\link{get_sim_grid}} \code{\link{resnik}}, \code{\link{lin}}
#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom ontologyIndex get_ancestors
#' @useDynLib ontologySimilarity
get_term_sim_mat <- function(
	ontology,
	information_content,
	method="lin",
	row_terms=names(information_content),
	col_terms=names(information_content)
) {
	if (!(method %in% c("lin","resnik")))
		stop("'method' argument must either be 'lin' or 'resnik'")

	lu_terms <- union(row_terms, col_terms)
	all_terms <- get_ancestors(ontology, lu_terms)

	if (length(setdiff(all_terms, names(information_content))) > 0)
		stop("Information content value is missing for some terms")

	sorted_ic <- sort(information_content[all_terms], decreasing=TRUE)
	ancs <- ontology$ancestors[lu_terms]
	anc_lens <- sapply(ancs, length)
	ancs_cpp <- unlist(use.names=FALSE, lapply(split(match(unlist(use.names=FALSE, ancs), names(sorted_ic))-1, rep(seq(length(ancs)), times=anc_lens)), sort))
	
	block_ends <- as.integer(cumsum(anc_lens))
	block_starts <- c(0L, block_ends[-length(block_ends)])

	resnik <- structure(
		calc_term_sim_mat(
			block_starts,
			block_ends,
			ancs_cpp,
			sorted_ic,
			match(row_terms, lu_terms)-1L,
			match(col_terms, lu_terms)-1L
		),
		dimnames=list(row_terms, col_terms))

	if (method == "lin") {
		raw <- 2 * resnik / outer(sorted_ic[row_terms], sorted_ic[col_terms], "+")
		matrix(ifelse(resnik == 0, 0, raw), nrow=nrow(raw), ncol=ncol(raw), dimnames=dimnames(raw))
	} else {
		resnik
	}
}

#' @title Get information content based on number of descendants each term has
#' @description Calculate information content of terms based on frequency with which it is an ancestor of other terms. Useful as a default if there is no population frequency information available as it captures the structure of the ontology.
#' @template ontology
#' @return Numeric vector of information contents named by term.
#' @export
#' @importFrom stats setNames
descendants_IC <- function(ontology) {
	tab <- table(factor(unlist(use.names=FALSE, ontology$ancestors), levels=ontology$id))
	setNames(nm=ontology$id, -log(as.integer(tab)/length(ontology$id)))
}

#' @importFrom ontologyIndex get_ancestors
components_for_calc_sim_from_ic <- function(ontology, information_content, list_of_term_set_lists) {
	if (!any(class(ontology) == "ontology_index"))
		stop("'ontology' must be an 'ontology_index'!")
	stopifnot(is.numeric(information_content))
	stopifnot(is.list(list_of_term_set_lists))

	all_terms <- union(unlist(use.names=TRUE, list_of_term_set_lists), get_ancestors(ontology, unlist(use.names=FALSE, list_of_term_set_lists)))
	if (!all(all_terms %in% ontology$id))
		stop("Term sets contain terms not present in ontology")
	if (length(setdiff(all_terms, names(information_content))) > 0)
		stop("Information content value is missing for some terms")

	sorted_ic <- sort(information_content[all_terms], decreasing=TRUE)
	ancs <- ontology$ancestors[names(sorted_ic)]
	anc_lens <- sapply(ancs, length)
	ancs_cpp <- unlist(use.names=FALSE, lapply(split(match(unlist(use.names=FALSE, ancs), names(sorted_ic))-1, rep(seq(length(ancs)), times=anc_lens)), sort))
	
	block_ends <- as.integer(cumsum(anc_lens))
	block_starts <- c(0L, block_ends[-length(block_ends)])

	list(
		start=block_starts,
		stop=block_ends,
		ancs=ancs_cpp,
		info=as.numeric(sorted_ic),
		term_sets=lapply(
			 list_of_term_set_lists,
			 function(term_sets) list(
				t=as.integer(match(unlist(use.names=FALSE, term_sets), names(sorted_ic)))-1,
				c=unlist(use.names=FALSE, mapply(SIMPLIFY=FALSE, FUN=rep, seq(from=0, length.out=length(term_sets)), sapply(term_sets, length))),
				n=length(term_sets)
			)
		)
	)
}

#' @title Create similarity index for list of term sets
#' @description Create light-weight similarity index for fast lookups of between term set similarity.
#' @template ontology
#' @template term_sets
#' @template information_content
#' @template term_sim_method
#' @template combine
#' @return Object of class \code{sim_index}.
#' @seealso \code{link{get_sim}} \code{\link{get_sim_p}} \code{\link{sample_group_sim}}
#' @export
create_sim_index <- function(ontology, term_sets, information_content=descendants_IC(ontology), term_sim_method="lin", combine="average") {
	if (missing(term_sets) | missing(ontology)) stop("Arguments 'ontology' and 'term_sets' must be specified")

	comps <- components_for_calc_sim_from_ic(ontology=ontology, information_content=information_content, list_of_term_set_lists=list(term_sets))

	if (!any(combine==c("average", "product"))) {
		stop("'combine' argument must be either \"average\" or \"product\"")
	}
	if (!any(term_sim_method==c("lin", "resnik"))) {
		stop("'term_sim_method' argument must be either \"lin\" or \"resnik\"")
	}

	structure(class="sim_index", c(
	  	list(N=length(term_sets), nm=names(term_sets), lin=term_sim_method=="lin", average_not_product=combine=="average"),
		comps[-which(names(comps)=="term_sets")],
		comps[["term_sets"]][[1]]
	))
}
