allowed_class_type_map <- setNames(nm=c("matrix", "numeric", "integer", "sim_index"), c("matrix","numeric","numeric","sim_index"))

check_pop_sim_class <- function(class_name, type) {
	if (!any(class_name %in% names(allowed_class_type_map)))
		stop(paste0("Function not applicable to objects of class ", paste0("'", class_name, "'", collapse=", ")))
	if (!(allowed_class_type_map[match(class_name, names(allowed_class_type_map))[1]]==type))
		stop(paste0("Class-type mismatch"))
}

get_c_group_inds <- function(group, nm, pop_size) {
	zero_inds <- if (is.character(group)) {
		if (is.null(nm)) stop("Term sets are not named, pass integer indices instead")
		match(group, nm)-1L 
	} else {
		group-1L
	}
	if (!all(zero_inds < pop_size & zero_inds >= 0L))
		stop("Group indices must be between 1 and size of population")

	as.integer(zero_inds)
}

#' Calculate the group similarity of a set of row/column indices
#'
#' Calculates the similarity of a group within a population by applying the function specified by \code{group_sim} to the pairwise similarities of group members.
#'
#' @template pop_sim
#' @template group
#' @template type
#' @template group_sim
#' @param ... Other arguments to be passed to \code{get_sim}.
#' @return Numeric value of group similarity
#' @export
get_sim <- function(pop_sim, ...) UseMethod("get_sim")

#' @rdname get_sim
#' @method get_sim integer
#' @seealso \code{\link{get_sim_p}} \code{\link{sample_group_sim}}
#' @name get_sim
NULL

#' @export
#' @rdname get_sim
get_sim.integer <- function(pop_sim, ...) {
	get_sim.numeric(pop_sim, ...)
}

#' @export
#' @rdname get_sim
#' @method get_sim numeric
get_sim.numeric <- function(pop_sim, group=seq(length(pop_sim)), ...) {
	get_sim.default(
		group=get_c_group_inds(group, names(pop_sim), length(pop_sim)),
		pop_sim=pop_sim,
		type="numeric",
		...
	)
}

#' @export
#' @rdname get_sim
#' @method get_sim matrix
get_sim.matrix <- function(pop_sim, group=seq(nrow(pop_sim)), ...) { 
	stopifnot(nrow(pop_sim) == ncol(pop_sim))
	stopifnot(identical(rownames(pop_sim), colnames(pop_sim)))
	stopifnot(nrow(pop_sim) > 0)
	get_sim.default(
		group=get_c_group_inds(group, rownames(pop_sim), nrow(pop_sim)),
		pop_sim=pop_sim,
		type="matrix",
		...
	)
}

#' @export
#' @rdname get_sim
#' @method get_sim sim_index
get_sim.sim_index <- function(pop_sim, group=seq(pop_sim[["N"]]), ...) {
	get_sim.default(
		group=get_c_group_inds(group, pop_sim[["nm"]], pop_sim[["N"]]),
		pop_sim=pop_sim,
		type="sim_index",
		...
	)
}

#' @export
#' @rdname get_sim
#' @method get_sim default
get_sim.default <- function(pop_sim, group, type, group_sim="average", ...) {
	check_pop_sim_class(class(pop_sim), type)
	stopifnot(group_sim %in% c("average","min"))
	stopifnot(is.integer(group))

	group_sim(
		type,
		pop_sim,
		group_sim=="average",
		group
	)
}

#' Get similarity p-value
#'
#' p-value of group similarity, calculated by estimating the proportion by random sampling of groups the same size as \code{group} which have at least as great group similarity than does \code{group}.
#'
#' @template pop_sim
#' @template group
#' @template type
#' @template min_its 
#' @template max_its
#' @template signif
#' @template log_dismiss
#' @template group_sim
#' @param ... Arguments for \code{get_sim_p}.
#' @return p-value. 
#' @seealso \code{\link{get_sim}} \code{\link{sample_group_sim}}
#' @name get_sim_p
NULL

#' @export
#' @rdname get_sim_p
get_sim_p <- function(pop_sim, ...) UseMethod("get_sim_p")

#' @export
#' @rdname get_sim_p
#' @method get_sim_p integer
get_sim_p.integer <- function(pop_sim, ...) {
	get_sim_p.numeric(pop_sim, ...)
}

#' @export
#' @rdname get_sim_p
#' @method get_sim_p numeric
get_sim_p.numeric <- function(pop_sim, group, ...) {
	get_sim_p.default(
		group=get_c_group_inds(group, names(pop_sim), length(pop_sim)),
		pop_sim=pop_sim,
		type="numeric",
		...
	)
}

#' @export
#' @rdname get_sim_p
#' @method get_sim_p matrix
get_sim_p.matrix <- function(pop_sim, group, ...) { 
	stopifnot(nrow(pop_sim) == ncol(pop_sim))
	stopifnot(identical(rownames(pop_sim), colnames(pop_sim)))
	stopifnot(nrow(pop_sim) > 0)
	get_sim_p.default(
		group=get_c_group_inds(group, rownames(pop_sim), nrow(pop_sim)),
		pop_sim=pop_sim,
		type="matrix",
		...
	)
}

#' @export
#' @rdname get_sim_p
#' @method get_sim_p sim_index
get_sim_p.sim_index <- function(pop_sim, group, ...) {
	get_sim_p.default(
		group=get_c_group_inds(group, pop_sim[["nm"]], pop_sim[["N"]]),
		pop_sim=pop_sim,
		type="sim_index",
		...
	)
}

#' @export
#' @rdname get_sim_p
#' @method get_sim_p default
#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
get_sim_p.default <- function(
	pop_sim, 
	group, 
	type, 
	min_its=1e3, 
	max_its=1e5, 
	signif=0.05, 
	log_dismiss=log(1e-6), 
	group_sim="average",
	...
) {
	check_pop_sim_class(class(pop_sim), type)
	stopifnot(group_sim %in% c("average","min"))
	stopifnot(is.integer(group))
	if (length(group) < 2)
		stop("'group' should contain at least two members")

	sim_p(
		type,
		pop_sim,
		group_sim=="average",
		group,
		min_its,
		max_its,
		signif,
		log_dismiss
	)
}

#' Get similarity p-value for subgroup of list of term sets
#'
#' @template ontology
#' @template term_sets
#' @template information_content
#' @template term_sim_method
#' @template combine
#' @param ... Other arguments to be passed to \code{\link{get_sim_p}}.
#' @return Numeric value.
#' @export
#' @seealso \code{\link{get_sim_p}} \code{\link{create_sim_index}}
get_sim_p_from_ontology <- function(
	ontology,
	term_sets,
	information_content=descendants_IC(ontology),
	term_sim_method="lin",
	combine="average",
	...
) {
	get_sim_p.sim_index(
		create_sim_index(
			ontology=ontology, 
			term_sets=term_sets, 
			information_content=information_content, 
			term_sim_method="lin",
			combine="average"
		),
		...
	)
}

#' Draw sample of group similarities of groups of given size
#'
#' @template pop_sim
#' @param group_size Integer giving the number of members of a group.
#' @param sample_size Number of samples to draw. 
#' @template type
#' @template group_sim
#' @param ... Other arguments to be passed to \code{sample_group_sim}.
#' @return Numeric vector of random group similarities.
#' @seealso \code{\link{get_sim}} \code{\link{get_sim_p}}
#' @name sample_group_sim
NULL

#' @export
#' @rdname sample_group_sim
sample_group_sim <- function(pop_sim, ...) UseMethod("sample_group_sim")

#' @export
#' @rdname sample_group_sim
#' @method sample_group_sim integer
sample_group_sim.integer <- function(pop_sim, ...) {
	sample_group_sim.numeric(pop_sim, ...)
}

#' @export
#' @rdname sample_group_sim
#' @method sample_group_sim numeric
sample_group_sim.numeric <- function(pop_sim, ...) {
	sample_group_sim.default(
		pop_sim=pop_sim,
		type="numeric",
		...
	)
}

#' @export
#' @rdname sample_group_sim
#' @method sample_group_sim matrix
sample_group_sim.matrix <- function(pop_sim, ...) { 
	stopifnot(nrow(pop_sim) == ncol(pop_sim))
	stopifnot(identical(rownames(pop_sim), colnames(pop_sim)))
	stopifnot(nrow(pop_sim) > 0)
	sample_group_sim.default(
		pop_sim=pop_sim,
		type="matrix",
		...
	)
}

#' @export
#' @rdname sample_group_sim
#' @method sample_group_sim sim_index
sample_group_sim.sim_index <- function(pop_sim, ...) {
	sample_group_sim.default(
		pop_sim=pop_sim,
		type="sim_index",
		...
	)
}

#' @export
#' @importFrom Rcpp evalCpp
#' @useDynLib ontologySimilarity
#' @rdname sample_group_sim
sample_group_sim.default <- function(pop_sim, type, group_size, group_sim="average", sample_size=1e4, ...) {
	check_pop_sim_class(class(pop_sim), type)
	stopifnot(group_sim %in% c("average","min"))

	sample_null(
		type,
		pop_sim,
		group_sim=="average",
		group_size,
		sample_size
	)
}

#' Draw sample of group similarities for groups of given size based on \code{ontology} argument
#'
#' @template ontology
#' @template term_sets
#' @template information_content
#' @template term_sim_method
#' @template combine
#' @param ... Other arguments to be passed to \code{\link{get_sim_p}}.
#' @return Numeric vector of group similarities.
#' @export
#' @seealso \code{\link{sample_group_sim}} \code{\link{create_sim_index}}
sample_group_sim_from_ontology <- function(
	ontology,
	term_sets,
	information_content=descendants_IC(ontology),
	term_sim_method="lin",
	combine="average",
	...
) {
	sample_group_sim.sim_index(
		create_sim_index(
			ontology=ontology, 
			term_sets=term_sets, 
			information_content=information_content, 
			term_sim_method="lin",
			combine="average"
		),
		...
	)
}

