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

#' @title Calculate the group similarity of a set of row/column indices
#' @description Calculates the similarity of a group within a population by applying the function specified by \code{group_sim} to the pairwise similarities of group members.
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

#' @title Get similarity p-value
#' @description p-value of group similarity, calculated by estimating the proportion by random sampling of groups the same size as \code{group} which have at least as great group similarity than does \code{group}.
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

#' @title Get similarity p-value for subgroup term sets
#' @description Compute a similarity p-value by permutation for subgroup of a list of term sets
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

#' @title Draw sample of group similarities
#' @description Draw sample of group similarities of groups of given size
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

#' @title Draw sample of group similarities
#' @description ample of group similarities for random groups of given drawn from the given \code{ontology} argument
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

#' @title Identify enriched terms in subgroup
#' @description Create a table of terms ranked by their significance of occurrence in a set of term sets amongst an enclosing set, with p-values computed by permutation. Terms are subselected so that only the minimal set of non-redundant terms at each level of frequency within the group are retained.
#' @template ontology
#' @template term_sets
#' @param group Integer/logical/character vector specifying indices/positions/names of subgroup for which to calculate a group similarity p-value.
#' @param permutations Number of permutations to test against, or if \code{NULL}, perform no permutations and return the unadjusted p-values for the occurrence of each term.
#' @param min_terms Minimum number of times a term should occur within the given group to be eligible for inclusion in the results.
#' @return \code{data.frame} containing columns: \code{term} (with the term ID); \code{name} (term readable name); \code{in_term} (number of sets in the given group of containing the term); \code{in_no_term} (number of sets in the given group not containing the term); \code{out_term} and \code{out_no_term} (equivalently for the sets not in the given group); \code{p} (the p-values calculated by permutation for seeing a term with such a strong association, measured using Fisher's exact test, in a group of term sets the size of the given group among \code{term_sets}). Rows ordered by significance (i.e. the \code{p} columns).
#' @param mc.cores If not null and greater than on, the number of cores use calculating permutations (passed to \code{mclapply}).
#' @export
#' @seealso \code{\link{sample_group_sim}} \code{\link{create_sim_index}}
#' @importFrom parallel mclapply
#' @importFrom stats fisher.test
#' @importFrom ontologyIndex minimal_set 
#' @export
group_term_enrichment <- function(
	ontology,
	term_sets,
	group,
	permutations=1000L,
	min_terms=2L,
	mc.cores=NULL
) {
	if (is.character(group) & !is.null(names(term_sets))) group <- names(term_sets) %in% group
	if (is.integer(group)) group <- seq_along(term_sets) %in% group
	if (!is.logical(group)) {
		stop("'group' argument must be integer/logical/character vector specifying indices/positions/names of subgroup within 'term_sets'")
	}
	stopifnot(length(group)==length(term_sets))
	w_ancs <- lapply(term_sets, get_ancestors, ontology=ontology)
	all_terms <- names(which(table(as.character(unlist(use.names=FALSE, w_ancs))) >= min_terms))
	tf <- factor(as.character(unlist(use.names=FALSE, w_ancs)), levels=all_terms)
	fn <- function(group, just_p=FALSE) {
		f <- factor(as.integer(group), 0:1)
		term_tab <- (function(x) x[,x[2,]>=min_terms,drop=FALSE])(as.matrix(table(
			gt=rep(f, times=lengths(w_ancs)),
			term=tf
		)))
		ft <- table(f)
		ps <- apply(term_tab, 2, function(x) fisher.test(rbind(x, ft-x))$p.value)
		if (just_p) min(ps) else list(ft=ft, p=ps, terms=term_tab)
	}
	r <- fn(group)
	m <- structure(as.matrix(r$terms), class="matrix")
	v <- t(rbind(m, as.integer(r$ft)-m))
	colnames(v) <- c("out_term","in_term","out_no_term","in_no_term")
	rep_fn <- function(x) replicate(x, fn(sample(group), TRUE), TRUE)
	df <- data.frame(
		stringsAsFactors=FALSE,
		check.names=FALSE,
		row.names=NULL,
		term=rownames(v),
		name=ontology$name[rownames(v)],
		v[,c(2,4,1,3)],
		p=if (is.null(permutations)) r$p else local({ ps <- sort(c(min(r$p), if (is.null(mc.cores) || mc.cores == 1) rep_fn(permutations) else { as.numeric(unlist(use.names=FALSE,mclapply(mc.cores=mc.cores,table(gl(n=mc.cores,k=1,length=permutations)),function(i) rep_fn(i)))) })); dedup <- ps[!duplicated(ps)]; cumsum(table(cut(x=ps, right=FALSE, labels=FALSE, breaks=c(dedup, Inf))))[cut(x=r$p, right=FALSE, labels=FALSE, breaks=c(dedup, Inf))]/(permutations+1) })
	)[order(r$p),]
	df[unsplit(f=df$`in_term`, value=lapply(split(df$term, f=df$`in_term`), function(x) x %in% minimal_set(ontology, terms=x))),]
}
