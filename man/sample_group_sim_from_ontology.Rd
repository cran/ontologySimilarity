% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group-sim.R
\name{sample_group_sim_from_ontology}
\alias{sample_group_sim_from_ontology}
\title{Draw sample of group similarities}
\usage{
sample_group_sim_from_ontology(
  ontology,
  term_sets,
  information_content = descendants_IC(ontology),
  term_sim_method = "lin",
  combine = "average",
  ...
)
}
\arguments{
\item{ontology}{\code{ontology_index} object.}

\item{term_sets}{List of character vectors of ontological term IDs.}

\item{information_content}{Numeric vector of information contents of terms (named by term)}

\item{term_sim_method}{Character string equalling either "lin" or "resnik" to use Lin or Resnik's expression for the similarity of terms.}

\item{combine}{Character string - either "average" or "product", indicating whether to use the best-match-product' method, or function accepting two arguments - the first, the similarity matrix obtained by averaging across term sets in \code{term_sets}, and the second averaging across those in \code{term_sets2}.}

\item{...}{Other arguments to be passed to \code{\link{get_sim_p}}.}
}
\value{
Numeric vector of group similarities.
}
\description{
ample of group similarities for random groups of given drawn from the given \code{ontology} argument
}
\seealso{
\code{\link{sample_group_sim}} \code{\link{create_sim_index}}
}
