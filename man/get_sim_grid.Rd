% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/term-set-sim.R
\name{get_sim_grid}
\alias{get_sim_grid}
\title{Get similarity matrix of pairwise similarities of term sets.}
\usage{
get_sim_grid(
  ontology,
  information_content,
  term_sim_method,
  term_sim_mat,
  term_sets,
  term_sets2 = term_sets,
  combine = "average"
)
}
\arguments{
\item{ontology}{\code{ontology_index} object.}

\item{information_content}{Numeric vector of information contents of terms (named by term)}

\item{term_sim_method}{Character string equalling either "lin" or "resnik" to use Lin or Resnik's expression for the similarity of terms.}

\item{term_sim_mat}{Numeric matrix with rows and columns corresponding to (and named by) term IDs, and cells containing the similarity between the row and column term}

\item{term_sets}{List of character vectors of ontological term IDs.}

\item{term_sets2}{Second set of term sets.}

\item{combine}{Character string - either "average" or "product", indicating whether to use the best-match-product' method, or function accepting two arguments - the first, the similarity matrix obtained by averaging across term sets in \code{term_sets}, and the second averaging across those in \code{term_sets2}.}
}
\value{
Numeric matrix of pairwise term set similarities.
}
\description{
Using either an \code{ontology_index} object and numeric vector of information content per term - or a matrix of between-term similarities (e.g. the output of \code{\link{get_term_sim_mat}}), create a numeric matrix of `between-term set' similarities. Either the `best-match-average' or `best-match-product' approach (i.e. where the 2 scores obtained by applying the asymmetric `best-match' similarity function to two term sets in each order are combined by taking the average or the product respectively). Either Lin's (default) or Resnik's definition of term similarity can be used. If \code{information_content} is not specified, a default value from \code{\link{descendants_IC}} is generated.
}
\details{
Note that if any term set within \code{term_sets} has 0 terms associated with it, it will get a similarity of 0 to any other set. If you do not want to compare term sets with no annotation, take care to filter out empty sets first, e.g. by `term_sets=term_sets[sapply(term_sets, length) > 0]`.
}
\examples{
library(ontologyIndex)
data(hpo)
term_sets <- list(
`case1`=c("HP:0001873", "HP:0011877"),
`case2`=c("HP:0001872", "HP:0001892"),
`case3`="HP:0001873")
get_sim_grid(ontology=hpo, term_sets=term_sets)
}
\seealso{
\code{\link{get_term_sim_mat}} \code{\link{get_sim_p}} \code{\link{get_asym_sim_grid}}
}
