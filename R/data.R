#' \code{list} object containing character vectors of term IDs of GO terms annotating each gene, named by gene. Users can simply select a list of annotation vectors for a subset of the annotated genes using a character vector of gene symbols, e.g. \code{gene_GO_terms[c("ACTN1", "TUBB1")]}, which can then be used in functions for calculating similarities, e.g. \code{\link{get_sim_grid}}. Note that these annotation vectors contain annotation from all major branches of the Gene Ontology, however one can simply extract the terms only relevant to one by calling the function in the \code{ontologyIndex} package: \code{intersection_with_branches}. 
#' 
#' @name gene_GO_terms 
#' @title Gene Ontology annotation of genes
#' @docType data
#' @format List of character vectors.
#' @references Annotation downloaded from Gene Ontology consortium website, http://geneontology.org/, dated 7/7/2016.
NULL

#' Numeric vector containing the information content of Gene Ontology terms based on frequencies of annotation data object \code{gene_GO_terms}. The object can be derived using the function \code{get_term_info_content} and data object \code{go} from the \code{ontologyIndex} package.
#' 
#' @name GO_IC 
#' @title Gene Ontology terms information content.
#' @docType data
#' @format List of character vectors.
NULL

