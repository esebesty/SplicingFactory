#' Calculate different diversity indexes for a matrix of transcripts.
#'
#' @param x An input \code{matrix}, or \code{data.frame} containing
#'   transcript-level expression values.
#' @param genes Character vector with equal length to the number of rows of the
#'   input dataset with transcript-level expression values. The values are
#'   grouped into genes, based on this vector.
#' @param method Method to use for splicing diversity calculation, including
#'   naive entropy (\code{naive}), Laplace entropy (\code{laplace}), Gini index
#'   (\code{gini}), Simpson index (\code{simpson}) and inverse Simpson index
#'   (\code{invsimpson}). The default method is Laplace entropy.
#' @param norm If \code{TRUE}, the entropy values are normalized to the number
#'   of transcripts for each gene. The normalized entropy values are always
#'   between 0 and 1. If \code{FALSE}, genes cannot be compared to each other,
#'   due to possibly different maximum entropy values.
#' @param verbose If \code{TRUE}, the function will print additional diagnostic
#'    messages.
#' @return Gene-level splicing diversity values in a \code{matrix}, where each
#'   row belongs to a gene and each column belongs to a sample from the data.
#' @details The function aggregates basic diversity calculations to a matrix of
#' transcript-level expression values, in order to calculate gene-level naive
#' and Laplace entropy values, Gini, Simpson or inverse-Simpson diversity
#' indexes.
#' @import stats
calculate_method <- function(x, genes, method, norm = TRUE, verbose = FALSE) {
    if (method == "naive") {
        x <- aggregate(x, by = list(genes), calculate_entropy, norm = norm)
    }
    if (method == "laplace") {
        x <- aggregate(x, by = list(genes), calculate_entropy, norm = norm,
                       pseudocount = 1)
    }
    if (method == "gini") {
        x <- aggregate(x, by = list(genes), calculate_gini)
    }
    if (method == "simpson") {
        x <- aggregate(x, by = list(genes), calculate_simpson)
    }
    if (method == "invsimpson") {
        x <- aggregate(x, by = list(genes), calculate_inverse_simpson)
    }

    y <- x[apply(x[2:ncol(x)], 1, function(X) all(!is.nan(X))), ]

    if (nrow(x) - nrow(y) > 0 && verbose == TRUE) {
        message(paste0("Note: There are ", nrow(x) - nrow(y), " genes with single isoforms,
    which will be exluded from the analysis."))
    }

    colnames(y)[1] <- "Gene"

    return(y)
}
