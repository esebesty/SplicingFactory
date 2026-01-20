#' Calculate diversity values for a matrix of transcripts.
#'
#' @param x An input \code{matrix}, or \code{data.frame} containing
#'   transcript-level expression values.
#' @param genes Character vector with equal length to the number of rows of the
#'   input dataset with transcript-level expression values. The values in
#'   \code{x} are grouped into genes based on this vector.
#' @param method Method to use for splicing diversity calculation, including
#'   naive entropy (\code{naive}), Laplace entropy (\code{laplace}), Tsallis entropy (\code{tsallis}),
#'   Gini index (\code{gini}), Simpson index (\code{simpson}) and inverse Simpson index
#'   (\code{invsimpson}). The default method is Laplace entropy.
#' @param norm If \code{TRUE}, the entropy values are normalized to the number
#'   of transcripts for each gene. The normalized entropy values are always
#'   between 0 and 1. If \code{FALSE}, genes cannot be compared to each other,
#'   due to possibly different maximum entropy values.
#' @param q Tsallis entropy parameter (q > 0, q != 1). Only used if method = "tsallis". Default is 2.
#' @param verbose If \code{TRUE}, the function will print additional diagnostic
#'   messages, besides the warnings and errors.
#' @return Gene-level splicing diversity values in a \code{data.frame}, where
#'   each row belongs to a gene and each column belongs to a sample from the
#'   data, in addition to the first column, containing gene names, given in the
#'   `genes` parameter.
#' @details The function calculates diversity values on a matrix of
#' transcript-level expression values, aggregated by the genes defined in the
#' \code{genes} parameter.
#' @import stats
calculate_method <- function(x, genes, method, norm = TRUE, verbose = FALSE, q = 2) {

    if (method == "naive") {
        x <- aggregate(x, by = list(genes), calculate_entropy, norm = norm)
    }

    if (method == "laplace") {
        x <- aggregate(x, by = list(genes), calculate_entropy, norm = norm,
                       pseudocount = 1)
    }

    if (method == "tsallis") {
        # For q as a vector, aggregate returns a matrix for each gene/sample
        x_list <- aggregate(x, by = list(genes), calculate_tsallis_entropy, norm = norm, q = q)
        # x_list is a data.frame, but the value columns are lists of named vectors
        # Convert to a data.frame with columns for each q
        qnames <- paste0("q=", q)
        # For each sample (column), extract a matrix of values for all genes and all q
        gene_names <- x_list[[1]]
        value_cols <- x_list[-1]
        # For each sample, build a matrix genes x q
        out_list <- lapply(value_cols, function(col) do.call(rbind, col))
        # If only one q, keep as vector
        if (length(q) == 1) {
            for (i in seq_along(out_list)) {
                out_list[[i]] <- as.vector(out_list[[i]])
                names(out_list[[i]]) <- NULL
            }
            out_df <- data.frame(Gene = gene_names, do.call(cbind, out_list), check.names = FALSE)
        } else {
            # For multiple q, flatten to wide format: one row per gene, columns for each sample/q
            out_df <- data.frame(Gene = gene_names)
            for (s in seq_along(out_list)) {
                mat <- out_list[[s]]
                colnames(mat) <- paste0(colnames(x)[s], "_", qnames)
                out_df <- cbind(out_df, mat)
            }
        }
        x <- out_df
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
