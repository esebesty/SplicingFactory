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
#' @param q Tsallis entropy parameter (q > 0). Only used if method = "tsallis". Default is 2. If q = 1, the function returns the Shannon entropy.
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
        # It is not possible to use aggregate here, because it expect to return
        # a single value per group, but calculate_tsallis_entropy can return
        # multiple values (if length(q) > 1)
        gene_levels <- unique(genes)
        coln <- as.vector(outer(colnames(x), q, function(s, qq) paste0(s, "_q=", qq)))
        rown <- gene_levels
        tsallis_row <- function(gene) {
            idx <- which(genes == gene)
            unlist(lapply(seq_len(ncol(x)), function(j) {
                v <- calculate_tsallis_entropy(x[idx, j], q = q, norm = norm)
                if (length(v) == length(q) && all(is.finite(v) | is.na(v))) v else setNames(rep(NA_real_, length(q)), paste0("q=", q))
            }))
        }
        result_mat <- t(vapply(gene_levels, tsallis_row, FUN.VALUE = setNames(numeric(length(coln)), coln)))
        colnames(result_mat) <- coln
        rownames(result_mat) <- rown
        out_df <- data.frame(Gene = rown, result_mat, check.names = FALSE)
        if (all(rowSums(!is.na(result_mat)) == 0)) {
            out_df <- data.frame(Gene=character(0))
            for (nm in coln) out_df[[nm]] <- numeric(0)
            x <- out_df
            return(x)
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
