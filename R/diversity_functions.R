#' Calculate entropy for a vector of transcript-level
#' expression values of one gene.
#'
#' @param x Vector of expression values.
#' @param norm If \code{TRUE}, the entropy values are normalized to the number
#' of transcripts for each gene. The normalized entropy values are always
#' between 0 and 1. If \code{FALSE}, genes cannot be compared to each other,
#' due to possibly different maximum entropy values.
#' @param pseudocount Pseudocount added to each transcript expression value.
#' Default is 0, while Laplace entropy uses a pseudocount of 1.
#' @export
#' @return A single gene-level entropy value.
#' @details
#' The function calculates an entropy value as part of different
#' diversity calculations. Given a vector of transcript-level expression values
#' of a gene, this function characterizes the diversity of splicing isoforms for
#' a gene. If there only a single transcript, the diversity value will be NaN,
#' as it cannot be calculated. If the expression of the given gene is 0,
#' the diversity value will be NA.
#' @examples
#' # read counts for the transcripts of a single gene with 5 transcripts
#' x <- rnbinom(5, size = 10, prob = 0.4)
#' # calculate non-normalized naive entropy value
#' entropy <- calculate_entropy(x, norm = FALSE)
#' # calculate Laplace-entropy, also normalized for transcript number
#' # (the default)
#' norm_laplace_entropy <- calculate_entropy(x, pseudocount = 1)
calculate_entropy <- function(x, norm = TRUE, pseudocount = 0) {
    if (sum(x) != 0 & length(x) > 1) {
        x <- (x + pseudocount) / sum(x + pseudocount)
        x_log = ifelse(is.finite(log(x, base = 2)), log(x, base = 2), 0)

        if (norm == FALSE) {
            x = -sum(x * x_log)
        }
        if (norm == TRUE) {
            x = -sum(x * x_log) / log2(length(x))
        }
    } else if (length(x) == 1) {
        x = NaN
    } else {
        x = NA
    }
    return(x)
}

#' Calculate Gini coefficient for a vector of transcript-level
#' expression values of one gene.
#'
#' @param x Vector of expression values.
#' @export
#' @return A single gene-level Gini coefficient.
#' @details
#' The function calculates a Gini coefficient as part of different
#' diversity calculations. Given a vector of transcript-level expression values
#' of a gene, this function characterize the diversity of splicing isoforms for
#' a gene. If there only one single transcript, the resulted index will be NaN,
#' as diversity cannot be calculated. If the expression of the given gene is 0,
#' the diversity index will be NA.
#' @examples
#' # read counts for the transcripts of a single gene with 5 transcripts
#' x <- rnbinom(5, size = 10, prob = 0.4)
#' # calculate Gini index
#' gini <- calculate_gini(x)
calculate_gini <- function(x) {
    if (sum(x) != 0 & length(x) > 1) {
        x <- sort(x)
        y <- 2 * sum(x * seq_len(length(x))) / sum(x) - (length(x) + 1L)
        y <- y / (length(x) - 1L)
    } else if (length(x) == 1) {
        y = NaN
    } else {
        y = NA
    }
    return(y)
}

#' #' Calculate Simpson index for a vector of transcript-level
#' expression values of one gene.
#'
#' @param x Vector of expression values.
#' @export
#' @return A single gene-level Simpson index.
#' @details
#' The function calculates a Simpson index as part of different
#' diversity calculations. Given a vector of transcript-level expression values
#' of a gene, this function characterize the diversity of splicing isoforms for
#' a gene. If there only one single transcript, the resulted index will be NaN,
#' as diversity cannot be calculated. If the expression of the given gene is 0,
#' the diversity index will be NA.
#' @examples
#' # read counts for the transcripts of a single gene with 5 transcripts
#' x <- rnbinom(5, size = 10, prob = 0.4)
#' # calculate Simpson index
#' simpson <- calculate_simpson(x)
calculate_simpson <- function(x) {
    if (sum(x) != 0 & length(x) > 1) {
        x <- x / sum(x)
        x <- 1 - sum(x * x)
    } else if (length(x) == 1) {
        x = NaN
    } else {
        x = NA
    }
    return(x)
}

#' #' Calculate inverse Simpson index for a vector of transcript-level
#' expression values of one gene.
#'
#' @param x Vector of expression values.
#' @export
#' @return A single gene-level inverse Simpson index.
#' @details
#' The function calculates an inverse Simpson index as part of different
#' diversity calculations. Given a vector of transcript-level expression values
#' of a gene, this function characterize the diversity of splicing isoforms for
#' a gene. If there only one single transcript, the resulted index will be NaN,
#' as diversity cannot be calculated. If the expression of the given gene is 0,
#' the diversity index will be NA.
#' @examples
#' # read counts for the transcripts of a single gene with 5 transcripts
#' x <- rnbinom(5, size = 10, prob = 0.4)
#' # calculate inverse Simpson index
#' invsimpson <- calculate_inverse_simpson(x)
calculate_inverse_simpson <- function(x) {
    if (sum(x) != 0 & length(x) > 1) {
        x <- x / sum(x)
        x <- 1 / sum(x * x)
    } else if (length(x) == 1) {
        x = NaN
    } else {
        x = NA
    }
    return(x)
}

#' Calculate Tsallis entropy for a vector of transcript-level
#' expression values of one gene.
#'
#' @param x Vector of expression values.
#' @param q Tsallis entropy parameter (q > 0). Default is 2. If q = 1, the function returns the Shannon entropy.
#' @param norm If \code{TRUE}, the entropy values are normalized to the number
#' of transcripts for each gene. The normalized entropy values are always
#' between 0 and 1. If \code{FALSE}, genes cannot be compared to each other,
#' due to possibly different maximum entropy values.
#' @export
#' @return A single gene-level Tsallis entropy value.
#' @details
#' The function calculates the Tsallis entropy, a generalization of Shannon entropy.
#' For q → 1, Tsallis entropy converges to Shannon entropy.
#' If there is only a single transcript, the value will be NaN.
#' If the expression of the given gene is 0, the value will be NA.
##' @examples
##' x <- rnbinom(5, size = 10, prob = 0.4)
##' tsallis <- calculate_tsallis_entropy(x, q = 2)
##'
##' @param q Tsallis entropy parameter (q > 0). Can be a single value or a numeric vector. Default is 2. If q = 1, the function returns the Shannon entropy.
calculate_tsallis_entropy <- function(x, q = 2, norm = TRUE) {
    if (!is.numeric(q)) stop("q must be numeric.")
    if (any(q <= 0)) stop("q must be greater than 0.")
    if (sum(x) != 0 & length(x) > 1) {
        p <- x / sum(x)
        tsallis_vec <- vapply(q, function(qi) {
            if (abs(qi - 1) < .Machine$double.eps^0.5) {
                # q == 1, return Shannon entropy
                shannon <- -sum(ifelse(p > 0, p * log2(p), 0))
                if (norm) {
                    shannon <- shannon / log2(length(x))
                }
                return(shannon)
            } else {
                tsallis <- (1 - sum(p^qi)) / (qi - 1)
                if (norm) {
                    max_tsallis <- (1 - length(x)^(1 - qi)) / (qi - 1)
                    tsallis <- tsallis / max_tsallis
                }
                return(tsallis)
            }
        }, numeric(1))
        if (length(q) == 1) {
            return(unname(tsallis_vec))
        } else {
            names(tsallis_vec) <- paste0("q=", q)
            return(tsallis_vec)
        }
    } else if (length(x) == 1) {
        return(rep(NaN, length(q)))
    } else {
        return(rep(NA, length(q)))
    }
}
