#' Calculate naive entropy for a vector of transcript-level expression values of
#' one gene.
#'
#' @param x Vector of expression values.
#' @param norm If \code{TRUE}, the entropy values are normalized to the number
#'   of transcripts for each gene. The normalized entropy values are always
#'   between 0 and 1. If \code{FALSE}, genes cannot be compared to each other,
#'   due to possibly different maximum entropy values.
#' @return A single gene-level naive entropy value.
#' @details The function calculates a naive entropy value as part of different
#' diversity calculations. Given a vector of transcript-level expression values
#' of a gene, this function characterize the diversity of splicing isoforms for
#' a gene. If there only one single transcript, the resulted index will be NaN,
#' as diversity cannot be calculated. If the expression of the given gene is 0,
#' the diversity index will be NA.
calculate_entropy <- function(x, norm = TRUE) {
    if (sum(x) != 0 & length(x) > 1) {
        x <- x/sum(x)
        x_log = ifelse(is.finite(log(x, base = 2)), log(x, base = 2), 0)

        if (norm == FALSE) {
            x = -sum(x * x_log)
        }
        if (norm == TRUE) {
            x = -sum(x * x_log)/log2(length(x))
        }
    } else if (length(x) == 1) {
        x = NaN
    } else {
        x = NA
    }
    return(x)
}

#' Calculate Laplace entropy for a vector of transcript-level
#' expression values of one gene.
#'
#' @param x Vector of expression values.
#' @param norm If \code{TRUE}, the entropy values are normalized to the number
#' of transcripts for each gene. The normalized entropy values are always
#' between 0 and 1. If \code{FALSE}, genes cannot be compared to each other,
#' due to possibly different maximum entropy values.
#' @return A single gene-level Laplace entropy value.
#' @details
#' The function calculates a Laplace entropy value as part of different
#' diversity calculations. Given a vector of transcript-level expression values
#' of a gene, this function characterize the diversity of splicing isoforms for
#' a gene. If there only one single transcript, the resulted index will be NaN,
#' as diversity cannot be calculated. If the expression of the given gene is 0,
#' the diversity index will be NA.
calculate_laplace_entropy <- function(x, norm = TRUE) {
    if (sum(x) != 0 & length(x) > 1) {
        x <- (x + 1)/sum(x + 1)
        x_log = ifelse(is.finite(log(x, base = 2)), log(x, base = 2), 0)

        if (norm == FALSE) {
            x = -sum(x * x_log)
        }
        if (norm == TRUE) {
            x = -sum(x * x_log)/log2(length(x))
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
#' @return A single gene-level Gini coefficient.
#' @details
#' The function calculates a Gini coefficient as part of different
#' diversity calculations. Given a vector of transcript-level expression values
#' of a gene, this function characterize the diversity of splicing isoforms for
#' a gene. If there only one single transcript, the resulted index will be NaN,
#' as diversity cannot be calculated. If the expression of the given gene is 0,
#' the diversity index will be NA.
calculate_gini <- function(x) {
    if (sum(x) != 0 & length(x) > 1) {
        x <- sort(x)
        y <- 2 * sum(x * seq_len(length(x)))/sum(x) - (length(x) + 1L)
        y <- y/(length(x) - 1L)
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
#' @return A single gene-level Simpson index.
#' @details
#' The function calculates a Simpson index as part of different
#' diversity calculations. Given a vector of transcript-level expression values
#' of a gene, this function characterize the diversity of splicing isoforms for
#' a gene. If there only one single transcript, the resulted index will be NaN,
#' as diversity cannot be calculated. If the expression of the given gene is 0,
#' the diversity index will be NA.
calculate_simpson <- function(x) {
    if (sum(x) != 0 & length(x) > 1) {
        x <- x/sum(x)
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
#' @return A single gene-level inverse Simpson index.
#' @details
#' The function calculates an inverse Simpson index as part of different
#' diversity calculations. Given a vector of transcript-level expression values
#' of a gene, this function characterize the diversity of splicing isoforms for
#' a gene. If there only one single transcript, the resulted index will be NaN,
#' as diversity cannot be calculated. If the expression of the given gene is 0,
#' the diversity index will be NA.
calculate_inverse_simpson <- function(x) {
    if (sum(x) != 0 & length(x) > 1) {
        x <- x/sum(x)
        x <- 1/sum(x * x)
    } else if (length(x) == 1) {
        x = NaN
    } else {
        x = NA
    }
    return(x)
}
