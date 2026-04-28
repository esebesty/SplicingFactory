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

#' Internal: Core Tsallis entropy calculation
#'
#' Centralized entropy calculation with support for Shannon (q=1), Tsallis (q≠1),
#' and species richness (q=0). Uses natural logarithm for mathematical rigor.
#'
#' @param proportions Numeric vector of species proportions (must sum to ~1)
#' @param q Numeric. Generalization parameter. Default: 1.0 (Shannon entropy)
#' @param norm Logical. Normalize by maximum entropy. Default: FALSE
#' @param log_base Numeric. Logarithm base. Default: exp(1) (natural log)
#' @param q_tol Numeric. Tolerance for detecting q=1 case. Default: 1e-6
#'
#' @return Numeric. Entropy value
#' @noRd
.entropy_core <- function(proportions, q = 1, norm = FALSE, log_base = exp(1), q_tol = 1e-6) {
    # Input validation
    if (!is.numeric(proportions) || length(proportions) == 0) {
        return(NA_real_)
    }

    if (any(is.na(proportions)) || any(is.infinite(proportions))) {
        return(NA_real_)
    }

    if (!is.numeric(q)) {
        stop("q must be numeric")
    }
    if (length(q) > 1) {
        stop("q must be a scalar (length 1), not a vector")
    }
    if (q < 0) {
        stop("q must be non-negative")
    }

    # Filter out zeros (standard in entropy)
    p_nonzero <- proportions[proportions > 1e-15]

    if (length(p_nonzero) == 0) {
        return(NA_real_)
    }

    # Normalize to sum to 1 (handle numerical errors)
    p <- p_nonzero / sum(p_nonzero)

    # Species richness (q=0): just count species
    if (q < q_tol) {
        H <- (log(length(p)) - 1) / log(log_base)
        return(H)
    }

    # Shannon entropy (q=1): use -sum(p*log(p))
    if (abs(q - 1) < q_tol) {
        H <- -sum(p * log(p)) / log(log_base)
    } else {
        # Tsallis entropy: (1 - sum(p^q)) / (q-1)
        # Note: Unlike Shannon, log_base is NOT applied to Tsallis
        H <- (1 - sum(p^q)) / (q - 1)
    }

    # Normalize by maximum entropy if requested
    if (norm) {
        n <- length(p)
        if (q < q_tol) {
            H_max <- (log(n) - 1) / log(log_base)
        } else if (abs(q - 1) < q_tol) {
            H_max <- log(n) / log(log_base)
        } else {
            # Tsallis: max entropy without log_base (unlike Shannon)
            H_max <- (1 - n^(1 - q)) / (q - 1)
        }

        if (!is.na(H_max) && !is.nan(H_max) && H_max > 0 && is.finite(H_max)) {
            H <- H / H_max
        }
    }

    return(H)
}

#' Internal: Vectorized Tsallis entropy calculation
#'
#' Compute entropy for multiple observations (rows = observations, cols = species)
#'
#' @param counts Matrix/data.frame where rows are observations, columns are species
#' @param q Numeric. Generalization parameter. Default: 1.0
#' @param norm Logical. Normalize by maximum entropy. Default: FALSE
#' @param log_base Numeric. Logarithm base. Default: exp(1)
#' @param pseudocount Numeric. Add to counts before normalization. Default: 0
#' @param q_tol Numeric. Tolerance for detecting q=1 case. Default: 1e-6
#'
#' @return Numeric vector of entropy values (one per observation)
#' @noRd
.entropy_vectorized <- function(counts, q = 1, norm = FALSE, log_base = exp(1), 
                                pseudocount = 0, q_tol = 1e-6) {
    counts <- as.matrix(counts)

    if (nrow(counts) == 0 || ncol(counts) == 0) {
        return(numeric(0))
    }

    # Apply to each row
    entropy_vals <- apply(counts, 1, function(row) {
        # Add pseudocount and normalize
        total <- sum(row, na.rm = TRUE) + length(row) * pseudocount
        if (total <= 0)
            return(NA_real_)

        p <- (row + pseudocount) / total
        .entropy_core(p, q = q, norm = norm, log_base = log_base, q_tol = q_tol)
    })

    return(unname(entropy_vals))
}

#' Internal: Maximum Tsallis entropy for n species
#'
#' Compute the theoretical maximum entropy for uniform distribution of n species
#'
#' @param n_species Integer. Number of species
#' @param q Numeric. Generalization parameter. Default: 1.0
#' @param log_base Numeric. Logarithm base. Default: exp(1)
#' @param q_tol Numeric. Tolerance for detecting q=1 case. Default: 1e-6
#'
#' @return Numeric. Maximum entropy value
#' @noRd
.entropy_max <- function(n_species, q = 1, log_base = exp(1), q_tol = 1e-6) {
    if (n_species < 1)
        return(NA_real_)

    if (q < q_tol) {
        # Species richness max: log(n)
        H_max <- (log(n_species) - 1) / log(log_base)
    } else if (abs(q - 1) < q_tol) {
        # Shannon max: log(n)
        H_max <- log(n_species) / log(log_base)
    } else {
        # Tsallis max: (1 - n^(1-q)) / (q-1) [log_base NOT applied to Tsallis]
        H_max <- (1 - n_species^(1 - q)) / (q - 1)
    }

    return(H_max)
}

#' Calculate Tsallis entropy for a vector of transcript-level
#' expression values of one gene.
#'
#' @param x Vector of expression values or matrix (rows = observations, cols = species).
#'   If matrix, entropy is calculated row-wise.
#' @param q Tsallis entropy parameter (q ≥ 0). Default is 2.
#'   Must be a single scalar value (not a vector).
#'   - q = 0: Species richness
#'   - q → 1: Shannon entropy
#'   - q ≠ 1: Tsallis entropy
#' @param norm If \code{TRUE}, the entropy values are normalized to [0,1].
#'   If \code{FALSE}, genes cannot be compared to each other.
#' @param log_base Logarithm base for Shannon entropy calculations.
#'   Default: exp(1) (natural logarithm). Not applied to Tsallis entropy.
#' @param pseudocount Pseudocount added to counts. Default: 0.
#'   Laplace entropy uses pseudocount = 1.
#' @param q_tol Tolerance for detecting q values close to special cases.
#'   Default: 1e-6.
#'
#' @export
#' @return 
#'   - If x is a vector: Numeric scalar
#'   - If x is a matrix: Numeric vector (one entropy per row)
#'
#' @details
#' The function calculates generalized Tsallis entropy with a fixed q parameter.
#' Only a single q value is accepted (for compatibility with statistical tests
#' like Wilcoxon). For multiple q values, use a loop or apply().
#' 
#' Includes:
#' - Species richness (q=0)
#' - Shannon entropy (q=1)
#' - Tsallis entropy (q≠1)
#' 
#' Uses natural logarithm by default for mathematical rigor.
#' If there is only a single species, the value will be NaN.
#' If all counts are zero, the value will be NA.
#'
#' @examples
#' # read counts for the transcripts of a single gene with 5 transcripts
#' x <- rnbinom(5, size = 10, prob = 0.4)
#' # calculate Tsallis entropy with q=2
#' tsallis <- calculate_tsallis_entropy(x, q = 2)
#' # calculate Shannon entropy (q=1)
#' shannon <- calculate_tsallis_entropy(x, q = 1)
calculate_tsallis_entropy <- function(x, q = 2, norm = TRUE, log_base = exp(1), 
                                      pseudocount = 0, q_tol = 1e-6) {
    
    # Input validation
    if (!is.numeric(q)) {
        stop("q must be numeric.")
    }
    if (length(q) != 1) {
        stop("q must be a single scalar value (length 1), not a vector. ",
             "Only one q parameter is accepted (required for statistical testing).")
    }
    if (q < 0) {
        stop("q must be non-negative.")
    }
    
    # Handle matrix input (vectorized)
    if (is.matrix(x) || is.data.frame(x)) {
        return(.entropy_vectorized(x, q = q, norm = norm, log_base = log_base, 
                                   pseudocount = pseudocount, q_tol = q_tol))
    }
    
    # Vector input
    x <- as.numeric(x)
    
    if (sum(x) != 0 & length(x) > 1) {
        # Calculate entropy for single q value
        p <- (x + pseudocount) / (sum(x) + length(x) * pseudocount)
        result <- .entropy_core(p, q = q, norm = norm, log_base = log_base, q_tol = q_tol)
        return(result)
    } else if (length(x) == 1) {
        return(NaN)
    } else {
        return(NA_real_)
    }
}
