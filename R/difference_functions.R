#' Compare Tsallis entropy across multiple q values between groups
#'
#' @param tsallis_matrix Matrix/data.frame: rows = genes/samples, columns = q values (e.g., q=1.1, q=1.5, ...)
#' @param group Vector: group assignment for each sample
#' @param paired Logical: whether samples are paired (default FALSE)
#' @param blocks Optional: subject IDs for pairing (required if paired = TRUE)
#' @param pcorr P-value correction method
#' @return Raw and corrected p-values for difference across q values between groups
#' @import stats
#' Calculate splicing diversity changes between two conditions.
#'
#' @param x A \code{matrix} with the splicing diversity values.
#' @param samples Character vector with an equal length to the number of columns
#'   in the input dataset, specifying the category of each sample.
#' @param control Name of the control sample category, defined in the
#'   \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control =
#'   'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#'   value in a condition. Can be \code{'mean'} or \code{'median'}.
#' @return A \code{data.frame} with mean or median value of splicing diversity
#'   across sample categories, the difference between these values and the log2
#'   fold change values.
#' @details The function uses a matrix of splicing diversity values in order to
#' calculate mean or median differences and log2 fold changes between two
#' conditions.
#' @import stats
calculate_fc <- function(x, samples, control, method = "mean") {
    if (method == "mean") {
        value <- aggregate(t(x), by = list(samples), mean, na.rm = TRUE)
    }

    if (method == "median") {
        value <- aggregate(t(x), by = list(samples), median, na.rm = TRUE)
    }

    sorted <- value[value$Group.1 != control, ]
    sorted[2, ] <- value[value$Group.1 == control, ]
    value <- t(sorted[, -1])
    value[is.nan(value[, 1]), c(1)] <- NA
    value[is.nan(value[, 2]), c(2)] <- NA

    result <- data.frame(value, ifelse(as.matrix(is.na(value[, 1]) | is.na(value[, 2])), NA, value[, 1] - value[,
        2]), ifelse(as.matrix(is.na(value[, 1]) | is.na(value[, 2])), NA, as.matrix(log(value[, 1] / value[, 2],
        base = 2))))
    colnames(result) <- c(paste0(sorted[1, 1], "_", method), paste0(sorted[2, 1], "_", method), paste0(method,
        "_difference"), "log2_fold_change")
    return(result)
}

#' Calculate p-values using Wilcoxon rank sum test.
#'
#' @param x A \code{matrix} with the splicing diversity values.
#' @param samples Character vector with an equal length to the number of columns
#' in the input dataset, specifying the category of each sample.
#' @param pcorr P-value correction method applied to the results, as defined in
#' the \code{p.adjust} function.
#' @param paired If \code{TRUE}, the Wilcox-test will be paired, and therefore
#' it will be a signed rank test instead of the rank sum test.
#' @param exact If \code{TRUE}, an exact p-value will be computed.
#' @return Raw and corrected p-values in a matrix.
#' @import stats
wilcoxon <- function(x, samples, pcorr = "BH", paired = FALSE, exact = FALSE) {
    p_values <- vector("list", nrow(x))

    for (i in seq_len(nrow(x))) {
        p_values[i] = wilcox.test(x[i, as.numeric(which(samples %in% unique(sort(samples))[1]))], x[i, as.numeric(which(samples %in%
            unique(sort(samples))[2]))], paired = paired, exact = exact)$p.value
    }

    raw_p_values <- ifelse(is.nan(vapply(p_values, c, numeric(1))), 1, vapply(p_values, c, numeric(1)))
    adjusted_p_values <- p.adjust(raw_p_values, method = pcorr)
    return(cbind(raw_p_values, adjusted_p_values))
}

#' Calculate p-values using label shuffling.
#'
#' @param x A \code{matrix} with the splicing diversity values.
#' @param samples Character vector with an equal length to the number of columns
#' in the input dataset, specifying the category of each sample.
#' @param control Name of the control sample category, defined in the
#' \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control = 'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#' value in a condition. Can be \code{'mean'} or \code{'median'}.
#' @param randomizations The number of random shuffles.
#' @param pcorr P-value correction method applied to the results, as defined in
#' the \code{p.adjust} function.
#' @return Raw and corrected p-values.
#' @import stats
label_shuffling <- function(x, samples, control, method, randomizations = 100,
                            pcorr = "BH") {

    log2_fc <- calculate_fc(x, samples, control, method)[, 4]

    shuffled <- replicate(randomizations, calculate_fc(x, sample(samples), control, method))
    shuffled <- vapply(shuffled, c, numeric(length(log2_fc)))[, seq(4, length(shuffled), by = 4)]

    p_values <- vector("list", nrow(x))

    for (i in seq_len(nrow(shuffled))) {
        if (ecdf(shuffled[i, ])(log2_fc[i]) >= 0.5) {
            p_values[i] <- 1 - ecdf(shuffled[i, ])(log2_fc[i])
        }
        if (ecdf(shuffled[i, ])(log2_fc[i]) < 0.5) {
            p_values[i] <- ecdf(shuffled[i, ])(log2_fc[i])
        }
        if (length(unique(shuffled[i, ])) == 1 | ecdf(shuffled[i, ])(log2_fc[i]) == 1) {
            p_values[i] <- 1
        }
        if (ecdf(shuffled[i, ])(log2_fc[i]) == 0) {
            p_values[i] <- 0 + .Machine$double.eps
        }
    }

    raw_p_values <- ifelse(is.nan(vapply(p_values, c, numeric(1))), 1, vapply(p_values, c, numeric(1)))
    adjusted_p_values <- p.adjust(raw_p_values, method = pcorr)
    return(cbind(raw_p_values, adjusted_p_values))
}

#' Compare Tsallis entropy across multiple q values between groups
#'
#' @param tsallis_matrix Matrix/data.frame: rows = genes/samples, columns = q values (e.g., q=1.1, q=1.5, ...)
#' @param group Vector: group assignment for each sample
#' @param paired Logical: whether samples are paired (default FALSE)
#' @param blocks Optional: subject IDs for pairing (required if paired = TRUE)
#' @param pcorr P-value correction method
#' @return Raw and corrected p-values for difference across q values between groups
#' @import stats
compare_q_values <- function(tsallis_matrix, group, paired = FALSE, blocks = NULL, pcorr = "BH") {
    # tsallis_matrix: samples x q values
    if (!paired) {
        # For each q, compare between groups (Wilcoxon)
        pvals <- sapply(seq_len(ncol(tsallis_matrix)), function(j) {
            wilcox.test(tsallis_matrix[group == unique(group)[1], j], tsallis_matrix[group == unique(group)[2], j])$p.value
        })
        adj_pvals <- p.adjust(pvals, method = pcorr)
        return(data.frame(q = colnames(tsallis_matrix), raw_p = pvals, adj_p = adj_pvals))
    } else {
        # For paired/repeated measures, use Friedman test across q
        if (is.null(blocks)) stop("blocks (subject IDs) required for paired analysis.")
        # Reshape to long format: value, q, group, block
        df <- data.frame(value = as.vector(tsallis_matrix),
                        q = rep(colnames(tsallis_matrix), each = nrow(tsallis_matrix)),
                        group = rep(group, times = ncol(tsallis_matrix)),
                        block = rep(blocks, times = ncol(tsallis_matrix)))
        # For each group, test across q
        results <- lapply(unique(group), function(g) {
            subdf <- df[df$group == g, ]
            p <- friedman.test(value ~ q | block, data = subdf)$p.value
            data.frame(group = g, raw_p = p, adj_p = p)
        })
        return(do.call(rbind, results))
    }
}
