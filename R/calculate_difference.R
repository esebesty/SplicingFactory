#' Calculate splicing diversity changes between two conditions.
#'
#' @param x A \code{data.frame} with gene names as the first column and splicing
#'   diversity values for each sample in additional columns.
#' @param samples Character vector with an equal length to the number of columns
#'   in the input dataset, specifying the category of each sample.
#' @param control Name of the control sample category, defined in the
#'   \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control =
#'   'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#'   value in a condition. Can be \code{'mean'} or \code{'median'}.
#' @param test Method to use for p-value calculation: use \code{'wilcoxon'} for
#'   Wilcoxon rank sum test or \code{'shuffle'} for a label shuffling test.
#' @param randomizations Number of random shuffles, used for the label shuffling
#'   test (default = 100).
#' @param ... Further arguments to be passed on for other methods.
#' @return A \code{data.frame} with the mean or median values of splicing
#'   diversity across sample categories and all samples, log2(fold change) of
#'   the two different conditions, raw and Benjamini-Hochberg corrected
#'   p-values.
#' @import methods ggplot2 tidyr
#' @export
#' @details The function calculates diversity changes between two sample
#' conditions, given by the user. It uses the output of the diversity
#' calculation function from this package, which is a \code{data.frame} of
#' splicing diversity values. The \code{data.frame} consists of the diversity
#' values, where the rows are genes and the columns are samples. Additionally,
#' the first column is the gene id. A vector of sample conditions also serves as
#' input, used for aggregating the samples by condition. It calculates the mean
#' or median of the splicing diversity data per sample condition, the difference
#' of these values and the log2 fold change of the two conditions.
#'
#' Furthermore, the user can select a statistical method to calculate the
#' significance of the changes. The p-values and adjusted p-values are
#' calculated by Wilcoxon sum rank test or label shuffling test.
#'
#' The function will exclude genes of low sample size from the significance
#' calculation, depending on which statistical test is applied.
#'
#' @examples
#' # data.frame with splicing diversity values
#' x <- data.frame(Genes = letters[seq_len(10)], matrix(runif(80), ncol = 8))
#'
#' # sample categories
#' samples <- c(rep('Healthy', 4), rep('Pathogenic', 4))
#'
#' # To calculate the difference of splicing diversity changes between the 'Healthy'
#' # and 'Pathogenic' condition together with the significance values, using mean
#' # and Wilcoxon rank sum test, use:
#' calculate_difference(x, samples, control = 'Healthy', method = 'mean', test = 'wilcoxon')
calculate_difference <- function(x, samples, control, method = "mean", test = "wilcoxon", randomizations = 100, 
    ...) {
    if (!is(x, "data.frame")) 
        stop("Input data type is not supported! Please use `?calculate_difference`
         to see the possible arguments and details.", 
            call. = FALSE)
    if (ncol(x) - 1 != length(samples)) 
        stop("The number of columns in the data.frame is not equal to the number of
         samples defined in the samples argument.", 
            call. = FALSE)
    if (length(levels(as.factor(samples))) > 2) 
        stop("The number of conditions are higher than two. Please use exactly two
         different sample conditions, e.g. healthy and pathogenic.", 
            call. = FALSE)
    if (length(levels(as.factor(samples))) < 2) 
        stop("The number of conditions are smaller than two. Please use exactly two
         different sample conditions, e.g. healthy and pathogenic.", 
            call. = FALSE)
    if (!(control %in% samples)) 
        stop("This control sample type cannot be found in your samples.")
    if (!(method %in% c("mean", "median"))) 
        stop("Invalid method. Please use `?calculate_diversity` to see the possible
         arguments and details.", 
            call. = FALSE)
    if (!(test %in% c("wilcoxon", "shuffle"))) 
        stop("Invalid test method. Please use `?calculate_diversity` to see the
         possible arguments and details.", 
            call. = FALSE)
    if (test == "wilcoxon") {
        if (randomizations != 100) 
            message("Note: The 'randomizations' argument is an option for label shuffling,
              it won't have any effect on the Wilcoxon rank sum test.", 
                call. = FALSE)
        if (length(grep(unique(samples)[1], samples)) < 3 | length(grep(unique(samples)[2], samples)) < 3 | length(samples) < 
            8) 
            warning("Low sample size. Wilcoxon rank sum test requires at least
      three samples in a given category and at least 8 samples overall for a
              theoretical p-value smaller than 0.05.", 
                call. = FALSE)
    }
    if (test == "shuffle") {
        if (length(samples) <= 5) 
            warning("Low sample size, not enough samples for label shuffling!", call. = FALSE)
        if (length(samples) > 5 & length(samples) < 10) 
            warning("Low sample size, label shuffling might not give informative and
              correct results.", 
                call. = FALSE)
    }
    x$cond_1 <- apply(x[grep(unique(samples)[1], samples) + 1], 1, function(x) sum(!is.na(x)))
    x$cond_2 <- apply(x[grep(unique(samples)[2], samples) + 1], 1, function(x) sum(!is.na(x)))
    if (test == "wilcoxon") {
        y <- x[x$cond_1 >= 3 & x$cond_2 >= 3 & sum(x$cond_1, x$cond_2) >= 8, ]
        x <- x[x$cond_1 < 3 | x$cond_2 < 3 | sum(x$cond_1, x$cond_2) < 8, ]
    }
    if (test == "shuffle") {
        y <- x[x$cond_1 + x$cond_2 >= 5, ]
        x <- x[x$cond_1 + x$cond_2 < 5, ]
    }
    genes_y <- y[, 1]
    y <- as.matrix(y[, c(-1, -ncol(y), -ncol(y) + 1)])
    genes_x <- x[, 1]
    x <- as.matrix(x[, c(-1, -ncol(x), -ncol(x) + 1)])
    if (nrow(y) != 0) {
        if (nrow(x) != 0) {
            message(paste0("Note: There are ", nrow(x), " genes with low sample size, which will be
    exluded from the statistical testing."))
        }
        if (test == "wilcoxon") {
            p_values <- wilcoxon(y, samples, ...)
        }
        if (test == "shuffle") {
            p_values <- label_shuffling(y, samples, control, method, randomizations)
        }
        y <- data.frame(genes = genes_y, calculate_fc(y, samples, control, method), p_values)
    }
    if (nrow(x) != 0) {
        x <- data.frame(genes = genes_x, calculate_fc(x, samples, control, method), raw_p_values = NA, adjusted_p_values = NA)
        if (nrow(y) != 0) {
            return(rbind(y, x))
        } else {
            return(x)
        }
    } else {
        return(y)
    }
}
