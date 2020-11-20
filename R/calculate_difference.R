#' Calculate splicing diversity changes between two conditions.
#'
#' @param x A \code{SummarizedExperiment} with splicing diversity values for
#'   each gene in each sample or a \code{data.frame} with gene names in the
#'   first column and splicing diversity values for each sample in additional
#'   columns.
#' @param samples A vector of length one, specifying the column name of the
#'   \code{colData} annotation column from the \code{SummarizedExperiment}
#'   object, that should be used as the category column or a character vector
#'   with an equal length to the number of columns in the input dataset,
#'   specifying the category of each sample in the case of a \code{data.frame}
#'   input.
#' @param control Name of the control sample category, defined in the
#'   \code{samples} vector, e.g. \code{control = 'Normal'} or \code{control =
#'   'WT'}.
#' @param method Method to use for calculating the average splicing diversity
#'   value in a condition. Can be \code{'mean'} or \code{'median'}.
#' @param test Method to use for p-value calculation: use \code{'wilcoxon'} for
#'   Wilcoxon rank sum test or \code{'shuffle'} for a label shuffling test.
#' @param randomizations Number of random shuffles, used for the label shuffling
#'   test (default = 100).
#' @param pcorr P-value correction method applied to the Wilcoxon rank sum test
#'   or label shuffling test results, as defined in the \code{p.adjust}
#'   function.
#' @param assayno An integer value. In case of multiple assays in a
#'    \code{SummarizedExperiment} input, the argument specifies the assay number
#'    to use for difference calculations.
#' @param verbose If \code{TRUE}, the function will print additional diagnostic
#'    messages.
#' @param ... Further arguments to be passed on for other methods.
#' @return A \code{data.frame} with the mean or median values of splicing
#'   diversity across sample categories and all samples, log2(fold change) of
#'   the two different conditions, raw and corrected p-values.
#' @import methods
#' @importFrom SummarizedExperiment SummarizedExperiment assays assay colData
#' @export
#' @details The function calculates diversity changes between two sample
#' conditions. It uses the output of the diversity calculation function, which
#' is a \code{SummarizedExperiment} object of splicing diversity values.
#' Additionally, it can use a \code{data.frame} as input, where the first column
#' contains gene names, and all additional columns contain splicing diversity
#' values for each sample. A vector of sample conditions also serves as input,
#' used for aggregating the samples by condition.
#'
#' It calculates the mean or median of the splicing diversity data per sample
#' condition, the difference of these values and the log2 fold change of the two
#' conditions. Furthermore, the user can select a statistical method to
#' calculate the significance of the changes. The p-values and adjusted p-values
#' are calculated using a Wilcoxon sum rank test or label shuffling test.
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
#' # To calculate the difference of splicing diversity changes between the
#' # 'Healthy' and 'Pathogenic' condition together with the significance values,
#' # using mean and Wilcoxon rank sum test, use:
#' calculate_difference(x, samples, control = 'Healthy', method = 'mean', test = 'wilcoxon')
calculate_difference <- function(x, samples, control, method = "mean",
                                 test = "wilcoxon", randomizations = 100,
                                 pcorr = "BH", assayno = 1, verbose = FALSE,
                                 ...) {
    if (!(is(x, "data.frame") || is(x, "RangedSummarizedExperiment") || is(x, "SummarizedExperiment"))) {
        stop("Input data type is not supported! Please use `?calculate_difference`
         to see the possible arguments and details.",
            call. = FALSE)
    }

    if (is(x, "RangedSummarizedExperiment") || is(x, "SummarizedExperiment")) {
      if (length(samples) != 1) {
        stop("In the case of SummarizedExperiment input, the samples argument
         must be a single character value, specifying the colData column name
         from the SummarizeExperiment object, that should be used as sample
         categories", call. = FALSE)
      }
      samples <- colData(x)[[samples]]

      if (!is.numeric(assayno) | length(SummarizedExperiment::assays(x)) < assayno) {
        stop("Please give a valid number to pick an assay from your data.", call. = FALSE)
      }
      else if (is.numeric(assayno)) {
        x <- as.data.frame(SummarizedExperiment::assays(x)[[assayno]])
      }
      genes <- rownames(x)
      rownames(x) <- NULL
      x <- cbind(genes, x)
    }

    if (ncol(x) - 1 != length(samples)) {
        stop("The number of columns in the data.frame is not equal to the number of
          samples defined in the samples argument.",
            call. = FALSE)
    }

    if (length(levels(as.factor(samples))) > 2) {
        stop("The number of conditions are higher than two. Please use exactly two
         different sample conditions, e.g. healthy and pathogenic.",
            call. = FALSE)
    }

    if (length(levels(as.factor(samples))) < 2) {
        stop("The number of conditions are smaller than two. Please use exactly two
         different sample conditions, e.g. healthy and pathogenic.",
            call. = FALSE)
    }

    if (!(control %in% samples)) {
        stop("This control sample type cannot be found in your samples.")
    }

    if (!(method %in% c("mean", "median"))) {
        stop("Invalid method. Please use `?calculate_difference` to see the possible
         arguments and details.",
            call. = FALSE)
    }

    if (!(test %in% c("wilcoxon", "shuffle"))) {
        stop("Invalid test method. Please use `?calculate_difference` to see the
         possible arguments and details.",
            call. = FALSE)
    }

    if (test == "wilcoxon") {
        if (randomizations != 100 && verbose == TRUE)
            message("Note: The 'randomizations' argument is an option for label shuffling,
              it won't have any effect on the Wilcoxon rank sum test.",
                call. = FALSE)
        if (length(grep(unique(samples)[1], samples)) < 3 |
            length(grep(unique(samples)[2], samples)) < 3 |
            length(samples) < 8)
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

    if (!(pcorr %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                       "fdr", "none"))) {
      stop("Invalid p-value correction method. Please use `?calculate_difference` to see the
         possible arguments and details.",
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
        if (nrow(x) != 0 && verbose == TRUE) {
            message(paste0("Note: There are ", nrow(x), " genes with low sample size, which will be
    exluded from the statistical testing."))
        }
        if (test == "wilcoxon") {
            p_values <- wilcoxon(y, samples, ...)
        }
        if (test == "shuffle") {
            p_values <- label_shuffling(y, samples, control, method,
                                        randomizations)
        }
        y <- data.frame(genes = genes_y, calculate_fc(y, samples, control,
                                                      method), p_values)
    }

    if (nrow(x) != 0) {
        x <- data.frame(genes = genes_x, calculate_fc(x, samples, control,
                                                      method), raw_p_values = NA,
                        adjusted_p_values = NA)
        if (nrow(y) != 0) {
            return(rbind(y, x))
        } else {
            return(x)
        }
    } else {
        return(y)
    }
}
