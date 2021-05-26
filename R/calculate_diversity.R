#' Main function for calculating splicing diversity
#'
#' @param x A numeric \code{matrix}, \code{data.frame}, \code{tximport} list,
#'   \code{DGEList}, \code{SummarizedExperiment} or \code{ExpressionSet}.
#' @param genes Character vector with equal length to the number of rows of the
#'   input dataset with transcript-level expression values. The values in
#'   \code{x} are grouped into genes based on this vector.
#' @param method Method to use for splicing diversity calculation, including
#'   naive entropy (\code{naive}), Laplace entropy (\code{laplace}), Gini index
#'   (\code{gini}), Simpson index (\code{simpson}) and inverse Simpson index
#'   (\code{invsimpson}). The default method is Laplace entropy.
#' @param norm If \code{TRUE}, the entropy values are normalized to the number
#'   of transcripts for each gene. The normalized entropy values are always
#'   between 0 and 1. If \code{FALSE}, genes cannot be compared to each other,
#'   due to possibly different maximum entropy values.
#' @param tpm In the case of a tximport list, TPM values or raw read counts can
#'   serve as an input. If \code{TRUE}, TPM values will be used, if
#'   \code{FALSE}, read counts will be used.
#' @param assayno An integer value. In case of multiple assays in a
#'    \code{SummarizedExperiment} input, the argument specifies the assay number
#'    to use for diversity calculations.
#' @param verbose If \code{TRUE}, the function will print additional diagnostic
#'    messages, besides the warnings and errors.
#' @return Gene-level splicing diversity values in a \code{SummarizedExperiment}
#'    object.
#' @import methods
#' @importFrom SummarizedExperiment SummarizedExperiment assays assay
#' @export
#' @details The function is intended to process transcript-level expression data
#' from RNA-seq or similar datasets.
#'
#' Given a N x M matrix or similar data structure, where the N rows are
#' transcripts and the M columns are samples, and a vector of gene ids, used for
#' aggregating the transcript level data, the function calculates transcript
#' diversity values for each gene in each sample. These diversity values can be
#' used to investigate the dominance of a specific transcript for a gene,
#' the diversity of transcripts in a gene, and analyze changes in diversity.
#'
#' There are a number of diversity values implemented in the package. These
#' include the following:
#' \itemize{
#'   \item Naive entropy: Shannon entropy using the transcript frequencies as
#'     probabilities. 0 entropy means a single dominant transcript, higher
#'     values mean a more diverse set of transcripts for a gene.
#'   \item Laplace entropy: Shannon entropy where the transcript frequencies are
#'     replaced by a Bayesian estimate, using Laplace's prior.
#'   \item Gini index: a measure of statistical dispersion originally used in
#'     economy. This measurement ranges from 0 (complete equality) to 1
#'     (complete inequality). A value of 1 (complete inequality) means a single
#'     dominant transcript.
#'   \item Simpson index: a measure of diversity, characterizing the number of
#'     different species (transcripts of a gene) in a dataset. Originally, this
#'     measurement calculates the probability that randomly selected individuals
#'     belong to different species. Simpson index ranges between 0 and 1; the
#'     higher the value, the higher the diversity.
#'   \item Inverse Simpson index: Similar concept as the Simpson index,
#'     although a higher inverse-Simpson index means greater diversity. It
#'     ranges between 1 and the total number of transcripts for a gene.
#' }
#'
#' The function can calculate the gene level diversity index using any kind of
#' expression measure, including raw read counts, FPKM, RPKM or TPM values,
#' although results may vary.
#'
#' @examples
#' # matrix with RNA-seq read counts
#' x <- matrix(rpois(60, 10), ncol = 6)
#' colnames(x) <- paste0("Sample", 1:6)
#'
#' # gene names used for grouping the transcript level data
#' gene <- c(rep("Gene1", 3), rep("Gene2", 2), rep("Gene3", 3), rep("Gene4", 2))
#'
#' # calculating normalized Laplace entropy
#' result <- calculate_diversity(x, gene, method = "laplace", norm = TRUE)
calculate_diversity <- function(x, genes = NULL, method = "laplace", norm = TRUE,
                                tpm = FALSE, assayno = 1, verbose = FALSE) {
  if (!(is.matrix(x) || is.data.frame(x) || is.list(x) || is(x, "DGEList") ||
    is(x, "RangedSummarizedExperiment") || is(x, "SummarizedExperiment"))) {
    stop("Input data type is not supported! Please use `?calculate_diversity`
\t to see the possible arguments and details.")
  }

  if (is(x, "data.frame")) {
    x <- as.matrix(x)
  }

  if (tpm == TRUE && !is.list(x) && verbose == TRUE) {
    message("Note: tpm as a logical argument is only interpreted in case of
            tximport lists.")
  }

  if (is.list(x)) {
    if (length(x) == 4 && "counts" %in% names(x)) {
      if (tpm == FALSE) {
        x <- as.matrix(x$counts)
      }
      if (tpm == TRUE) {
        x <- as.matrix(x$abundance)
      }
    } else if (is(x, "DGEList")) {
      x <- as.matrix(x$counts)
      if (verbose == TRUE) {
        message("Note: calculate_diversity methods are only applicable if your
                DGEList contains transcript-level expression data.")
      }
      if (tpm == TRUE && verbose == TRUE) {
        message("Note: tpm as a logical argument is only interpreted in case of
              tximport lists.")
      }
    } else {
      stop("The package cannot find any expression data in your input.", call. = FALSE)
    }
  }

  if (is(x, "RangedSummarizedExperiment") || is(x, "SummarizedExperiment")) {
    if (!is.numeric(assayno) | length(SummarizedExperiment::assays(x)) < assayno) {
      stop("Please give a valid number to pick an assay from your data.", call. = FALSE)
    }
    else if (is.numeric(assayno)) {
      x <- as.matrix(SummarizedExperiment::assays(x)[[assayno]])
    }
    if (is.null(genes)) {
      genes <- rownames(x)
      rownames(x) <- NULL
      if (is.null(genes)) {
        stop("Please construct a valid gene set for your SummarizedExperiment.",
          call. = FALSE
        )
      }
    }
  }

  if (!is.numeric(x)) {
    stop("Input data  must be numeric!", call. = FALSE)
  }

  if (any(is.na(x))) {
    stop("The data contains NA as expression values. NAs are not allowed in the
         input.", call. = FALSE)
  }

  if (nrow(x) != length(genes)) {
    stop("The number of rows is not equal to the given gene set.", call. = FALSE)
  }

  if (!(method %in% c("naive", "laplace", "gini", "simpson", "invsimpson"))) {
    stop("Invalid method. Please use `?calculate_diversity` to see the possible
         arguments and details.",
      call. = FALSE
    )
  }

  if (method == "gini" && norm == FALSE && verbose == TRUE) {
    message("Gini coefficient ranges between 0 (complete equality) and 1 (complete
    inequality). The 'norm' logical argument does not have any effect on the
            calculation.", call. = FALSE)
  }

  if (method == "simpson" && norm == FALSE && verbose == TRUE) {
    message("Simpson index ranges between 0 to 1. The 'norm' logical argument
            does not have any effect on the calculation.",
      call. = FALSE
    )
  }

  if (method == "invsimpson" && norm == FALSE && verbose == TRUE) {
    message("Inverse Simpson index does not use the 'norm' argument, and it won't
            have any effect on the calculation.", call. = FALSE)
  }

  result <- calculate_method(x, genes, method, norm, verbose = verbose)

  result_assay <- result[, -1, drop = FALSE]
  result_rowData <- data.frame(genes = result[, 1], row.names = result[, 1])
  result_colData <- data.frame(samples = colnames(x), row.names = colnames(x))
  result_metadata <- list(method = method, norm = norm)

  result <- SummarizedExperiment(assays = list(diversity = result_assay),
                                 rowData = result_rowData,
                                 colData = result_colData,
                                 metadata = result_metadata)

  return(result)
}
