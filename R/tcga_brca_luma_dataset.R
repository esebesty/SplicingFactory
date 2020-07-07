#' TCGA Luminal A breast cancer dataset
#'
#' Data from The Cancer Genome Atlas, downloaded on 12th April, 2020. It
#' contains transcript level read counts of 20 patient with Luminal A type
#' breast cancer (primary tumor and solid normal samples).
#'
#' @docType data
#'
#' @usage data(example_dataset)
#'
#' @format A data frame with 1054 rows and 41 columns. The first 3 columns
#' contain the following:
#' \describe{
#'   \item{genes}{Gene names}
#'   \item{TCGA-A7-A0CH_N}{TCGA-A7-A0CH patient normal tissue sample}
#'   \item{TCGA-A7-A0CH_T}{TCGA-A7-A0CH patient tumor tissue sample}
#' }.
#'
#' @keywords datasets
#'
#' @references The Cancer Genome Atlas Network (2012) Nature 490, 61–70
#'   (\href{https://doi.org/10.1038/nature11412}{doi:10.1038/nature11412})
#'
#' @source \href{https://portal.gdc.cancer.gov/legacy-archive}{TCGA Legacy}
#'
"tcga_brca_luma_dataset"
