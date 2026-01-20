context("Main diversity calculation")

test_that("Basic input error handling is working.", {

    for (method in c("naive", "laplace", "gini", "simpson", "invsimpson")) {

        count_matrix <- data.frame(A = letters[1:8], B = letters[1:8])
        genes <- c("A", "B", "B", "C", "C", "C", "D", "D")
        norm <- TRUE
        tpm <- FALSE

        expect_error(calculate_diversity(count_matrix, genes, method, norm, tpm), "Input data  must be numeric!")

        count_matrix <- data.frame(A = c(1, 2, 3, rep(NA, 5)))

        expect_error(calculate_diversity(count_matrix, genes, method, norm, tpm), "The data contains NA as expression values. NAs are not allowed in the
         input.")

        count_matrix <- matrix(rpois(60, 10), ncol = 6)
        colnames(count_matrix) <- paste0("Sample", 1:6)

        expect_error(calculate_diversity(count_matrix, genes, method, norm, tpm), "The number of rows is not equal to the given gene set.")

        genes <- c(rep("Gene1", 3), rep("Gene2", 2), rep("Gene3", 3), rep("Gene4", 2))

        expect_message(calculate_diversity(count_matrix, genes, method, norm, tpm = TRUE, verbose = TRUE), "Note: tpm as a logical argument is only interpreted in case of
            tximport lists.")

        expect_error(calculate_diversity(count_matrix, genes, "calculation", norm, tpm), "Invalid method. Please use `?calculate_diversity` to see the possible
         arguments and details.",
            fixed = TRUE)

        for (method in c("gini", "simpson", "invsimpson")) {

            expect_message(calculate_diversity(count_matrix, genes, method = "gini", norm = FALSE, tpm, verbose = TRUE))
        }
    }
})

test_that("calculate_diversity supports q as a vector and returns correct metadata", {
    x <- matrix(c(0, 0, 5, 4, 1, 2, 2, 2, 2, 2), ncol = 2)
    colnames(x) <- c("Sample1", "Sample2")
    gene <- c("Gene1", "Gene1", "Gene1", "Gene1", "Gene1")
    qvec <- c(1.1, 1.5, 2)
    result <- calculate_diversity(x, gene, method = "tsallis", norm = TRUE, q = qvec)
    # Assay should have columns for each q and sample
    assay_names <- colnames(assay(result))
    for (qi in qvec) {
        expect_true(any(grepl(paste0("q=", qi), assay_names)))
    }
    # Metadata should contain the q vector
    expect_equal(metadata(result)$q, qvec)
})
test_that("calculate_diversity passes q parameter for Tsallis entropy", {
    x <- matrix(c(0, 0, 5, 4, 1, 2, 2, 2, 2, 2), ncol = 2)
    colnames(x) <- c("Sample1", "Sample2")
    gene <- c("Gene1", "Gene1", "Gene1", "Gene1", "Gene1")
    # Calculate with q = 2
    result_q2 <- calculate_diversity(x, gene, method = "tsallis", norm = TRUE, q = 2)
    # Calculate with q = 1.5
    result_q15 <- calculate_diversity(x, gene, method = "tsallis", norm = TRUE, q = 1.5)
    # Extract values
    val_q2 <- assay(result_q2)[1, 1]
    val_q15 <- assay(result_q15)[1, 1]
    expect_true(is.numeric(val_q2))
    expect_true(is.numeric(val_q15))
    expect_false(is.na(val_q2))
    expect_false(is.na(val_q15))
    expect_false(abs(val_q2 - val_q15) < 1e-8) # Should be different for different q
})
