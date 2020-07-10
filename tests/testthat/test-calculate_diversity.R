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
