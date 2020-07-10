context("Aggregation of methods")

test_that("Aggregation of diversity calculation is working", {

    read_count_matrix <- rbind(matrix(rpois(36, 6), ncol = 6), matrix(0, nrow = 2, ncol = 6))
    genes <- c("A", "B", "B", "C", "C", "C", "D", "D")

    for (method in c("naive", "laplace", "gini", "simpson", "invsimpson")) {

        for (norm in c(TRUE, FALSE)) {

            diversity <- calculate_method(read_count_matrix, genes, method, norm)

            expect_true(is.data.frame(diversity))
            expect_length(diversity, 7)
            expect_equal(diversity$Gene, c("B", "C", "D"))
            expect_equal(as.character(diversity[3, 2:7]), rep("NA", 6))

            expect_message(calculate_method(read_count_matrix, genes, method, norm, verbose = TRUE), "Note: There are 1 genes with single isoforms,
    which will be exluded from the analysis.")
        }
    }
})
