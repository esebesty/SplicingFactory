context("Basic diversity calculations")

test_that("Basic diversity calculation is working", {

    read_counts <- c(0, 0, 5, 4, 1)

    naive_entropy <- calculate_entropy(read_counts)

    expect_true(is.numeric(naive_entropy))
    expect_length(naive_entropy, 1)

    naive_entropy_norm <- calculate_entropy(read_counts, norm = FALSE)

    expect_true(is.numeric(naive_entropy_norm))
    expect_length(naive_entropy_norm, 1)
    expect_equal(naive_entropy_norm, 1.360964, tolerance = 0.001, scale = 1)

    laplace_entropy <- calculate_entropy(read_counts, pseudocount = 1)

    expect_true(is.numeric(laplace_entropy))
    expect_length(laplace_entropy, 1)
    expect_equal(laplace_entropy, 0.8465362, tolerance = 0.001, scale = 1)

    laplace_entropy_norm <- calculate_entropy(read_counts, norm = FALSE,
                                              pseudocount = 1)

    expect_true(is.numeric(laplace_entropy_norm))
    expect_length(laplace_entropy_norm, 1)
    expect_equal(laplace_entropy_norm, 1.965596, tolerance = 0.001, scale = 1)

    gini_coef <- calculate_gini(read_counts)

    expect_true(is.numeric(gini_coef))
    expect_length(gini_coef, 1)
    expect_equal(gini_coef, 0.7, tolerance = 0.001, scale = 1)

    simpson_index <- calculate_simpson(read_counts)

    expect_true(is.numeric(simpson_index))
    expect_length(simpson_index, 1)
    expect_equal(simpson_index, 0.58, tolerance = 0.001, scale = 1)

    invsimpson_index <- calculate_inverse_simpson(read_counts)

    expect_true(is.numeric(invsimpson_index))
    expect_length(invsimpson_index, 1)
    expect_equal(invsimpson_index, 2.380952, tolerance = 0.001, scale = 1)

})

test_that("Diversity calculation of single isoforms", {

    read_counts_single <- c(1)

    naive_entropy_single <- calculate_entropy(read_counts_single)

    expect_true(is.nan(naive_entropy_single))

})

test_that("Diversity calculation of zero expression", {

    read_counts_zero <- c(0, 0)

    naive_entropy_zero <- calculate_entropy(read_counts_zero)

    expect_true(is.na(naive_entropy_zero))

})

test_that("Tsallis entropy calculation is working", {
    read_counts <- c(0, 0, 5, 4, 1)
    tsallis_entropy <- calculate_tsallis_entropy(read_counts, q = 2)
    expect_true(is.numeric(tsallis_entropy))
    expect_length(tsallis_entropy, 1)
    # For q=2, expected value can be calculated manually
    # p = c(0,0,0.5,0.4,0.1), sum(p^2) = 0.25+0.16+0.01=0.42, (1-0.42)/(2-1)=0.58
    expect_equal(tsallis_entropy, 0.58, tolerance = 0.001, scale = 1)
    tsallis_entropy_norm <- calculate_tsallis_entropy(read_counts, q = 2, norm = FALSE)
    expect_true(is.numeric(tsallis_entropy_norm))
    expect_length(tsallis_entropy_norm, 1)
    expect_equal(tsallis_entropy_norm, 0.58, tolerance = 0.001, scale = 1)
    # q=1 should error
    expect_error(calculate_tsallis_entropy(read_counts, q = 1))
    # Single isoform
    expect_true(is.nan(calculate_tsallis_entropy(c(1), q = 2)))
    # Zero expression
    expect_true(is.na(calculate_tsallis_entropy(c(0,0), q = 2)))
})
