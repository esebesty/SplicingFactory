context("Difference calculation")

test_that("Difference calculation methods are correct", {

    diversity <- data.frame(Genes = letters[1:10], matrix(runif(80), ncol = 8))
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    control <- "Healthy"
    test <- "wilcoxon"
    pcorr <- "BH"

    expect_error(calculate_difference(diversity, samples, control, "Unknown method", test), "Invalid method. Please use `?calculate_difference` to see the possible
         arguments and details.",
        fixed = TRUE)

    expect_error(calculate_difference(diversity, samples, control, "mean", "bootstrap"), "Invalid test method. Please use `?calculate_difference` to see the
         possible arguments and details.",
        fixed = TRUE)

    expect_error(calculate_difference(diversity, samples, control, "mean", test, 100, "Unknown correction"), "Invalid p-value correction method. Please use `?calculate_difference` to see the
         possible arguments and details.",
        fixed = TRUE)
})

test_that("Difference calculation input handling is working.", {

    for (method in c("mean", "median")) {

        for (test in c("wilcoxon", "shuffle")) {

            diversity <- matrix(rpois(60, 10), ncol = 6)
            samples <- c(rep("Healthy", 4), rep("Pathogenic", 5))
            control <- "Healthy"

            expect_error(calculate_difference(diversity, samples, control, method, test), "Input data type is not supported! Please use `?calculate_difference`
         to see the possible arguments and details.",
                fixed = TRUE)

            diversity <- data.frame(Genes = letters[1:10], matrix(runif(80), ncol = 8))


            expect_error(calculate_difference(diversity, samples, control, method, test), "The number of columns in the data.frame is not equal to the number of
          samples defined in the samples argument.",
                fixed = TRUE)

            samples <- c(rep("Healthy", 4), rep("Pathogenic", 2), rep("Completely other type of biological condition.",
                2))

            expect_error(calculate_difference(diversity, samples, control, method, test), "The number of conditions are higher than two. Please use exactly two
         different sample conditions, e.g. healthy and pathogenic.",
                fixed = TRUE)

            samples <- c(rep("Healthy", 8))

            expect_error(calculate_difference(diversity, samples, control, method, test), "The number of conditions are smaller than two. Please use exactly two
         different sample conditions, e.g. healthy and pathogenic.",
                fixed = TRUE)

            samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))

            expect_error(calculate_difference(diversity, samples, "Healthy control", method, test), "This control sample type cannot be found in your samples.",
                fixed = TRUE)
        }
    }
})

test_that("Sample size warnings are working.", {

    for (method in c("mean", "median")) {

        diversity <- data.frame(Genes = letters[1:10], matrix(runif(80), ncol = 8))
        samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
        control <- "Healthy"
        test <- "wilcoxon"

        expect_message(calculate_difference(diversity, samples, control, method, test, 1000, verbose = TRUE), "Note: The 'randomizations' argument is an option for label shuffling,
              it won't have any effect on the Wilcoxon rank sum test.")

        diversity <- data.frame(Genes = letters[1:10], matrix(runif(40), ncol = 4))
        samples <- c(rep("Healthy", 2), rep("Pathogenic", 2))

        expect_warning(calculate_difference(diversity, samples, control, method, test), "Low sample size. Wilcoxon rank sum test requires at least
      three samples in a given category and at least 8 samples overall for a
              theoretical p-value smaller than 0.05.",
            fixed = TRUE)

        test <- "shuffle"

        expect_warning(calculate_difference(diversity, samples, control, method, test), "Low sample size, not enough samples for label shuffling!")
    }
})

test_that("Calculate difference output is correct.", {

    diversity <- data.frame(Genes = letters[1], S1 = 0.1, S2 = 0.2, S3 = 0.3, S4 = 0.4, S5 = 0.5, S6 = 0.6, S7 = 0.7,
        S8 = 0.8)
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    control = "Healthy"

    result <- calculate_difference(diversity, samples, control)

    expect_true(is.data.frame(result))
    expect_length(result, 7)
    expect_equal(mean(result$Pathogenic_mean), 0.65, tolerance = 0.001, scale = 1)
    expect_equal(mean(result$Healthy_mean), 0.25, tolerance = 0.001, scale = 1)
})

test_that("Tsallis entropy works with calculate_difference and Wilcoxon test.", {

    # Create simple diversity data from tsallis entropy
    diversity <- data.frame(
        Genes = letters[1:5],
        S1 = c(0.3, 0.4, 0.2, 0.5, 0.35),
        S2 = c(0.32, 0.42, 0.22, 0.48, 0.33),
        S3 = c(0.31, 0.41, 0.21, 0.49, 0.34),
        S4 = c(0.33, 0.43, 0.23, 0.47, 0.36),
        S5 = c(0.6, 0.7, 0.5, 0.8, 0.65),
        S6 = c(0.62, 0.72, 0.52, 0.78, 0.67),
        S7 = c(0.61, 0.71, 0.51, 0.79, 0.66),
        S8 = c(0.63, 0.73, 0.53, 0.77, 0.68)
    )
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    control <- "Healthy"

    result <- calculate_difference(diversity, samples, control, method = "mean", test = "wilcoxon")

    expect_true(is.data.frame(result))
    expect_true(all(c("genes", "Healthy_mean", "Pathogenic_mean", "log2_fold_change", "raw_p_values", "adjusted_p_values") %in% colnames(result)))
    expect_equal(nrow(result), 5)
    # Pathogenic should have higher diversity values on average
    expect_gt(mean(result$Pathogenic_mean), mean(result$Healthy_mean))
    # P-values should be numeric and non-negative
    expect_true(all(result$raw_p_values >= 0 & result$raw_p_values <= 1, na.rm = TRUE))
})

test_that("Tsallis entropy works with calculate_difference and median method.", {

    diversity <- data.frame(
        Genes = letters[1:10],
        S1 = runif(10, 0.2, 0.4),
        S2 = runif(10, 0.2, 0.4),
        S3 = runif(10, 0.2, 0.4),
        S4 = runif(10, 0.2, 0.4),
        S5 = runif(10, 0.6, 0.8),
        S6 = runif(10, 0.6, 0.8),
        S7 = runif(10, 0.6, 0.8),
        S8 = runif(10, 0.6, 0.8)
    )
    samples <- c(rep("Control", 4), rep("Treatment", 4))
    control <- "Control"

    result <- calculate_difference(diversity, samples, control, method = "median", test = "wilcoxon")

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 10)
    expect_true(all(c("genes", "Control_median", "Treatment_median") %in% colnames(result)))
    expect_true(all(!is.na(result$raw_p_values)))
})

test_that("Tsallis entropy works with calculate_difference and label shuffling.", {

    # Use larger sample size (5+5=10) to avoid low sample size warning (requires at least 8 total)
    diversity <- data.frame(
        Genes = letters[1:8],
        S1 = c(0.25, 0.35, 0.45, 0.55, 0.30, 0.40, 0.50, 0.60),
        S2 = c(0.28, 0.38, 0.48, 0.58, 0.32, 0.42, 0.52, 0.62),
        S3 = c(0.27, 0.37, 0.47, 0.57, 0.31, 0.41, 0.51, 0.61),
        S4 = c(0.29, 0.39, 0.49, 0.59, 0.33, 0.43, 0.53, 0.63),
        S5 = c(0.70, 0.80, 0.90, 1.0, 0.75, 0.85, 0.95, 1.0),
        S6 = c(0.72, 0.82, 0.92, 0.98, 0.77, 0.87, 0.97, 0.99),
        S7 = c(0.75, 0.85, 0.95, 0.96, 0.80, 0.90, 1.0, 1.0),
        S8 = c(0.68, 0.78, 0.88, 0.99, 0.73, 0.83, 0.93, 1.0),
        S9 = c(0.26, 0.36, 0.46, 0.56, 0.31, 0.41, 0.51, 0.61),
        S10 = c(0.71, 0.81, 0.91, 0.99, 0.76, 0.86, 0.96, 1.0)
    )
    samples <- c(rep("WT", 5), rep("Mutant", 5))
    control <- "WT"

    result <- calculate_difference(
        diversity, samples, control,
        method = "mean", test = "shuffle", randomizations = 50
    )

    expect_true(is.data.frame(result))
    expect_true(all(c("genes", "WT_mean", "Mutant_mean", "log2_fold_change", "raw_p_values", "adjusted_p_values") %in% colnames(result)))
    expect_equal(nrow(result), 8)
})
