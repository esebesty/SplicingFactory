context("Difference calculation")

test_that("Difference calculation methods are correct", {

    diversity <- data.frame(Genes = letters[1:10], matrix(runif(80), ncol = 8))
    samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
    control <- "Healthy"

    expect_error(calculate_difference(diversity, samples, control, "Unknown method", test), "Invalid method. Please use `?calculate_diversity` to see the possible
         arguments and details.",
        fixed = TRUE)

    expect_error(calculate_difference(diversity, samples, control, "mean", "bootstrap"), "Invalid test method. Please use `?calculate_diversity` to see the
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
