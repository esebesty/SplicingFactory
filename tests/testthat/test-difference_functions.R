context("Basic difference calculations")

diversity_1 <- matrix(runif(80), ncol = 8)
diversity_2 <- data.frame(S1 = 0.1, S2 = 0.2, S3 = 0.3, S4 = 0.4, S5 = 0.5, S6 = 0.6, S7 = 0.7, S8 = 0.8)
samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
control <- "Healthy"

test_that("Fold change calculation is correct", {
    
    for (method in c("mean", "median")) {
        
        fold_change <- calculate_fc(diversity_1, samples, control, "mean")
        
        expect_length(fold_change, 4)
        expect_true(is.data.frame(fold_change))
        
        fold_change <- calculate_fc(as.matrix(diversity_2), samples, control, "mean")
        
        expect_equal(fold_change$Pathogenic_mean, 0.65, tolerance = 0.001, scale = 1)
        expect_equal(fold_change$Healthy_mean, 0.25, tolerance = 0.001, scale = 1)
        expect_equal(fold_change$mean_difference, 0.4, tolerance = 0.001, scale = 1)
        expect_equal(fold_change$log2_fold_change, 1.378512, tolerance = 0.001, scale = 1)
    }
})

test_that("Wilcoxon sum rank test is correct", {
    
    wilcoxon_result <- wilcoxon(diversity_1, samples)
    
    expect_length(wilcoxon_result, 20)
    expect_true(is.matrix(wilcoxon_result))
    
    wilcoxon_result <- wilcoxon(as.matrix(diversity_2), samples)
    
    expect_equal(as.numeric(wilcoxon_result[1, 1]), 0.03038282, tolerance = 0.001, scale = 1)
    expect_equal(as.numeric(wilcoxon_result[1, 2]), 0.03038282, tolerance = 0.001, scale = 1)
    
})

test_that("Label shuffling test is correct", {
    
    shuffling_result <- label_shuffling(diversity_1, samples, control, "mean")
    
    expect_length(shuffling_result, 20)
    expect_true(is.matrix(shuffling_result))
    
    diversity_2 <- rbind(diversity_2, data.frame(S1 = 0.2, S2 = 0.3, S3 = 0.4, S4 = 0.5, S5 = 0.6, S6 = 0.7, 
        S7 = 0.8, S8 = 0.9))
    
    shuffling_result <- label_shuffling(as.matrix(diversity_2), samples, control, "mean")
    
    expect_equal(as.numeric(shuffling_result[1, 1]), 1, tolerance = 0.001, scale = 1)
    expect_equal(as.numeric(shuffling_result[1, 2]), 1, tolerance = 0.001, scale = 1)
})
