library(testthat)
set.seed(20191015)

sim1_weights <- run_sim(patterns = diag(10), frequency = 1:10,
                        duration = 10:1, lrate_onset = 0.05,
                        lrate_drop_time = 2,
                        lrate_drop_perc = 0)$weight_matrix
sim1_weights_before <- as.matrix(read.csv("sim1_weights.csv"))

test_that("sim 1 produces same results as before", {
            expect_equal(sim1_weights, sim1_weights_before,
                         check.attributes = F)
          })

duration <- c(4, 2, 1, 8, 4, 2, 12, 6, 3)
frequency <- c(2, 4, 8, 2, 4, 8, 2, 4, 8)
lrate_drop_perc <- seq(0, 1, 0.04)
sim4 <- lapply(lrate_drop_perc, function(x)
  run_exp(frequency, duration, 0.05, 2, x, diag(9), 10, 0.1))

sim4 <- plyr::ldply(sim4, "data.frame")
sim4 <- cbind(sim4, lrate_drop_perc)

sim4_before <- read.csv("sim4.csv")
test_that("sim 4 produces same results as before", {
            expect_equal(sim4, sim4_before, check.attributes = F)
          })
