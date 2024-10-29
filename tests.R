library(tidyverse)
MC_runs <- 100000

source("grid_approximation.R")

# Test 1 - normal-normal works

test1_grid <- grid_approximation(
    grid_length = MC_runs,
    grid_min = -100, grid_max = 100,
    prior = c(0, 5),
    prior_type = "dnorm",
    likelihood = c(4, 2),
    likelihood_type = "dnorm",
    seed = 2023
)
summary(test1_grid$posterior)
sd(test1_grid$posterior)

test1_conjugate <- combine_normals(
    0, 5, 4, 2
)
test1_conjugate

# Test 2 - normal-normal works (changing grid length)
test2_grid <- grid_approximation(
    grid_length = MC_runs,
    grid_min = -10, grid_max = 10,
    prior = c(0, 5),
    prior_type = "dnorm",
    likelihood = c(4, 2),
    likelihood_type = "dnorm",
    seed = 2023
)
summary(test2_grid$posterior)

# test 3 - using generated data
test3_grid <- grid_approximation(
    grid_length = MC_runs,
    grid_min = -10, grid_max = 10,
    prior = rnorm(MC_runs, 0, 5),
    prior_type = "data",
    likelihood = rnorm(MC_runs, 4, 2),
    likelihood_type = "data",
    seed = 2023
)
summary(test3_grid$posterior)

# Compare the weightings

test1_grid$weights
test2_grid$weights
test3_grid$weights

calculate_weight <- function(X, Y, Z) {
    if (X == Y) {
        return("X and Y are equal, weight can be any value between 0 and 1.")
    } else {
        w <- (Z - Y) / (X - Y)
        return(w)
    }
}

calculate_weight(0, 5, test1_conjugate$mu_post)

# Test 4 - using skewed data

# Generate two skewed distributions
distribution_1 = rgamma(MC_runs, shape = 2, scale = 2)  # Shape and scale adjusted for positive skew
distribution_2 = rgamma(MC_runs, shape = 2, scale = 3)   # Different scale for the second distribution


test4_grid <- grid_approximation(
    grid_length = MC_runs,
    grid_min = -1, grid_max = 100,
    prior = distribution_1,
    prior_type = "data",
    likelihood = distribution_2,
    likelihood_type = "data",
    seed = 2023
)
mean(test4_grid$posterior)
mean(distribution_1)
mean(distribution_2)

test4_grid$graph

# weighting can't work in this situation
calculate_weight(mean(distribution_1), mean(distribution_2), mean(test4_grid$posterior))

# but this weighting can
test4_grid$weights
