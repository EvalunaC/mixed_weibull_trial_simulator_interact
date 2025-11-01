# Test script to verify median accuracy improvements
# Tests the simulate_progression and simulate_survival functions

library(survival)

# Source the helper functions from app.R
calculate_scale_from_median <- function(median_time, shape) {
  median_time / ((log(2))^(1/shape))
}

# Test simulate_progression function
simulate_progression <- function(n, nv12, lamb12, target_median, max_followup) {
  tolerance <- 0.1
  max_iterations <- 100

  if (!is.na(lamb12) && lamb12 > 0) {
    scale_param <- (1/lamb12)^(1/nv12)
  } else {
    scale_param <- calculate_scale_from_median(target_median, nv12)
  }

  for (iter in 1:max_iterations) {
    prog_time <- rweibull(n, shape = nv12, scale = scale_param)
    cens_time <- runif(n, 0, max_followup)
    prog_observed <- pmin(prog_time, cens_time)
    prog_event <- as.numeric(prog_time <= cens_time)

    fit <- survfit(Surv(prog_observed, prog_event) ~ 1)
    actual_median <- summary(fit)$table["median"]

    if (is.na(actual_median) || actual_median == 0) {
      actual_median <- median(prog_observed)
    }

    difference <- abs(actual_median - target_median)

    if (difference <= tolerance) {
      cat(sprintf("Progression converged in %d iterations (median: %.2f, target: %.2f)\n",
                  iter, actual_median, target_median))
      break
    }

    adjustment_factor <- target_median / actual_median
    if (adjustment_factor > 1.5) adjustment_factor <- 1.5
    if (adjustment_factor < 0.67) adjustment_factor <- 0.67
    scale_param <- scale_param * adjustment_factor
  }

  list(time = prog_observed, event = prog_event, true_time = prog_time,
       iterations = iter, final_scale = scale_param, final_median = actual_median)
}

# Test simulate_survival function
simulate_survival <- function(n, nv13, lamb13, nv23, lamb23, target_median, prog_times, max_followup) {
  tolerance <- 0.1
  max_iterations <- 100

  cens_time <- runif(n, 0, max_followup)

  if (!is.na(lamb13) && lamb13 > 0) {
    scale_param_13 <- (1/lamb13)^(1/nv13)
  } else {
    scale_param_13 <- calculate_scale_from_median(target_median * 1.5, nv13)
  }

  if (!is.na(lamb23) && lamb23 > 0) {
    scale_param_23 <- (1/lamb23)^(1/nv23)
  } else {
    scale_param_23 <- calculate_scale_from_median(target_median * 0.5, nv23)
  }

  for (iter in 1:max_iterations) {
    death_without_prog <- rweibull(n, shape = nv13, scale = scale_param_13)
    death_after_prog_time <- rweibull(n, shape = nv23, scale = scale_param_23)
    death_after_prog <- prog_times + death_after_prog_time
    os_time <- pmin(death_without_prog, death_after_prog)
    os_observed <- pmin(os_time, cens_time)
    os_event <- as.numeric(os_time <= cens_time)

    fit <- survfit(Surv(os_observed, os_event) ~ 1)
    actual_median <- summary(fit)$table["median"]

    if (is.na(actual_median) || actual_median == 0) {
      actual_median <- median(os_observed)
    }

    difference <- abs(actual_median - target_median)

    if (difference <= tolerance) {
      cat(sprintf("Survival converged in %d iterations (median: %.2f, target: %.2f)\n",
                  iter, actual_median, target_median))
      break
    }

    adjustment_factor <- target_median / actual_median
    if (adjustment_factor > 1.5) adjustment_factor <- 1.5
    if (adjustment_factor < 0.67) adjustment_factor <- 0.67
    scale_param_13 <- scale_param_13 * adjustment_factor
    scale_param_23 <- scale_param_23 * adjustment_factor

    if (scale_param_13 <= 0 || scale_param_13 > 10000 ||
        scale_param_23 <= 0 || scale_param_23 > 10000) {
      cat(sprintf("Warning: Scale parameters out of bounds at iteration %d\n", iter))
      break
    }
  }

  list(time = os_observed, event = os_event, true_time = os_time,
       iterations = iter, final_median = actual_median)
}

# Run tests
cat("=== Testing Progression Simulation ===\n")
set.seed(123)
prog_result <- simulate_progression(
  n = 200,
  nv12 = 1.5,
  lamb12 = NA,
  target_median = 6,
  max_followup = 24
)
cat(sprintf("Final progression median: %.2f (target: 6.00)\n\n", prog_result$final_median))

cat("=== Testing Survival Simulation ===\n")
set.seed(123)
# Use the progression times from previous simulation
surv_result <- simulate_survival(
  n = 200,
  nv13 = 1.5,
  lamb13 = NA,
  nv23 = 1.2,
  lamb23 = NA,
  target_median = 12,
  prog_times = prog_result$true_time,
  max_followup = 24
)
cat(sprintf("Final survival median: %.2f (target: 12.00)\n\n", surv_result$final_median))

cat("=== Test Complete ===\n")
cat(sprintf("Progression: %.2f months (target: 6.00, error: %.2f)\n",
            prog_result$final_median, abs(prog_result$final_median - 6)))
cat(sprintf("Survival: %.2f months (target: 12.00, error: %.2f)\n",
            surv_result$final_median, abs(surv_result$final_median - 12)))
