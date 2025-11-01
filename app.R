#
# Shiny App for Mixed Weibull Trial Simulator
# Simulates progression-free survival (PFS) and overall survival (OS) data
# with Weibull mixture distributions
#

library(shiny)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(DT)

# UI definition
ui <- fluidPage(
  titlePanel("Mixed Weibull Trial Simulator"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Simulation Parameters"),
      
      # Sample size and allocation
      numericInput("n", "Total Sample Size:", value = 200, min = 10, max = 10000),
      numericInput("allocation_ratio", "Treatment:Control Allocation Ratio:", value = 1, min = 0.1, max = 10, step = 0.1),
      helpText("Example: 1 = 1:1 (50% each), 2 = 2:1 (67% treatment, 33% control)"),
      
      # Progression parameters - Control
      h4("Progression Parameters - Control"),
      numericInput("prog_median_control", "Control Progression Median Time:", value = 6, min = 0.1, max = 100, step = 0.1),
      
      # Progression parameters - Treatment
      h4("Progression Parameters - Treatment"),
      numericInput("prog_median_treatment", "Treatment Progression Median Time:", value = 8, min = 0.1, max = 100, step = 0.1),
      
      # Shared progression Weibull parameters
      numericInput("nv12", "Weibull Shape (nv12) for Progression:", value = 1.5, min = 0.1, max = 10, step = 0.1),
      numericInput("lamb12", "Weibull Scale (lambda12) - Optional (0=auto):", value = 0, min = 0, max = 10, step = 0.001),
      helpText("Set to 0 to auto-calculate from median. Setting manually will change the resulting median."),
      
      # Survival parameters - Control
      h4("Survival Parameters - Control"),
      numericInput("surv_median_control", "Control Survival Median Time:", value = 12, min = 0.1, max = 100, step = 0.1),
      
      # Survival parameters - Treatment
      h4("Survival Parameters - Treatment"),
      numericInput("surv_median_treatment", "Treatment Survival Median Time:", value = 16, min = 0.1, max = 100, step = 0.1),
      
      # Shared survival Weibull parameters
      h5("Death After Progression (State 2→3)"),
      numericInput("nv23", "Weibull Shape (nv23) for Death After Progression:", value = 1.2, min = 0.1, max = 10, step = 0.1),
      numericInput("lamb23", "Weibull Scale (lambda23) - Optional (0=auto):", value = 0, min = 0, max = 10, step = 0.001),
      helpText("Set to 0 to auto-calculate from median."),
      h5("Death Without Progression (State 1→3)"),
      numericInput("nv13", "Weibull Shape (nv13) for Death Without Progression:", value = 1.5, min = 0.1, max = 10, step = 0.1),
      numericInput("lamb13", "Weibull Scale (lambda13) - Optional (0=auto):", value = 0, min = 0, max = 10, step = 0.001),
      helpText("Set to 0 to auto-calculate from median."),
      p("Survival has two pathways: death without progression (nv13) and death after progression (nv23). Overall survival median will be adjusted to match target."),
      
      # Censoring
      h4("Censoring"),
      numericInput("max_followup", "Maximum Follow-up Time:", value = 24, min = 1, max = 200, step = 1),
      
      # Action button
      actionButton("simulate", "Generate Simulation", class = "btn-primary"),
      
      br(), br(),
      h4("Notes"),
      p("This tool simulates multi-state survival data with two pathways:",
        tags$ul(
          tags$li("State 1 (At Risk) → State 2 (Progressed) with Weibull(nv12, lambda12)"),
          tags$li("Pathway 1: State 1 → State 3 (Death without progression) with Weibull(nv13, lambda13)"),
          tags$li("Pathway 2: State 1 → State 2 → State 3 (Death after progression) with progression_time + Weibull(nv23, lambda23)"),
          tags$li("Overall survival = min(death_without_prog, progression_time + death_after_prog)")
        )
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Progression-Free Survival",
                 plotOutput("pfs_plot", height = "600px"),
                 h4("Summary Statistics"),
                 verbatimTextOutput("pfs_summary")
        ),
        tabPanel("Overall Survival",
                 plotOutput("os_plot", height = "600px"),
                 h4("Summary Statistics"),
                 verbatimTextOutput("os_summary")
        ),
        tabPanel("Simulated Data",
                 h4("First 100 rows of simulated data"),
                 DT::dataTableOutput("data_table"),
                 downloadButton("download_data", "Download Full Dataset", class = "btn-default")
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Calculate Weibull scale parameter from median and shape
  # For Weibull with shape k and scale s: S(t) = exp(-(t/s)^k)
  # Median occurs when S(t) = 0.5, so: 0.5 = exp(-(median/s)^k)
  # Therefore: median = s * (log(2))^(1/k), so s = median / (log(2))^(1/k)
  
  calculate_scale_from_median <- function(median_time, shape) {
    median_time / ((log(2))^(1/shape))
  }
  
  # Calculate lambda (hazard scale parameter) from median and shape
  # For hazard h(t) = k * lambda * t^(k-1), the survival is S(t) = exp(-lambda * t^k)
  # Median: 0.5 = exp(-lambda * median^k), so lambda = log(2) / median^k
  calculate_lambda_from_median <- function(median_time, shape) {
    log(2) / (median_time^shape)
  }
  
  # Calculate hazard for transition 1->2 (progression)
  # hazd12 = nv12 * lambda12 * T^(nv12-1) * exp(modmat %*% beta12)
  calculate_hazard_12 <- function(t, nv12, lamb12, beta12 = 0, modmat = matrix(1)) {
    nv12 * lamb12 * (t^(nv12 - 1)) * exp(modmat %*% beta12)
  }
  
  # Calculate hazard for transition 2->3 (death after progression)
  calculate_hazard_23 <- function(t, nv23, lamb23, beta23 = 0, modmat = matrix(1)) {
    nv23 * lamb23 * (t^(nv23 - 1)) * exp(modmat %*% beta23)
  }
  
  # Simulate progression time using Weibull with iterative refinement
  simulate_progression <- function(n, nv12, lamb12, target_median, max_followup) {
    tolerance <- 0.1  # Tighter tolerance for better accuracy
    max_iterations <- 100  # More iterations allowed

    # Initialize scale parameter
    if (!is.na(lamb12) && lamb12 > 0) {
      # Use provided lambda to calculate scale for R's rweibull
      scale_param <- (1/lamb12)^(1/nv12)
    } else {
      # Calculate initial scale to match target median
      scale_param <- calculate_scale_from_median(target_median, nv12)
    }

    # Initialize variables
    actual_median <- target_median
    iter <- 0

    # Iterative refinement to ensure observed median is within tolerance
    for (iter in 1:max_iterations) {
      # Simulate using Weibull distribution
      prog_time <- rweibull(n, shape = nv12, scale = scale_param)

      # Apply censoring
      cens_time <- runif(n, 0, max_followup)
      prog_observed <- pmin(prog_time, cens_time)
      prog_event <- as.numeric(prog_time <= cens_time)

      # Calculate observed median using Kaplan-Meier estimator (correct method for censored data)
      fit <- survfit(Surv(prog_observed, prog_event) ~ 1)
      actual_median <- summary(fit)$table["median"]

      # Handle case where median cannot be estimated (e.g., too much censoring)
      if (is.na(actual_median) || is.null(actual_median) || actual_median == 0) {
        # Fall back to simple median as last resort
        actual_median <- median(prog_observed)
      }

      difference <- abs(actual_median - target_median)

      # Check if within tolerance
      if (difference <= tolerance) {
        cat(sprintf("Progression converged in %d iterations (median: %.2f, target: %.2f)\n",
                    iter, actual_median, target_median))
        break
      }

      # Adjust scale parameter to move median toward target
      # If actual median is too low, increase scale; if too high, decrease scale
      adjustment_factor <- target_median / actual_median

      # Prevent extreme adjustments for stability
      if (adjustment_factor > 1.5) adjustment_factor <- 1.5
      if (adjustment_factor < 0.67) adjustment_factor <- 0.67

      scale_param <- scale_param * adjustment_factor
    }

    if (iter == max_iterations) {
      cat(sprintf("Warning: Progression did not fully converge after %d iterations (median: %.2f, target: %.2f)\n",
                  max_iterations, actual_median, target_median))
    }

    # Final simulation with refined parameters
    prog_time <- rweibull(n, shape = nv12, scale = scale_param)
    cens_time <- runif(n, 0, max_followup)
    prog_observed <- pmin(prog_time, cens_time)
    prog_event <- as.numeric(prog_time <= cens_time)

    list(time = prog_observed, event = prog_event, true_time = prog_time, iterations = iter, final_scale = scale_param)
  }
  
  # Simulate survival time using multi-pathway Weibull with iterative refinement
  # Pathways:
  # 1. State 1 → State 3 (death without progression): Weibull(nv13, lambda13)
  # 2. State 1 → State 2 → State 3 (death after progression): progression_time + Weibull(nv23, lambda23)
  simulate_survival <- function(n, nv13, lamb13, nv23, lamb23, target_median, prog_times, max_followup) {
    tolerance <- 0.2  # Relaxed tolerance for survival (harder to converge due to competing risks)
    max_iterations <- 100  # More iterations allowed

    # Generate censoring times once (will be reused for consistency)
    cens_time <- runif(n, 0, max_followup)

    # Initialize scale parameters with better starting estimates
    # For death without progression (State 1 → State 3)
    if (!is.na(lamb13) && lamb13 > 0) {
      scale_param_13 <- (1/lamb13)^(1/nv13)
    } else {
      # Use target_median directly as starting point for better convergence
      scale_param_13 <- calculate_scale_from_median(target_median, nv13)
    }

    # For death after progression (State 2 → State 3)
    if (!is.na(lamb23) && lamb23 > 0) {
      scale_param_23 <- (1/lamb23)^(1/nv23)
    } else {
      # Use target_median directly as starting point
      scale_param_23 <- calculate_scale_from_median(target_median, nv23)
    }

    # Initialize variables
    actual_median <- target_median
    iter <- 0
    prev_median <- target_median

    # Iterative refinement to ensure observed median is within tolerance
    for (iter in 1:max_iterations) {
      # Simulate death without progression (State 1 → State 3)
      death_without_prog <- rweibull(n, shape = nv13, scale = scale_param_13)

      # Simulate death after progression (State 2 → State 3)
      # Time from progression to death
      death_after_prog_time <- rweibull(n, shape = nv23, scale = scale_param_23)
      # Total survival time if progression occurs first
      death_after_prog <- prog_times + death_after_prog_time

      # Overall survival is the minimum of the two pathways (competing risks)
      os_time <- pmin(death_without_prog, death_after_prog)

      # Apply censoring (using pre-generated censoring times)
      os_observed <- pmin(os_time, cens_time)
      os_event <- as.numeric(os_time <= cens_time)

      # Calculate observed median using Kaplan-Meier estimator (correct method for censored data)
      fit <- survfit(Surv(os_observed, os_event) ~ 1)
      actual_median <- summary(fit)$table["median"]

      # Handle case where median cannot be estimated (e.g., too much censoring)
      if (is.na(actual_median) || is.null(actual_median) || actual_median == 0) {
        # Fall back to simple median as last resort
        actual_median <- median(os_observed)
        if (is.na(actual_median) || actual_median == 0) {
          cat(sprintf("Warning: Cannot estimate median at iteration %d, using target\n", iter))
          actual_median <- target_median
          break
        }
      }

      difference <- abs(actual_median - target_median)

      # Check if within tolerance
      if (difference <= tolerance) {
        cat(sprintf("Survival converged in %d iterations (median: %.2f, target: %.2f)\n",
                    iter, actual_median, target_median))
        break
      }

      # Check for oscillation (median bouncing around target)
      if (iter > 10 && abs(prev_median - actual_median) < 0.05 && difference < 0.5) {
        cat(sprintf("Survival converged (oscillating) in %d iterations (median: %.2f, target: %.2f)\n",
                    iter, actual_median, target_median))
        break
      }

      # Adjust scale parameters to move median toward target
      # Adjust both parameters proportionally
      adjustment_factor <- target_median / actual_median

      # More conservative adjustments for stability (especially for survival with competing risks)
      if (adjustment_factor > 1.3) adjustment_factor <- 1.3
      if (adjustment_factor < 0.77) adjustment_factor <- 0.77

      scale_param_13 <- scale_param_13 * adjustment_factor
      scale_param_23 <- scale_param_23 * adjustment_factor

      prev_median <- actual_median
    }

    if (iter == max_iterations) {
      cat(sprintf("Warning: Survival did not fully converge after %d iterations (median: %.2f, target: %.2f, diff: %.2f)\n",
                  max_iterations, actual_median, target_median, abs(actual_median - target_median)))
    }

    # Final simulation with refined parameters
    death_without_prog <- rweibull(n, shape = nv13, scale = scale_param_13)
    death_after_prog_time <- rweibull(n, shape = nv23, scale = scale_param_23)
    death_after_prog <- prog_times + death_after_prog_time
    os_time <- pmin(death_without_prog, death_after_prog)

    os_observed <- pmin(os_time, cens_time)
    # Event = 1 if death occurred before censoring, 0 if censored
    os_event <- as.numeric(os_time <= cens_time)

    # Verify os_event was created correctly
    if (length(os_event) != n || any(is.na(os_event))) {
      stop("Error: os_event not properly generated")
    }

    list(time = os_observed, event = os_event, true_time = os_time, iterations = iter,
         final_scale_13 = scale_param_13, final_scale_23 = scale_param_23)
  }
  
  # Reactive simulation
  simulated_data <- eventReactive(input$simulate, {
    n_total <- input$n
    allocation_ratio <- input$allocation_ratio
    
    # Calculate sample sizes for each arm
    n_treatment <- round(n_total * allocation_ratio / (1 + allocation_ratio))
    n_control <- n_total - n_treatment
    
    nv12 <- input$nv12
    nv13 <- input$nv13
    nv23 <- input$nv23
    lamb12 <- if(input$lamb12 <= 0) NA else input$lamb12
    lamb13 <- if(input$lamb13 <= 0) NA else input$lamb13
    lamb23 <- if(input$lamb23 <= 0) NA else input$lamb23
    max_followup <- input$max_followup
    
    # Control arm parameters
    prog_median_control <- input$prog_median_control
    surv_median_control <- input$surv_median_control  # Total OS median (not post-progression)

    # Treatment arm parameters
    prog_median_treatment <- input$prog_median_treatment
    surv_median_treatment <- input$surv_median_treatment  # Total OS median (not post-progression)
    
    # Simulate control arm
    prog_control <- simulate_progression(n_control, nv12, lamb12, prog_median_control, max_followup)
    surv_control <- simulate_survival(n_control, nv13, lamb13, nv23, lamb23, surv_median_control, prog_control$true_time, max_followup)
    
    # Simulate treatment arm
    prog_treatment <- simulate_progression(n_treatment, nv12, lamb12, prog_median_treatment, max_followup)
    surv_treatment <- simulate_survival(n_treatment, nv13, lamb13, nv23, lamb23, surv_median_treatment, prog_treatment$true_time, max_followup)
    
    # Combine data frames
    data_control <- data.frame(
      patient_id = 1:n_control,
      arm = "Control",
      pfs_time = prog_control$time,
      pfs_event = prog_control$event,
      os_time = surv_control$time,
      os_event = surv_control$event,
      progression_time = prog_control$true_time,
      survival_time = surv_control$true_time
    )
    
    data_treatment <- data.frame(
      patient_id = (n_control + 1):n_total,
      arm = "Treatment",
      pfs_time = prog_treatment$time,
      pfs_event = prog_treatment$event,
      os_time = surv_treatment$time,
      os_event = surv_treatment$event,
      progression_time = prog_treatment$true_time,
      survival_time = surv_treatment$true_time
    )
    
    # Combine and randomize order (to simulate randomization)
    combined_data <- rbind(data_control, data_treatment)
    combined_data <- combined_data[sample(nrow(combined_data)), ]
    combined_data$patient_id <- 1:nrow(combined_data)
    
    # Verify all required columns are present
    required_cols <- c("patient_id", "arm", "pfs_time", "pfs_event", "os_time", "os_event")
    if (!all(required_cols %in% names(combined_data))) {
      missing_cols <- setdiff(required_cols, names(combined_data))
      stop(paste("Missing columns in output data:", paste(missing_cols, collapse = ", ")))
    }
    
    # Verify os_event contains valid values (0s and 1s)
    if (!all(combined_data$os_event %in% c(0, 1))) {
      stop("os_event contains invalid values (should be 0 or 1)")
    }
    
    combined_data
  })
  
  # PFS Kaplan-Meier plot
  output$pfs_plot <- renderPlot({
    data <- simulated_data()
    
    if (is.null(data) || nrow(data) == 0) {
      return(NULL)
    }
    
    # Fit KM curves by treatment arm
    fit_pfs <- survfit(Surv(pfs_time, pfs_event) ~ arm, data = data)
    
    # Calculate medians for each arm
    median_data <- data.frame(
      arm = c("Control", "Treatment"),
      median = NA
    )
    
    for (i in 1:length(fit_pfs$strata)) {
      arm_name <- names(fit_pfs$strata)[i]
      arm_level <- gsub("arm=", "", arm_name)
      arm_subset <- data[data$arm == arm_level, ]
      fit_arm <- survfit(Surv(pfs_time, pfs_event) ~ 1, data = arm_subset)
      median_val <- summary(fit_arm)$table["median"]
      median_data$median[median_data$arm == arm_level] <- ifelse(is.na(median_val), NA, median_val)
    }
    
    # Plot with separate curves by arm
    p <- ggsurvplot(
      fit_pfs,
      data = data,
      title = "Kaplan-Meier Curve: Progression-Free Survival by Treatment Arm",
      xlab = "Time",
      ylab = "Survival Probability",
      risk.table = TRUE,
      conf.int = TRUE,
      legend.title = "Treatment Arm",
      legend.labs = c("Control", "Treatment"),
      palette = c("#2196F3", "#F44336"),  # Blue (Control) and Red (Treatment) colors
      pval = TRUE,
      pval.method = TRUE
    )
    
    # Add median lines and labels for each arm
    y_positions <- list("Control" = 0.45, "Treatment" = 0.35)
    colors <- list("Control" = "#2196F3", "Treatment" = "#F44336")
    
    for (i in 1:nrow(median_data)) {
      if (!is.na(median_data$median[i])) {
        arm_name <- median_data$arm[i]
        median_val <- median_data$median[i]
        median_label <- paste0(arm_name, " Median: ", round(median_val, 2))
        
        p$plot <- p$plot + 
          geom_vline(xintercept = median_val, linetype = "dashed", 
                    color = colors[[arm_name]], linewidth = 1) +
          annotate("text", x = median_val, y = y_positions[[arm_name]], 
                  label = median_label, 
                  color = colors[[arm_name]], hjust = -0.1, size = 4, fontface = "bold")
      }
    }
    
    print(p)
  })
  
  # OS Kaplan-Meier plot
  output$os_plot <- renderPlot({
    data <- simulated_data()
    
    if (is.null(data) || nrow(data) == 0) {
      return(NULL)
    }
    
    # Fit KM curves by treatment arm
    fit_os <- survfit(Surv(os_time, os_event) ~ arm, data = data)
    
    # Calculate medians for each arm
    median_data <- data.frame(
      arm = c("Control", "Treatment"),
      median = NA
    )
    
    for (i in 1:length(fit_os$strata)) {
      arm_name <- names(fit_os$strata)[i]
      arm_level <- gsub("arm=", "", arm_name)
      arm_subset <- data[data$arm == arm_level, ]
      fit_arm <- survfit(Surv(os_time, os_event) ~ 1, data = arm_subset)
      median_val <- summary(fit_arm)$table["median"]
      median_data$median[median_data$arm == arm_level] <- ifelse(is.na(median_val), NA, median_val)
    }
    
    # Plot with separate curves by arm
    p <- ggsurvplot(
      fit_os,
      data = data,
      title = "Kaplan-Meier Curve: Overall Survival by Treatment Arm",
      xlab = "Time",
      ylab = "Survival Probability",
      risk.table = TRUE,
      conf.int = TRUE,
      legend.title = "Treatment Arm",
      legend.labs = c("Control", "Treatment"),
      palette = c("#2196F3", "#F44336"),  # Blue (Control) and Red (Treatment) colors
      pval = TRUE,
      pval.method = TRUE
    )
    
    # Add median lines and labels for each arm
    y_positions <- list("Control" = 0.45, "Treatment" = 0.35)
    colors <- list("Control" = "#2196F3", "Treatment" = "#F44336")
    
    for (i in 1:nrow(median_data)) {
      if (!is.na(median_data$median[i])) {
        arm_name <- median_data$arm[i]
        median_val <- median_data$median[i]
        median_label <- paste0(arm_name, " Median: ", round(median_val, 2))
        
        p$plot <- p$plot + 
          geom_vline(xintercept = median_val, linetype = "dashed", 
                    color = colors[[arm_name]], linewidth = 1) +
          annotate("text", x = median_val, y = y_positions[[arm_name]], 
                  label = median_label, 
                  color = colors[[arm_name]], hjust = -0.1, size = 4, fontface = "bold")
      }
    }
    
    print(p)
  })
  
  # PFS summary
  output$pfs_summary <- renderPrint({
    data <- simulated_data()
    
    if (is.null(data) || nrow(data) == 0) {
      return(cat("Click 'Generate Simulation' to create data"))
    }
    
    cat("=== Progression-Free Survival Summary by Treatment Arm ===\n\n")
    fit_pfs <- survfit(Surv(pfs_time, pfs_event) ~ arm, data = data)
    
    # Overall summary
    cat("Overall Summary:\n")
    print(summary(fit_pfs))
    
    cat("\n\n=== By Treatment Arm ===\n")
    # Summary by arm
    for (arm in unique(data$arm)) {
      cat("\n", arm, "Arm:\n")
      arm_data <- data[data$arm == arm, ]
      fit_arm <- survfit(Surv(pfs_time, pfs_event) ~ 1, data = arm_data)
      print(summary(fit_arm))
    }
  })
  
  # OS summary
  output$os_summary <- renderPrint({
    data <- simulated_data()
    
    if (is.null(data) || nrow(data) == 0) {
      return(cat("Click 'Generate Simulation' to create data"))
    }
    
    cat("=== Overall Survival Summary by Treatment Arm ===\n\n")
    fit_os <- survfit(Surv(os_time, os_event) ~ arm, data = data)
    
    # Overall summary
    cat("Overall Summary:\n")
    print(summary(fit_os))
    
    cat("\n\n=== By Treatment Arm ===\n")
    # Summary by arm
    for (arm in unique(data$arm)) {
      cat("\n", arm, "Arm:\n")
      arm_data <- data[data$arm == arm, ]
      fit_arm <- survfit(Surv(os_time, os_event) ~ 1, data = arm_data)
      print(summary(fit_arm))
    }
  })
  
  # Data table
  output$data_table <- DT::renderDataTable({
    data <- simulated_data()
    
    if (is.null(data) || nrow(data) == 0) {
      return(data.frame(message = "Click 'Generate Simulation' to create data"))
    }
    
    head(data, 100)
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # Download handler
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("simulated_trial_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(simulated_data(), file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)

