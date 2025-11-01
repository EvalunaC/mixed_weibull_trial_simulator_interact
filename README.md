# Mixed Weibull Trial Simulator

A Shiny R application for generating simulated clinical trial data with progression-free survival (PFS) and overall survival (OS) endpoints, using Weibull mixture distributions for state transitions.

## Features

- **Simulate multi-state survival data**: Models transitions between disease states (at risk → progressed → death) using Weibull distributions
- **Customizable parameters**: 
  - Sample size
  - Target median times for progression and survival
  - Weibull shape (nv) and scale (lambda) parameters for each transition
  - Maximum follow-up time for censoring
- **Visualization**: 
  - Kaplan-Meier curves for PFS and OS
  - Median survival times labeled on plots
  - 95% confidence intervals
  - Risk tables showing number at risk over time
- **Data export**: Download simulated datasets as CSV files

## Installation

1. Ensure you have R installed (version 4.0 or higher recommended)
2. Install required R packages:

```r
install.packages(c("shiny", "survival", "survminer", "dplyr", "ggplot2", "DT"))
```

## Usage

1. Open R or RStudio
2. Navigate to the directory containing `app.R`
3. Run the application:

```r
shiny::runApp("app.R")
```

Or if you're in the directory:

```r
shiny::runApp()
```

4. In the Shiny interface:
   - Set your desired parameters in the sidebar
   - Click "Generate Simulation" to create new data
   - View the Kaplan-Meier curves in the tabs
   - Download the simulated data if needed

## Parameters

### Progression Parameters
- **Progression Median Time**: Target median for progression-free survival
- **nv12**: Weibull shape parameter for progression transition (State 1 → State 2)
- **lambda12**: Weibull scale parameter for progression transition

### Overall Survival Parameters
- **Survival Median Time**: Target median for overall survival
- **nv23**: Weibull shape parameter for death transition (State 2 → State 3)
- **lambda23**: Weibull scale parameter for death transition

### Weibull Hazard Function

The hazard function used for transitions follows the form:

```
hazd12 = nv12 * lambda12 * T^(nv12-1) * exp(modmat %*% beta12)
```

Where:
- `nv12` is the Weibull shape parameter
- `lambda12` is the Weibull scale parameter
- `T` is time
- `modmat %*% beta12` represents covariate effects (set to 0 in this simplified version)

## Output

The application generates:
- **Progression-Free Survival (PFS)**: Time from entry to disease progression or death
- **Overall Survival (OS)**: Time from entry to death from any cause

Both endpoints include:
- Observed times (with censoring applied)
- Event indicators (1 = event occurred, 0 = censored)
- True underlying event times (before censoring)

## Notes

- The simulation uses Weibull distributions to model time-to-event outcomes
- Censoring is applied uniformly across the maximum follow-up period
- The median times are matched by adjusting the Weibull scale parameters
- For a true multi-state model, survival after progression could be modeled separately; currently, OS is calculated as the maximum of progression time and survival time from entry

## License

This tool is provided for research and educational purposes.

