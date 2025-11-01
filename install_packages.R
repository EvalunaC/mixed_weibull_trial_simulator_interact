# Script to install required packages for Mixed Weibull Trial Simulator

# Check if packages are installed, and install if missing
required_packages <- c("shiny", "survival", "survminer", "dplyr", "ggplot2", "DT")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("All required packages are installed!\n")
cat("You can now run the Shiny app with: shiny::runApp('app.R')\n")

