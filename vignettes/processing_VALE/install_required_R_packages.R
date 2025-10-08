# A list of the required packages
required_packages <- c(
  "cowplot", "ggplot2", "ggrepel", "plyr", "dplyr", "R.utils",
  "lubridate", "tidyr", "DEoptim",
  "abind", "combinat", "vcfR", "purrr", "stringr",
  "doParallel", "rtern", "fields", "RColorBrewer", "ggpubr", "gridExtra", 
  "rprojroot", "yaml"
)

# Set the CRAN mirror
options(repos = "https://cran.r-project.org")

# Function to check and install packages
install_if_missing <- function(packages) {
  # Get a list of already installed packages
  installed_packages <- installed.packages()[, "Package"]
  
  # Find packages that are not yet installed
  missing_packages <- packages[!packages %in% installed_packages]
  
  if (length(missing_packages) > 0) {
    # Install the missing packages
    install.packages(missing_packages, dependencies = TRUE)
  }
}

# The main function to install the packages
main <- function() {
  # Note: R version is not a direct argument to the install.packages() function.
  # The packages will be installed from the CRAN repository that corresponds to the R version being run.
  
  # Check and install the required packages
  install_if_missing(required_packages)
  
  print("All packages have been checked and installed if they were missing.")
}

# Run the main function
main() # make sure R version 4.2.3 is used!