# A vector containing all the paths to be created
paths <- c(
  "results/vale/vale_fits_pass/annotated",
  "results/vale/vale_fits_from_dataframe/annotated",
  "results/vale/vale_plots_pass",
  "results/vale/vale_pass_ploidy_mix",
  "results/vale/vale_pass_ploidy_mix_new/annotated",
  "results/vale/vale_pass_ploidy_mix_new/final",
  "results/vale/vale_scatterplots_pass/clinical_scatters/significance",
  "results/vale/vale_scatterplots_pass/clinical_scatters/significance_untransformed"
)

# Loop through each path
for (path in paths) {
  # Check if the directory already exists
  if (!dir.exists(path)) {
    # If it doesn't exist, create it recursively
    dir.create(path, recursive = TRUE)
    cat("Created new directory:", path, "\n")
  } else {
    # If it already exists, do nothing
    cat("Directory already exists, skipping:", path, "\n")
  }
}

cat("All specified directories have been processed.\n")