library(dplyr)
library(purrr)
library(mclust)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

calculate_strict_ari_batch <- function(file_list, col_reference, cols_to_compare) {
  
  results <- purrr::map_dfr(file_list, function(file_path) {
    
    tryCatch({
      if (grepl("\\.Rds$", file_path, ignore.case = TRUE)) {
        data <- readRDS(file_path)
      } else {
        data <- read_csv(file_path, show_col_types = FALSE)
      }
      
      if (!col_reference %in% names(data)) {
        return(tibble(
          Filename = basename(file_path),
          Comparison = "Reference Missing",
          ARI = NA_real_,
          Status = paste("Error: Ref column", col_reference, "not found")
        ))
      }
      
      file_results <- purrr::map_dfr(cols_to_compare, function(current_col) {
        
        if (!current_col %in% names(data)) {
          return(tibble(
            Filename = basename(file_path),
            Comparison = paste(current_col, "vs", col_reference),
            ARI = NA_real_,
            Status = paste("Error: Column", current_col, "not found")
          ))
        }
        
        ari_score <- mclust::adjustedRandIndex(data[[col_reference]], data[[current_col]])
        
        return(tibble(
          Filename = basename(file_path),
          Comparison = paste(current_col, "vs", col_reference),
          ARI = ari_score,
          Status = "Success"
        ))
      })
      
      return(file_results)
      
    }, error = function(e) {
      return(tibble(
        Filename = basename(file_path),
        Comparison = "File Error",
        ARI = NA_real_,
        Status = paste("Error reading file:", e$message)
      ))
    })
  })
  
  valid_files <- results %>%
    filter(Status == "Success") %>%
    group_by(Filename) %>%
    dplyr::summarise(success_count = n()) %>%
    filter(success_count == length(cols_to_compare)) %>%
    pull(Filename)
  
  final_table <- results %>%
    filter(Filename %in% valid_files)
  
  message(sprintf("Found %d files. Retained %d files where all %d methods succeeded.", 
                  length(unique(results$Filename)), length(valid_files), length(cols_to_compare)))
  
  return(final_table)
}

my_files <- list.files(
  path = "/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/sampling_simulations/dataframes", 
  pattern = "\\.Rds$", 
  full.names = TRUE
)

final_table <- calculate_strict_ari_batch(
  file_list = my_files,
  col_reference = "Type",
  cols_to_compare = c("mobster_cluster", "lobster_cluster")
)

final_table <- final_table %>% 
  filter(!str_detect(Filename, "0.05"))

wide_view <- final_table %>%
  select(Filename, Comparison, ARI) %>%
  pivot_wider(names_from = Comparison, values_from = ARI) %>%
  tidyr::extract(
    col = Filename,
    into = c("subclonal_ccf", "purity"),
    regex = "^([A-Za-z_0-9.]*)_purity_([0-9.]+).Rds",
    remove = FALSE,  
    convert = TRUE   
  ) %>% 
  dplyr::mutate(subclonal_ccf = str_remove(subclonal_ccf, "subclone_ccf_")) %>%
  dplyr::rename("Mobster cluster vs. simulated ground truth" = "mobster_cluster vs Type",
                "Lobster cluster vs. simulated ground truth" = "lobster_cluster vs Type")

mean_mobster_ari = mean(wide_view$`Mobster cluster vs. simulated ground truth`)
mean_lobster_ari = mean(wide_view$`Lobster cluster vs. simulated ground truth`)

median_mobster_ari = median(wide_view$`Mobster cluster vs. simulated ground truth`)
median_lobster_ari = median(wide_view$`Lobster cluster vs. simulated ground truth`)



plot_df <- wide_view %>%
  pivot_longer(
    cols = c("Mobster cluster vs. simulated ground truth", "Lobster cluster vs. simulated ground truth"), 
    names_to = "Method",       
    values_to = "ARI"          
  )

text_box_df <- data.frame(
  subclonal_ccf = 'clone',   
  purity = 1.15,        
  ARI = 0.5,            
  label = "Adjusted Rand Index (ARI):\n1.0 - Perfect agreement\n0.0 - Random clustering\n A score of 0.0 either indicates that the\nfitting failed or that it fit to one cluster\nMobster fitting failed for each instance\nwhere purity<=0.1\nMedian Lobster ARI: 0.58\nMedian Mobster ARI: 0.46"
)

purity_breaks <- sort(unique(plot_df$purity))
my_plot = ggplot(plot_df, aes(x = purity, 
                    y = ARI, 
                    color = Method, 
                    group = Method)) + 
  
  geom_text(data = text_box_df, 
            aes(x = purity, y = ARI, label = label), 
            inherit.aes = FALSE, 
            hjust = 0,          
            color = "grey30", 
            size = 3.5,
            lineheight = 1.2) + 
  
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  
  geom_point(size = 3, alpha = 0.7) + 
  geom_line(linewidth = 1) + 
  scale_x_continuous(breaks = purity_breaks) +
  
  facet_wrap(~ subclonal_ccf, labeller = as_labeller(function(x) paste("Subclonal CCF:", x))) + 
  
  theme_minimal() +
  labs(
    title = 'Single timepoint results',
    x = "Purity",
    y = "ARI Score",
    color = "Scoring Method"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12) 
  )

print(my_plot)
ggsave('/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/singlesample_simulation_ARI.pdf',my_plot, height = 10, width=10)



