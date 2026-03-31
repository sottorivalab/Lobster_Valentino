library(dplyr)
library(purrr)
library(mclust)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork) 
library(ggdist)
library(grid) 

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
  path = "/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/sampling_simulations/dataframes_timeseries/", 
  pattern = "\\.Rds$", 
  full.names = TRUE
)

final_table <- calculate_strict_ari_batch(
  file_list = my_files,
  col_reference = "Type",
  cols_to_compare = c("lobster_cluster")
)

wide_view_timeseries <- final_table %>%
  select(Filename, Comparison, ARI) %>%
  pivot_wider(names_from = Comparison, values_from = ARI)

mean_lobster_ari = mean(wide_view_timeseries$`lobster_cluster vs Type`)
median_lobster_ari = median(wide_view_timeseries$`lobster_cluster vs Type`)

df <- wide_view_timeseries %>%
  tidyr::extract(
    col = Filename,
    into = c("subclonal_ccf", "purity1", "purity2", "purity3"),
    regex = "subclonal_ccf_([0-9.]+)_purity1_([0-9.]+)_purity2_([0-9.]+)_purity3_([0-9.]+).Rds",
    remove = FALSE,  
    convert = TRUE   
  ) %>%
  dplyr::mutate(
    mean_purity = (purity1 + purity2 + purity3) / 3,
    geom_mean_purity = (purity1 * purity2 * purity3)^(1/3),
    min_purity = pmin(purity1, purity2, purity3)
  ) 

plot_df <- df %>%
  mutate(subclonal_ccf = as.factor(subclonal_ccf)) %>% 
  dplyr::rename(lobster_cluster_vs_Type = 'lobster_cluster vs Type') %>% 
  pivot_longer(
    cols = c(geom_mean_purity, min_purity),  
    names_to = "purity_metric",        
    values_to = "purity_value"         
  ) %>%
  mutate(
    purity_metric = case_when(
      purity_metric == "geom_mean_purity" ~ "A. Geometric mean purity across timepoints",
      purity_metric == "min_purity" ~ "B. Minimum purity across timepoints"
    ),
    subclonal_ccf = if_else(subclonal_ccf == 0.0, "Clonal", as.character(subclonal_ccf))
  ) 

shared_x_limits <- c(0.0, 0.65)

plot_df_min = plot_df %>% 
  dplyr::filter(grepl("Min", purity_metric))
plot_df_geom_mean = plot_df %>% 
  dplyr::filter(grepl("Geometric", purity_metric))

plot_A = ggplot(plot_df_geom_mean, aes(x = purity_value, 
                         y = lobster_cluster_vs_Type, 
                         color = subclonal_ccf)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(aes(fill = subclonal_ccf), fun.data = mean_se, geom = "ribbon", alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.5) + 
  ylim(0, 1) +
  scale_x_continuous(limits = shared_x_limits) +
  theme_minimal() +
  labs(
    x = "Purity", y = NULL,
    color = "True Subclonal CCF"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold", hjust = 0), 
    panel.spacing = unit(1, "lines") 
  )
print(plot_A)

plot_B = ggplot(plot_df_min, aes(x = purity_value, 
                         y = lobster_cluster_vs_Type, 
                         color = subclonal_ccf)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(aes(fill = subclonal_ccf), fun.data = mean_se, geom = "ribbon", alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.5) + 
  ylim(0, 1) +
  scale_x_continuous(limits = shared_x_limits) +
  theme_minimal() +
  labs(
    x = "Purity", y = NULL,
    color = "True Subclonal CCF"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold", hjust = 0), 
    panel.spacing = unit(1, "lines") 
  )
print(plot_B)

vale_purities = read.csv('/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/results/vale/vale_sample_purity.csv')
vale_purities = vale_purities %>%
  dplyr::mutate(
    mean_purity = (purity.1 + purity.2 + purity.3) / 3,
    geom_mean_purity = (purity.1 * purity.2 * purity.3)^(1/3),
    min_purity = pmin(purity.1, purity.2, purity.3)
  ) 

violin_A <- ggplot(vale_purities, aes(x = geom_mean_purity)) + 
  
  geom_density(fill = "grey70", color = NA, alpha = 0.6) +
  
  geom_boxplot(aes(y = 0.5), width = 1, fill = "white", outlier.shape = NA) + 
  
  scale_x_continuous(limits = shared_x_limits) +
  labs(title = "Geometric mean purity across timepoints",
       subtitle = "Distribution of geometric mean purity\nacross timepoints of VALE patients") +
  theme_void() + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.margin = margin(b = -10),
    plot.subtitle = element_text(size = 8, face = "italic", color = "grey30", hjust = 1)
  ) 
  

violin_B <- ggplot(vale_purities, aes(x = min_purity)) + 
  
  geom_density(fill = "grey70", color = NA, alpha = 0.6) +
  
  geom_boxplot(aes(y = 0.5), width = 1, fill = "white", outlier.shape = NA) + 
  
  scale_x_continuous(limits = shared_x_limits) +
  labs(title = "Minimum purity across timepoints",
       subtitle = "Distribution of minimum purity\nacross timepoints of VALE patients") +
  theme_void() + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.margin = margin(b = -10),
    plot.subtitle = element_text(size = 8, face = "italic", color = "grey30", hjust = 1)
  )

main_column <- violin_A / 
  plot_A / 
  violin_B / 
  plot_B + 
  
  plot_layout(
    heights = c(1, 5, 1, 5), 
    guides = "collect"       
  )

global_y <- textGrob("ARI (Lobster cluster vs. simulated ground truth)", 
                     rot = 90, 
                     gp = gpar(fontsize = 14))

final_combined_plot <- wrap_elements(global_y) + main_column + 
  plot_layout(
    widths = c(1, 30),
    guides = "collect" 
  ) + 
  plot_annotation(
    title = "Longitudinal fitting",
    caption = "Adjusted Rand Index (ARI): 1.0 - perfect agreement, 0.0 - random clustering or clustered to one component\nGrey distributions represent the density of the purity metric for Valentino patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom" 
    )
  ) & 
  theme(legend.position = "bottom") 

print(final_combined_plot)

ggsave('/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/timeseries_simulation_ARI.pdf',final_combined_plot, height = 8, width=8)

library(dplyr)
library(purrr)
library(mclust)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork) 
library(ggdist)
library(grid) 

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
  
  # mobster fit failed for purities<=0.1
  valid_files <- results %>%
    dplyr::filter(Status == "Success") %>% 
    group_by(Filename) %>%
    dplyr::summarise(success_count = n()) %>%
    filter(success_count == length(cols_to_compare)) %>%
    pull(Filename)
  
  final_table <- results %>%
    filter(Filename %in% valid_files)
  
  message(sprintf("Found %d files. Retained %d files where all %d methods returned results.",  
                  length(unique(results$Filename)), length(valid_files), length(cols_to_compare)))
  
  return(final_table)
}

my_files <- list.files(
  path = "/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/sampling_simulations/dataframes_timeseries/", 
  pattern = "\\.Rds$", 
  full.names = TRUE
)

final_table <- calculate_strict_ari_batch(
  file_list = my_files,
  col_reference = "Type",
  cols_to_compare = c("lobster_cluster")
)

wide_view_timeseries <- final_table %>%
  select(Filename, Comparison, ARI) %>%
  pivot_wider(names_from = Comparison, values_from = ARI)

mean_lobster_ari = mean(wide_view_timeseries$`lobster_cluster vs Type`)
median_lobster_ari = median(wide_view_timeseries$`lobster_cluster vs Type`)

df <- wide_view_timeseries %>%
  tidyr::extract(
    col = Filename,
    into = c("subclonal_ccf", "purity1", "purity2", "purity3"),
    regex = "subclonal_ccf_([0-9.]+)_purity1_([0-9.]+)_purity2_([0-9.]+)_purity3_([0-9.]+).Rds",
    remove = FALSE,  
    convert = TRUE   
  ) %>%
  dplyr::mutate(
    mean_purity = (purity1 + purity2 + purity3) / 3,
    geom_mean_purity = (purity1 * purity2 * purity3)^(1/3),
    min_purity = pmin(purity1, purity2, purity3)
  ) 

plot_df <- df %>%
  mutate(subclonal_ccf = as.factor(subclonal_ccf)) %>% 
  dplyr::rename(lobster_cluster_vs_Type = 'lobster_cluster vs Type') %>% 
  pivot_longer(
    cols = c(geom_mean_purity, min_purity),  
    names_to = "purity_metric",        
    values_to = "purity_value"         
  ) %>%
  mutate(
    purity_metric = case_when(
      purity_metric == "geom_mean_purity" ~ "A. Geometric mean purity across timepoints",
      purity_metric == "min_purity" ~ "B. Minimum purity across timepoints"
    ),
    subclonal_ccf = if_else(subclonal_ccf == 0.0, "Clonal", as.character(subclonal_ccf))
  ) 

shared_x_limits <- c(0.0, 0.65)

plot_df_min = plot_df %>% 
  dplyr::filter(grepl("Min", purity_metric))
plot_df_geom_mean = plot_df %>% 
  dplyr::filter(grepl("Geometric", purity_metric))

plot_A = ggplot(plot_df_geom_mean, aes(x = purity_value, 
                         y = lobster_cluster_vs_Type, 
                         color = subclonal_ccf)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(aes(fill = subclonal_ccf), fun.data = mean_se, geom = "ribbon", alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.5) + 
  ylim(0, 1) +
  scale_x_continuous(limits = shared_x_limits) +
  theme_minimal() +
  labs(
    x = "Purity", y = NULL,
    color = "True Subclonal CCF"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold", hjust = 0), 
    panel.spacing = unit(1, "lines") 
  )
print(plot_A)

plot_B = ggplot(plot_df_min, aes(x = purity_value, 
                         y = lobster_cluster_vs_Type, 
                         color = subclonal_ccf)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) + 
  stat_summary(aes(fill = subclonal_ccf), fun.data = mean_se, geom = "ribbon", alpha = 0.4, color = NA, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.5) + 
  ylim(0, 1) +
  scale_x_continuous(limits = shared_x_limits) +
  theme_minimal() +
  labs(
    x = "Purity", y = NULL,
    color = "True Subclonal CCF"
  ) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold", hjust = 0), 
    panel.spacing = unit(1, "lines") 
  )
print(plot_B)

vale_purities = read.csv('/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/results/vale/vale_sample_purity.csv')
vale_purities = vale_purities %>%
  dplyr::mutate(
    mean_purity = (purity.1 + purity.2 + purity.3) / 3,
    geom_mean_purity = (purity.1 * purity.2 * purity.3)^(1/3),
    min_purity = pmin(purity.1, purity.2, purity.3)
  ) 

violin_A <- ggplot(vale_purities, aes(x = geom_mean_purity)) + 
  
  geom_density(fill = "grey70", color = NA, alpha = 0.6) +
  
  geom_boxplot(aes(y = 0.5), width = 1, fill = "white", outlier.shape = NA) + 
  
  scale_x_continuous(limits = shared_x_limits) +
  labs(title = "Geometric mean purity across timepoints",
       subtitle = "Distribution of geometric mean purity\nacross timepoints of VALE patients") +
  theme_void() + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.margin = margin(b = -10),
    plot.subtitle = element_text(size = 8, face = "italic", color = "grey30", hjust = 1)
  ) 
  

violin_B <- ggplot(vale_purities, aes(x = min_purity)) + 
  
  geom_density(fill = "grey70", color = NA, alpha = 0.6) +
  
  geom_boxplot(aes(y = 0.5), width = 1, fill = "white", outlier.shape = NA) + 
  
  scale_x_continuous(limits = shared_x_limits) +
  labs(title = "Minimum purity across timepoints",
       subtitle = "Distribution of minimum purity\nacross timepoints of VALE patients") +
  theme_void() + 
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    plot.margin = margin(b = -10),
    plot.subtitle = element_text(size = 8, face = "italic", color = "grey30", hjust = 1)
  )

col_1 <- (violin_A / plot_A) + plot_layout(heights = c(1, 5))
col_2 <- (violin_B / plot_B) + plot_layout(heights = c(1, 5))

main_grid <- col_1 | col_2

global_y <- textGrob("ARI (Lobster cluster vs.\nsimulated ground truth)", 
                     rot = 90, 
                     gp = gpar(fontsize = 14))

final_combined_plot <- wrap_elements(global_y) + main_grid + 
  plot_layout(
    widths = c(2, 40), 
    guides = "collect" 
  ) + 
  plot_annotation(
    title = "Longitudinal fitting",
    caption = "Adjusted Rand Index (ARI): 1.0 - perfect agreement, 0.0 - random clustering or clustered to one component\nGrey distributions represent the density of the purity metric for Valentino patients",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "bottom"
    )
  ) & 
  theme(legend.position = "bottom")

print(final_combined_plot)

ggsave('/scratch/madeleine.dale/lobster_paper_fresh/lobster4publishing.VALE/vignettes/fitting_validation/timeseries_simulation_ARI.pdf', final_combined_plot, height = 7, width = 12)
