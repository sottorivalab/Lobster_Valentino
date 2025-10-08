
wrap_wide_table <- function(df, max_cols = 7) {
  num_cols <- ncol(df)
  
  if (num_cols > max_cols) {
    # If the number of columns exceeds the limit, split the data frame
    df_part1 <- df %>% select(1:max_cols)
    df_part2 <- df %>% select((max_cols + 1):num_cols)
    
    # Return a list of the two data frames
    return(list(part1 = df_part1, part2 = df_part2))
  } else {
    # If it's within the limit, return a list containing the original data frame
    return(list(full = df))
  }
}

get_driver_escape_table <- function(fit, cohort){
  
  # -------------------------------------------------
  ## TABLE OF DRIVER MUTATIONS AND THEIR SIGNIFICANCE
  
  drivers = fit %>% filter(IS_DRIVER==TRUE) %>% select(GENE, significance)
  groups <- unique(drivers$significance)
  driver_matrix <- matrix("", nrow = nrow(drivers), ncol = length(groups))
  colnames(driver_matrix) = gsub(" ", "\n", groups)
  # Populate the matrix with row names where the Group matches
  if (nrow(drivers)>0){
    for (i in 1:nrow(drivers)){
      group <- gsub(" ", "\n", as.vector(drivers$significance)[i])
      # put in first empty space
      j = which(driver_matrix[,group] == "")[1]
      driver_matrix[j, group] <- as.vector(drivers$GENE)[[i]]
    }
    # Convert the matrix to a data frame
    driver_df <- as.data.frame(driver_matrix, stringsAsFactors = FALSE)
    driver_df <- driver_df[apply(driver_df, 1, function(row) !all(row == "")),,drop=FALSE]
    
    # --- APPLY VERTICAL WRAPPING TO LONG COLUMNS ---
    for (col in names(driver_df)) {
      driver_df <- split_long_column(driver_df, col, 25)
    }
    
    # Sort the data frame columns
    driver_df <- driver_df[, order(names(driver_df)), drop = FALSE]
    
    # --- APPLY HORIZONTAL WRAPPING TO WIDE TABLES ---
    driver_df_list <- wrap_wide_table(driver_df, max_cols = 7)
    
  } else {
    driver_df_list = NULL
  }
  
  #--------------------------
  
  if (cohort!='KATS'){
    escapes = fit %>% filter(IS_ESCAPE==TRUE) %>% select(GENE, significance)
    groups <- unique(escapes$significance)
    escape_matrix <- matrix("", nrow = nrow(escapes), ncol = length(groups))
    colnames(escape_matrix) = gsub(" ", "\n", groups)
    # Populate the matrix with row names where the Group matches
    if (nrow(escapes)>0){
      for (i in 1:nrow(escapes)){
        if (cohort=='KATS'){
          group <- gsub(" ", "\n", as.vector(escapes$clusterNA)[i])
        } else {
          group <- gsub(" ", "\n", as.vector(escapes$significance)[i])
        }
        # put in first empty space
        j = which(escape_matrix[,group] == "")[1]
        escape_matrix[j, group] <- as.vector(escapes$GENE)[[i]]
      }
      # Convert the matrix to a data frame
      escape_df <- as.data.frame(escape_matrix, stringsAsFactors = FALSE)
      escape_df <- escape_df[apply(escape_df, 1, function(row) !all(row == "")),,drop=FALSE]
      
      # --- APPLY VERTICAL WRAPPING TO LONG COLUMNS ---
      for (col in names(escape_df)) {
        escape_df <- split_long_column(escape_df, col, 25)
      }
      
      # Sort the data frame columns
      escape_df <- escape_df[, order(names(escape_df)), drop = FALSE]
      
      # --- APPLY HORIZONTAL WRAPPING TO WIDE TABLES ---
      escape_df_list <- wrap_wide_table(escape_df, max_cols = 7)
      
    } else {
      escape_df_list = NULL
    }
  } else {
    escape_df_list = NULL
  }
  
  # Define the theme for all tables
  my_theme <- ttheme(
    base_size = 6, 
    colnames.style = list(fg_params = list(fontsize = 8, fontface = "bold")), 
    tbody.style = list(fg_params = list(fontsize = 6)),
    padding = unit(c(4, 2), "mm")
  )
  
  driver_plot <- NULL
  escape_plot <- NULL
  
  if (!is.null(driver_df_list)) {
    driver_tables <- lapply(driver_df_list, function(df) {
      ggtexttable(df, rows = NULL, theme = my_theme)
    })
    driver_plot <- grid.arrange(grobs = driver_tables, ncol = 1, padding=0)
  }
  
  if (!is.null(escape_df_list)) {
    escape_tables <- lapply(escape_df_list, function(df) {
      ggtexttable(df, rows = NULL, theme = my_theme)
    })
    escape_plot <- grid.arrange(grobs = escape_tables, ncol = 1, padding=0)
  }
  
  return(list(driver_table = driver_plot, escape_table = escape_plot))
}


split_long_column <- function(data, col_name, max_length) {
  
  # Extract the column data
  col_data <- data[[col_name]]
  
  # Calculate the number of required columns
  num_columns <- ceiling(sum(col_data!="") / max_length)
  
  # Check if the column is too long
  if (length(col_data) > max_length) {
    
    for (i in 1:num_columns) {
      start_index <- (i - 1) * max_length + 1
      end_index <- min(i * max_length, sum(col_data!=""))
      new_colname <- paste(col_name, sprintf("(%d)",i), sep = "\n")
      
      # Add new column to the dataframe
      if (is.null(nrow(df))){
        df = data.frame(rep("",max_length))  
        colnames(df) = new_colname
      } else {
        new_df = data.frame(rep("",max_length))
        colnames(new_df) = new_colname
        df = cbind(df,new_df)
      }
      df[[new_colname]][1:(end_index - start_index + 1)] <- col_data[start_index:end_index]
    }
    
    # rest of dataframe
    reduced_df = data %>% select(-!!sym(col_name))
    if (ncol(reduced_df)>0){
      reduced_df <- reduced_df[apply(reduced_df, 1, function(row) !all(row == "")),,drop=FALSE]
      if (nrow(reduced_df) < max_length){
        n_rows = max_length - nrow(reduced_df)
        NA_rows = as.data.frame(matrix(rep("", n_rows*ncol(reduced_df)), nrow = n_rows, ncol = ncol(reduced_df)))
        colnames(NA_rows) = colnames(reduced_df)
        reduced_df = rbind(reduced_df, NA_rows)
      } else if (nrow(reduced_df) > max_length){
        n_rows = nrow(reduced_df) - max_length
        NA_rows = matrix(rep("", n_rows*ncol(df)), nrow = n_rows, ncol = ncol(df))
        colnames(NA_rows) = colnames(df)
        df = rbind(df, NA_rows)
      }
      data = cbind(df,reduced_df)
    } else {
      data = df
    }
    
  }
  
  return(data)
}