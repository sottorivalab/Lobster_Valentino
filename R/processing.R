# Functions for pre- and post-processing of lobster data for longitudinal cohorts MISSONI, MAYA and VALENTINO
#------------------------------------------------------------------------------------------------------------
# CONTENTS:
# 1) pre-processing cohort data
# 2) ploidies
# 3) manipulating vcfs
# 4) gene and driver labelling

library(rprojroot)
library(yaml)

#####################################################################
# 1) pre-processing cohort data

select_fit <- function(fit,patient_str,ploidy,cohort,user_dir){
  ploidy = as.numeric(ploidy)
  fit_select = read.csv(sprintf('%s/%s/%s_chosen_fits.csv',user_dir,tolower(cohort),tolower(cohort)),header = TRUE,sep=';')
  if (!is.na(ploidy)){
    num = fit_select %>% 
      dplyr::filter(PLOIDY==ploidy,SAMPLE==gsub("_"," ",patient_str)) %>% 
      select(chosen.fit) %>% unlist()
  } else {
    num = fit_select %>% 
      dplyr::filter(is.na(PLOIDY),SAMPLE==gsub("_"," ",patient_str)) %>% 
      select(chosen.fit) %>% unlist()
  }
  
  if (length(num)==0){
    selected_fit = fit$best
  } else if (is.na(num)){
    selected_fit = fit$best
  } else {
    selected_fit = fit$runs[[as.character(num)]]
  }
  
  return(selected_fit)
}

my_time_split <- function(df){
  # group split that returns empty dataframe in case of missing number
  group_range <- 1:4
  
  # Perform group split and create empty dataframes for missing groups
  grouped_dataframes <- lapply(group_range, function(group) {
    if (group %in% unique(df$encoded_timepoint)) {
      subset(df, df$encoded_timepoint==group)
    } else {
      data.frame(group = numeric(0), value = numeric(0))
    }
  })
  
  return (grouped_dataframes)
}

append_VAF_time <- function(patient_time,i){
  # joins together timepoints for patient i and makes df neat  
  non_empty_indices <- sapply(patient_time[[i]], function(x) nrow(x) > 0)
  result_indices <- which(non_empty_indices)

  for (j in 1:length(result_indices)){
    df = patient_time[[i]][[result_indices[j]]] %>% 
         dplyr::select(-sample,-encoded_timepoint,-patient) %>%
         dplyr::rename("VAF_{j}" := VAF) %>% 
         dplyr::rename("FILTER_{j}" := FILTER) %>% 
         dplyr::rename("ploidy_{j}" := ploidy) %>% 
         dplyr::rename("NV_{j}" := NV) %>% 
         dplyr::rename("DP_{j}" := DP) %>% 
         dplyr::mutate("t_{j}" := result_indices[j])
    if (j==1){ joined = df } else {
      joined = joined %>% dplyr::full_join(df,by=c("mappability","CHROM","POS","REF","ALT","MUTATION_TYPE","STRENGTH","GENE","IS_DRIVER","VariantType","Variant_Classification"))
    }
  }
  
  VAF_df = joined %>% dplyr::select(contains("VAF_")) %>% tibble()
  FILTER_df = joined %>% dplyr::select(contains("FILTER_")) %>% tibble()
  times_df = joined %>% dplyr::select(contains("t_") & !contains("Variant")) %>% tibble()
  ploidy_df = joined %>% dplyr::select(contains("ploidy_")) %>% tibble()
  NV_df = joined %>% dplyr::select(contains("NV_")) %>% tibble()
  DP_df = joined %>% dplyr::select(contains("DP_")) %>% tibble()
  joined = joined %>% 
    dplyr::select(-contains("VAF_")) %>%
    dplyr::select(-contains("FILTER_")) %>%
    dplyr::select(-(contains("t_") & !contains("Variant"))) %>%
    dplyr::select(-contains("ploidy_")) %>%
    dplyr::select(-contains("NV_")) %>%
    dplyr::select(-contains("DP_")) %>%
    dplyr::mutate(VAF = VAF_df) %>%
    dplyr::mutate(FILTER = FILTER_df) %>%
    dplyr::mutate(times = times_df) %>%
    dplyr::mutate(DP = DP_df) %>%
    dplyr::mutate(NV = NV_df) %>%
    dplyr::mutate(ploidy = ploidy_df)
  names(joined$VAF) = gsub(pattern = "VAF_", replacement = "", x = names(joined$VAF))
  names(joined$FILTER) = gsub(pattern = "FILTER_", replacement = "", x = names(joined$FILTER))
  names(joined$times) = gsub(pattern = "t_", replacement = "", x = names(joined$times))
  names(joined$ploidy) = gsub(pattern = "ploidy_", replacement = "", x = names(joined$ploidy))
  names(joined$NV) = gsub(pattern = "NV_", replacement = "", x = names(joined$NV))
  names(joined$DP) = gsub(pattern = "DP_", replacement = "", x = names(joined$DP))
  
  joined$times = joined$times %>% fill(names(joined$times), .direction = 'downup')

  return(joined)
}

read_maf <- function(wes_path){
  
    wes_maf = read.csv(wes_path,header = TRUE,sep = "\t", skip=1) %>% tibble::as_tibble() 
    wes_df_filtered = wes_maf %>% 
      select(Chromosome,Start_Position,Reference_Allele,Allele,Hugo_Symbol,Consequence,IMPACT,FILTER,MAPPABILITY,t_alt_count,t_ref_count,Tumor_Sample_Barcode,Variant_Type, Variant_Classification) %>%
      dplyr::rename(CHROM = Chromosome, POS = Start_Position, REF = Reference_Allele, ALT = Allele, GENE = Hugo_Symbol, NV = t_alt_count, NR = t_ref_count, sample = Tumor_Sample_Barcode, MUTATION_TYPE = Consequence, STRENGTH = IMPACT, mappability = MAPPABILITY, VariantType = Variant_Type) %>% 
      dplyr::mutate(DP = NV+NR) %>% # there are some cases in maf where depth!=NV+NR 
      dplyr::mutate(VAF = NV/DP) %>%
      dplyr::select(-NR)
    
    return(wes_df_filtered)
}

filter_maf <- function(df, filtering){
  
    if (filtering=='pass'){
      df = df %>% dplyr::filter(FILTER %in% c('PASS'))
    } else {
      df = df %>% dplyr::filter(FILTER %in% c('PASS','strand_bias'))
    }
    df = df %>% dplyr::filter(mappability>=0.6)
  
  return(df)
}


return_config_dict <- function(project_root=NULL){
  if (is.null(project_root)){
    project_root <- rprojroot::find_root_file(criterion = rprojroot::has_file("DESCRIPTION"), path = ".")
    print(sprintf("%s detected as project root",project_root))
  }
  config <- yaml::read_yaml(sprintf("%s/user_config.yml",project_root))
  
  return(config)
}


format_vale_wes <- function(i, filtering){

  project_root <- rprojroot::find_root_file(criterion = rprojroot::has_file("DESCRIPTION"), path = ".")
  print(sprintf("%s detected as project root",project_root))
  vale_df = read.csv(sprintf('%s/vale/VALENTINO_samplesheet.reduced.csv',project_root))
  sample_names = vale_df$sample
  patient_names = vale_df$case
  
  samples = sample_names[patient_names==(unique(patient_names)[i])]
  patient = unique(patient_names)[i]
  print(samples)
  
  config = return_config_dict(project_root = project_root)
  maf_dir <- config$paths$maf_dir
  print(maf_dir)
  
  for (j in 1:length(samples)){
    my_sample = samples[j]
    print(my_sample)
    wes_path = sprintf("%s/%s.maf",maf_dir,my_sample)
    
    encoded_timepoint = vale_df %>% filter(sample==my_sample) %>% select(encoded_timepoint) %>% unlist()
    print(encoded_timepoint)
    
    wes_df = read_maf(wes_path)
    wes_df_filtered = filter_maf(wes_df, filtering)
    
    wes_df_filtered$encoded_timepoint = rep(encoded_timepoint,nrow(wes_df_filtered))
    
    cat(sprintf('%s: %d\n',my_sample,nrow(wes_df_filtered)), file = sprintf("%s/vale/VALE_mut_nums.txt", project_root), append = TRUE)
    
    if (j==1){ full_df = wes_df_filtered } else { full_df = rbind(full_df,wes_df_filtered)}     
  } 
  
  full_df <- MC_GEL_driver_labelling_posthoc(full_df) # also appropriate for vale
  # mask to NA all mutations with VAF of 0
  full_df$VAF[full_df$VAF==0] = NA
  
  full_df$patient <- gsub("_(BL|BR|PD)$", "", full_df$sample)
  full_df = append_mutation_ploidies_vale(full_df)
  
  patient_split <- full_df %>% group_split(patient)
  double_split <- lapply(patient_split,my_time_split)
  
  return(list(patient_split=patient_split, patient_time_split=double_split))
}

#####################################################################
# 2) ploidies

split_ploidies <- function(processed_list){
  
  set.seed(1)
  patient_time = processed_list$patient_time_split

  df = append_VAF_time(patient_time,1)
  
  # create diff dataframes of max 5000 mutations for different and constant ploidies
	df_ploidy_NA = data.frame() 
	df_ploidy = list()
	
  if (!all(is.na(df$ploidy))){
    # NA ploidies (out of sequence bounds)
    na_rows_idx <- which(apply(df$ploidy, 1, function(row) {all(is.na(row))}))
    na_rows = df[na_rows_idx,]
    if (nrow(na_rows)>0){
      df_ploidy[['ploidy_mix']] = na_rows
    }
    
    # ploidy changes between timepoints
    different_rows_idx <- which(apply(df$ploidy, 1, function(row) {length(unique(row))>1}))
    different_rows = df[different_rows_idx,]
    if (nrow(different_rows)>0){
      if (!is.null(df_ploidy[['ploidy_mix']])){
        df_ploidy[['ploidy_mix']] = rbind(df_ploidy[['ploidy_mix']], different_rows)
      } else {
        df_ploidy[['ploidy_mix']] = different_rows
      }
    }
    
    # ploidy doesn't change 
    existing_ploidies = unique(unlist(df$ploidy))        
    for (ploidy in existing_ploidies){
        df_ploid_idx   = which(apply(df$ploidy, 1, function(row) all(row == ploidy))); # for fit
        df_ploid = df[df_ploid_idx,]
        ploidy_str = sprintf("ploidy_%d",ploidy)
        
        # proposal dataframe (dataframe that would be fed into lobster fitting)
        proposal = df_ploid[df_ploid$VariantType=='SNP',]
        coverage_limit = 10
        proposal = proposal[!apply(proposal$DP, 1, function(row) any(row < coverage_limit)),]
        
        if (nrow(proposal)>=10){
          df_ploid <- df_ploid[sample(nrow(df_ploid), min(5000,nrow(df_ploid)), replace = FALSE), ]
          df_ploidy[[ploidy_str]] = df_ploid
        } else { # if not 10 mutations, don't create fit for this ploidy - add to mixed ploidies df
          if (!is.null(df_ploidy[['ploidy_mix']])){
            df_ploidy[['ploidy_mix']] = rbind(df_ploidy[['ploidy_mix']], df_ploid)
          } else {
            df_ploidy[['ploidy_mix']] = df_ploid
          }
        }
        
        print(sprintf('nrow(%s): %d',ploidy_str,nrow(df_ploidy[[ploidy_str]])))
    }
    
    print(sprintf('nrow(ploidy_mix): %d',nrow(df_ploidy[['ploidy_mix']])))
    
  } else {
    df_ploidy_NA = df[sample(nrow(df), min(5000,nrow(df)), replace = FALSE), ]
    df_ploidy[['ploidy_NA']] = df_ploidy_NA
  }
	
	total_df = do.call(rbind, df_ploidy)
	print(sprintf('total rows at end: %d',nrow(total_df)))
	missing_rows = dplyr::anti_join(total_df,df)
	
  return(df_ploidy)
}

append_mutation_ploidies_vale <- function(df){
  patient = unique(df$patient)
  samples = unique(df$sample)
  
  df$ploidy = rep(NA, nrow(df))
  df$POS = df$POS %>% as.integer()
  
  config = return_config_dict()
  seg_path = config$paths$ichor_seg_dir
    
  if (patient!='VALE_38'){
            
    for (j in 1:length(unique(samples))){
    
      sample = unique(samples)[j]
      file_str = sprintf('%s/%s.seg.txt',seg_path,sample)
      cn_df = read.csv(file_str,sep='\t') %>% select(ID, chrom, start, end, Corrected_Copy_Number) %>% rename(sample = ID, CHROM = chrom)
      counts <- table(cn_df$CHROM)
      uniquely_appearing_values <- names(counts[counts == 1])
      not_unique = names(counts)[!(names(counts) %in% uniquely_appearing_values)]
      
      # if appears in file uniquely then it applies to the whole chromosome
      for (num in uniquely_appearing_values){
        tryCatch({
          df[(df$sample==sample)&(df$CHROM==sprintf("chr%s",num)),]$ploidy = cn_df %>% filter(CHROM==num) %>% select(Corrected_Copy_Number) %>% unlist()   
        }, error = function(e){})
      } # apply separately with position if not unique
      for (num in not_unique){
        lines = cn_df %>% filter(CHROM==num) 
        for (lineno in 1:nrow(lines)){
          line = lines[lineno,]
          start = line %>% select(start) %>% unlist()
          end =  line %>% select(end) %>% unlist()
          print(line %>% select(Corrected_Copy_Number) %>% unlist())
          tryCatch({
            df[(df$sample==sample)&(df$CHROM==sprintf("chr%s",num))&
                 (df$POS<=end)&(df$POS>=start),]$ploidy = line %>% select(Corrected_Copy_Number) %>% unlist()   
          }, error = function(e){})
        }
      }
    }
  }
  # exclude X and Y chromosomes
  df = df[!(df$CHROM %in% c('chrX','chrY')),]
  
  return(df)
} 

#####################################################################
# 3) manipulating vcfs

select_first_element <- function(vec) {
  matching_elements <- vec[grep("^CSQ", vec)]
  if (length(matching_elements) > 0) {
    return(matching_elements[1])
  } else {
    return(NA)
  }
}

extractGenotypeData <- function(vcf) {
  gt_all <- extract_gt_tidy(vcf)
  gt_field = gt_all %>% 
    tidyr::separate(gt_AD, sep = ',', into = c("NR", "NV")) %>%
    dplyr::mutate(NR = as.numeric(NR), NV = as.numeric(NV), DP = NV + NR, VAF = NV/DP) %>%
    dplyr::rename(sample = Indiv) 
  na_muts <- gt_field %>% dplyr::filter(is.na(VAF))
  if(nrow(na_muts) > 0) {
    cat("There are n =", nrow(na_muts), "mutations with VAF that is NA; they will NOT be removed.\n")
  }
  gt_field
}

#####################################################################
# 4) gene and driver labelling

MC_GEL_driver_labelling_posthoc <- function(df){
  config = return_config_dict()
  user_dir = config$paths$user_dir
  MC_GEL_genelist = read.csv(sprintf("%s/drivers_martincor_gel.txt",user_dir), header=FALSE) %>% unlist()
  escape_genes = readLines(sprintf('%s/escapeGenesLuizZapata_renamed.txt',user_dir))
  
  # must be in gene list, must not be synonymous or intron variant
  df$IS_DRIVER = (df$GENE %in% MC_GEL_genelist)&(!(df$Variant_Classification %in% c('Intron','Silent', "5'Flank", "3'Flank", 'Splice_Region')))
  df$IS_ESCAPE = (df$GENE %in% escape_genes)&(!(df$Variant_Classification %in% c('Intron','Silent', "5'Flank", "3'Flank", 'Splice_Region')))
  
  return(df)
}
