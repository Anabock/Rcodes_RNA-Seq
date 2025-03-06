#########################

####This script automates the following processes involved in RNA-Seq data processing:
  ###Removal of low read counts (Default cut off: >= 10 reads) 
  ###Removal of astrocytic and microglial genes for NuTrap data processing only
  ###Batch correction for both the RNA-Seq count matrix and the removed genes count matrix, as an internal control for batch correction
  ###Creates PCA plots for both RNA-Seq corrected matrix and astrocytic and microglial genes matrix, for comparison


#######################

##Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(sva)
  library(tidyr)
  library(ggplot2)
  library(ggforce)
  library(ggfortify)
  library(glue)
  library(viridis)
  library(RColorBrewer)
})


#######################

##Crucial steps:
  #1.Load your matrix to the environment as "My_data"
    ###If CBDS' HPC/RNA-Seq pipeline was used to align and quantify the data, this input is likely named "gene_names_counts.csv"
    My_data <- as.data.frame(gene_names_counts)

  #2.Load the reference lists used for enrichment 
    ### The code will only use the lists for "Astrocytic" and "Microglial" genes here
    astrocytic <- data.frame(gene_symbol = Table_S2_Enrichment_genes_lists$Astrocytic)
      astrocytic <- astrocytic[!is.na(astrocytic$gene_symbol), , drop = FALSE]
        rownames(astrocytic) <- astrocytic$gene_symbol
      
    microglial <- data.frame(gene_symbol = Table_S2_Enrichment_genes_lists$Microglial)
      rownames(microglial) <- microglial$gene_symbol
    
  #3.Load your batch correction project sheet as "project_sheet"
    ###This will be used in "Script 4" and "Script 5"
    project_sheet <- as.data.frame(Proj_sheet)
      

##################################################################################################################################################

##Script 1:
  ###This script removes low count reads from an RNA-Seq count matrix, after genome alignment and quantification 
    #(Default cut off: >= 10 reads)
 
         
#######################
      
###Often times, the matrix has repeated gene names which will keep the filter from running 
###This first step prevents the issue by eliminating the repeated genes with the lowest counts
      
#Inspect for repeated elements and eliminate duplicated genes
  #Detect genes and list them
  duplicated_genes <- duplicated(My_data$gene_symbol) | duplicated(My_data$gene_symbol, fromLast = TRUE) #Identify duplicated gene names
    duplicated_rows <- My_data[duplicated_genes, ] #Extract duplicated gene names and corresponding rows
      duplicated_gene_names <- unique(duplicated_rows$gene_symbol) #List the duplicated gene names
      
  #Identify the sample columns (starting from the 3rd column)
  sample_columns <- colnames(My_data)[3:ncol(My_data)]
      
  #For each duplicated gene, sum the reads across all samples and retain the row with the highest sum
  My_data_dedup <- My_data %>%
    group_by(gene_symbol) %>%
      filter(rowSums(across(all_of(sample_columns))) == max(rowSums(across(all_of(sample_columns))))) %>%
        ungroup()
      
    
#######################    

###This is the actual low count removal
  #(Default cut off: >= 10 reads)      

#Remove low counts     
  #Calculate the number of sample columns starting from column 3
  num_samples <- length(sample_columns)
      
  #Create a logical vector indicating rows with counts lower than your defined cut off
  rows_to_keep <- rowSums(My_data_dedup[, sample_columns] >= 10) == num_samples  # <------------- CHOOSE CUT OFF VALUE HERE
      
  #Subset the count matrix to keep only rows with counts >= defined cut off
  My_data_clean <- My_data_dedup[rows_to_keep, ]
    
#save the data in a csv file
    write.csv(My_data_clean, "gene_counts_clean.csv", row.names = FALSE) #rename your My_data to include the selected filter
    #Double check that ", row.names = FALSE" is written on your script, or the file will have an extra column.
        

##################################################################################################################################################

##Script 2:
  ###This script screens Script 1's output "My_data_clean" and removes from it commonly known astrocytic and microglial genes
  ###This cleaning step was implemented in the NuTRAP analysis to prevent detection of gene expression changes not related to cell-specific analysis  

#Search astrocytic genes in your data set
My_data_clean <- as.data.frame(My_data_clean)
  rownames(My_data_clean) <- My_data_clean$gene_symbol
    found_astrocytic <- My_data_clean[ My_data_clean$gene_symbol %in% astrocytic$gene_symbol, ]
      found_astrocytic <- found_astrocytic[,c("gene_symbol")]

#Find the indices of rows to keep
rows_to_keep <- !(rownames(My_data_clean) %in% rownames(astrocytic))
  rows_out <- rownames(My_data_clean) %in% rownames(astrocytic)

#Subset the count My_data to keep only the rows (genes) that are not in the gene list
astrocytic_filtered <- My_data_clean[rows_to_keep, ]

#confirm that the genes are removed from your data
found_astrocytic_after_clean <- astrocytic_filtered[ astrocytic_filtered$gene_symbol %in% astrocytic$gene_symbol, ] 
###If the code ran properly, this should output "0 obs"
  

#######################

#Search microglial genes in your data set
found_microglial  <- astrocytic_filtered[ astrocytic_filtered$gene_symbol %in% microglial$gene_symbol, ]
  found_microglial  <- found_microglial[,c("gene_symbol")]

#Find the indices of rows to keep
rows_to_keep <- !(rownames(astrocytic_filtered) %in% rownames(microglial))
  rows_out <- rownames(astrocytic_filtered) %in% rownames(microglial)

#Subset the count My_data to keep only the rows (genes) that are not in the gene list
astr_micr_filtered <- astrocytic_filtered[rows_to_keep, ]

#Confirm that the genes are removed from your data
found_microglial_after_clean <- astr_micr_filtered[ astr_micr_filtered$gene_symbol %in% microglial$gene_symbol, ]
###If the code ran properly, this should output "0 obs"

#save the data in a csv file
    write.csv(astr_micr_filtered,"filtered_matrix.csv", row.names = FALSE)
      

##################################################################################################################################################

##Script 3:
  ###This script screens Script 1's output "My_data_clean" and creates a matrix with astrocytic and microglial genes only 
  ###The output matrix in this code will be subjected to the same batch corrections as the "filtered_matrix" and then plotted with PCA
  ###This step works as a control for the batch correction of the "filtered_matrix"

#Search astrocytic genes in your dataset
found_astrocytic2 <- My_data_clean[ My_data_clean$gene_symbol %in% astrocytic$gene_symbol, ]

#Search microglial genes in your dataset
found_microglial2 <- My_data_clean[ My_data_clean$gene_symbol %in% microglial$gene_symbol, ]

#Check the column names to ensure they match
if (!all(colnames(found_astrocytic2) == colnames(found_microglial2))) {
  stop("The sample names (column names) in the matrices do not match.")
}

#Combine the matrices by row binding (genes)
astr_micr_matrix <- bind_rows(found_astrocytic2, found_microglial2)

#Ensure the combined matrix has unique row names
astr_micr_matrix <- astr_micr_matrix[!duplicated(rownames(astr_micr_matrix)), ]

#save the data in a csv file
    write.csv(astr_micr_matrix, "astrocytic_and_microglial_genes.csv", row.names = FALSE)

          
##################################################################################################################################################
 
###Let's clean up the workspace before moving forward to batch correction...
      
#Remove unecessary files
    remove("astrocytic", 
           "astrocytic_filtered",
           "duplicated_rows", 
           "found_astrocytic_after_clean", 
           "found_astrocytic2", 
           "found_microglial_after_clean",
           "found_microglial2", 
           "microglial",
           "My_data",
           "My_data_clean",
           "My_data_dedup",
           "duplicated_gene_names",
           "duplicated_genes",
           "found_astrocytic",
           "found_microglial",
           "rows_out" ,
           "rows_to_keep"
           )
  
        
##################################################################################################################################################

##Script 4:
  ###This script performs batch corrections of an RNA-Seq data set, accounting for covariates such as age and sex
  ###This script will also reiterate the code to run multiple batch corrections (as many batch corrections as there are in the project_sheet)
  ###This will use as input the matrix "astr_micr_filtered"

    
#######################    
    
#Preparing the project sheet for ComBat
  #Remove the first 2 lines of the project sheet
  ###This code assumes that the first column in the project sheet has the names for the FastQ files
  project_sheet <-  project_sheet[-c(1,2),]
    project_sheet <- project_sheet[, -1] ###Turn this off if the project sheet does not have the FastQ names in the first column
      
  #Determine the number of batch columns starting from column 3
  start_col <- 3
    end_col <- ncol(project_sheet)
      num_batches <- end_col - start_col + 1

  #Generate batch names dynamically
  batch_names <- paste0("batch", seq_len(num_batches))
      
  #Combine the fixed column names with the dynamically generated batch names
  all_names <- c("samples", "condition", batch_names)
      
  #Assign the new names to the project_sheet data frame
  names(project_sheet) <-  all_names
  
  #Separate combined covariates into individual columns (In this case individuals were identified by "age_sex", ex.: "old_male")
  project_sheet <- project_sheet %>%
    separate(condition, into = c("age", "sex"), sep = "_")
    
  #Identify the batch columns starting from the 4th column
  batch_columns <- names(project_sheet)[4:ncol(project_sheet)] #Added a column to the "batch start column", because the "condition" column is now two columns ("age" and "sex")

  
#######################
        
#Preparing the count matrix for ComBat    
  #Convert astr_micr_filtered to a matrix, excluding the first two non-numeric columns
  astr_micr_filtered2 <- as.matrix(astr_micr_filtered[, -c(1, 2)])

  #Check and handle negative values
  if (any(astr_micr_filtered2 < 0)) {
    astr_micr_filtered2[astr_micr_filtered2 < 0] <- 0
  }

  #Match sample names
  temp <- names(astr_micr_filtered)[-c(1, 2)]
    matched_samples <- project_sheet[match(temp, project_sheet$samples), ]

          
#######################
      
#Detecting covariates from project sheet and creating a covariate reference
  #Check if the 'sex' column has only one unique value
  if (length(unique(matched_samples$sex)) == 1) {
    # Only one unique value, exclude 'sex' from the design matrix
    design <- model.matrix(~age, data = matched_samples)
  } else {  
    # Create the design matrix for covariates
    ### Include only the relevant covariates (excluding batch)
    design <- model.matrix(~age + sex, data = matched_samples) 
  }
  
  # Remove the intercept column from the design matrix if it's the only column
  if (ncol(design) == 1) {
    covar_mod <- NULL
  } else {
    covar_mod <- as.matrix(design[, -1])  # Exclude the intercept column
  }
 
     
#######################

#Iterate over the batch columns and run ComBat_seq for each one, updating the data matrix each time
  for (batch_col in batch_columns) {
    # Extract the batch information for the current batch column
    batch_info <- as.matrix(matched_samples[, batch_col])
    
    # Check dimensions
    stopifnot(ncol(astr_micr_filtered2) == nrow(batch_info))
    
    #Run ComBat for batch correction
    if (!is.null(covar_mod)) {
    astr_micr_filtered2 <- ComBat_seq(astr_micr_filtered2, batch=batch_info, group=NULL, covar_mod = covar_mod)
    } else {
    astr_micr_filtered2 <- ComBat_seq(astr_micr_filtered2, batch=batch_info, group=NULL)  
    }
  }
 
#Combine adjusted matrix with the first two columns of My_data  
adjusted_genes <- cbind(astr_micr_filtered[,c(1,2)],astr_micr_filtered2)  
 
#save the data in a csv file
    write.csv(adjusted_genes, "batch_corrected_matrix.csv", row.names = FALSE) 
      

##################################################################################################################################################

##Script 5:
  ###This script performs the same batch correction with the astrocytic and microglial genes list
  ###This uses "astr_micr_matrix" as input

#Convert astr_micr_matrix to a matrix, excluding the first two non-numeric columns
astr_micr_matrix2 <- as.matrix(astr_micr_matrix[, -c(1, 2)])

#Check and handle negative values
if (any(astr_micr_matrix2 < 0)) {
  astr_micr_matrix2[astr_micr_matrix2 < 0] <- 0
}

#Match sample names
temp_astr_micr <- names(astr_micr_matrix)[-c(1, 2)]
  matched_samples_astr_micr <- project_sheet[match(temp_astr_micr, project_sheet$samples), ] #Notice this is the already modified "project_sheet"

  
#######################

#Iterate over the batch columns and run ComBat_seq for each one, updating the data matrix each time
  for (batch_col in batch_columns) {
    # Extract the batch information for the current batch column
    batch_info <- as.matrix(matched_samples_astr_micr[, batch_col])
    
    #Check dimensions
    stopifnot(ncol(astr_micr_matrix2) == nrow(project_sheet))
    
    #Run ComBat for batch correction
    if (!is.null(covar_mod)) {
      astr_micr_adjusted <- ComBat_seq(astr_micr_matrix2, batch=batch_info, group=NULL, covar_mod = covar_mod)
    } else {
      astr_micr_adjusted <- ComBat_seq(astr_micr_matrix2, batch=batch_info, group=NULL)  
    }
  }
  
#Combine adjusted matrix with the first two columns
astr_micr_matrix_adjusted <- cbind(astr_micr_matrix[,c(1,2)],astr_micr_adjusted)

#save the data in a csv file
    write.csv(astr_micr_matrix_adjusted, "astr_micr_batch_corrected.csv", row.names = FALSE) 

    
##################################################################################################################################################

##Script 6:
  ###This script creates PCA plots from the corrected matrix
  ###This uses "adjusted_gene" as input (the final batch corrected matrix)

#Hardcoded values
label_points <- TRUE
WIDTH <- 20 # plot width
HEIGHT <- 10 # plot height
    
#Preparing "project_sheet" to assign group names
project_sheet <- as.data.frame(project_sheet, skip = 1, col_names = FALSE)

#Creating groups
sample_info <- as.data.frame(project_sheet)
    
  if (! "group" %in% colnames(sample_info)) {
      
    # If genotype and/or age, use it (or combination of them) for "group"
    genotype_column_present <- "genotype" %in% colnames(sample_info)
    condition_column_present <- "age" %in% colnames(sample_info)
    both_present <- genotype_column_present && condition_column_present
      
      if (both_present) {
        sample_info <- sample_info %>%
          dplyr::mutate(group = glue("{.$genotype}_{.$age}"))
      } else if (genotype_column_present) {
        sample_info <- sample_info %>%
          dplyr::mutate(group = .$genotype)
      } else if (condition_column_present) {
        sample_info <- sample_info %>%
          dplyr::mutate(group = .$age)
      } else {
        cat("Error. No group identified")
        q()
      }
  }
    
#Preparing "adjusted_genes" for PCA
counts <- as.data.frame(adjusted_genes) # Convert tibble to data frame
  counts_numeric <- counts [,sapply(counts,is.numeric)] # Select numeric columns
    
# Change samples from columns to rows
t_counts <- t(counts_numeric)

      
#######################
    
##Building PNG files with PCA data
    
#Open the drawing device.
png("PCA_corrected.png", width = 20, height = 10, units = "in", res = 300) #open the drawing device  
 
  par(mfrow = c(2, 1))
    pca_plot <- prcomp(t_counts, scale. = TRUE) %>% # calc principle components
      autoplot() +                        # plot
      ggtitle(aes("PCA")) + # Add title
      geom_point(aes(col = sample_info$group) , size = 4) + # color by group
      #geom_mark_ellipse(aes(col = sample_info$group) , expand = unit(0.5, "mm")) +
      theme(legend.title = element_blank() , # no legend title
            legend.text  = element_text (size = 12)) +
      scale_color_viridis(discrete = TRUE) # color-blind friendly
    
  #Add labels if desired
  if (label_points) {
      
    #Put labels above the points
    nudge <- position_nudge(y = 0.03)
      
    pca_plot <- pca_plot +
      geom_text(
        aes(label = sample_columns),
        position = nudge ,
        size = 4.5
        )
  }
    
#Print plot
print(pca_plot)
    
dev.off() #close PNG device
    
    
##################################################################################################################################################
    
##Script 7:
  ###This script creates PCA plots from the astrocytic and microglial genes list for comparison, it should overlap the corrected data's PCA    
  ###This will use "astr_micr_matrix_adjusted" as input
 
#Hardcoded values
label_points <- TRUE
WIDTH <- 20 # plot width
HEIGHT <- 10 # plot height
  
#Preparing "astr_micr_matrix_adjusted" for PCA
counts_astr_micr <- as.data.frame(astr_micr_matrix_adjusted) # Convert tibble to data frame
  counts_numeric_astr_micr <- counts_astr_micr [, sapply(counts, is.numeric)] # Select numeric columns
    
#Change samples from columns to rows
t_counts_astr_micr <- t(counts_numeric_astr_micr)

      
#######################
    
##Building PNG files with PCA data
    
#Open the drawing device.
png("PCA_astr_micr.png", width = 20, height = 10, units = "in", res = 300) #open the drawing device  
    
  par(mfrow = c(2, 1))
  pca_plot <- prcomp(t_counts, scale. = TRUE) %>% # calc principle components
    autoplot() +                        # plot
    ggtitle(aes("PCA")) + # Add title
    geom_point(aes(col = sample_info$group) , size = 4) + # color by group
    #geom_mark_ellipse(aes(col = sample_info$group) , expand = unit(0.5, "mm")) +
    theme(legend.title = element_blank() , # no legend title
          legend.text  = element_text (size = 12)) +
    scale_color_viridis(discrete = TRUE) # color-blind friendly
    
  #Add labels if desired
  if (label_points) {
      
    #Put labels above the points
    nudge <- position_nudge(y = 0.03)
      
    pca_plot <- pca_plot +
      geom_text(
        aes(label = sample_columns),
        position = nudge ,
        size = 4.5
      )
  }
    
#Print plot
print(pca_plot)
    
dev.off() #close PNG device
