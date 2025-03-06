#########################

####This script generates "Kyoto Encyclopedia of Genes and Genomes" (KEGG) analysis from a DESeq2 output matrix
####This script also does the following:
  ###Generates a filtered KEGG result, based on pre-selected pathways (i.e., neuronal related pathways)
  ###Plots (Dotplot) pre-selected KEGG pathways analysis

#######################

##Libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)#Use org.Hs.eg.db for human KEGG analysis
  library(DOSE)
  library(pathview)
  library(enrichplot)
})


########################

#Loading files:
DESeq2 <- as.data.frame(Old_vs_Young_10_DESeq2)
relevant_pathways <- as.data.frame(Neuron_Related_KEGG_Pathways)


##################################################################################################################################################

##Getting "DESeq2" data frame ready for the analysis:
  ###Eliminates any potential issues with the table, such as low counts, "0", "NA"
  ###Order the genes based on "Adjp" and then based on "log2FoldChange"
  ###Creates DEG list based on chosen cut off
  ###Creates DEG list with only 3 variables: "gene_symbol", "Adjp" and "log2FoldChange"

#######################

#Remove rows with sums lower than determined cut off 
###(Default cut off = 1)
sample_columns <- colnames(DESeq2)[12:ncol(DESeq2)]  
  num_samples <- ncol(DESeq2)-11
    rows_to_keep <- rowSums(DESeq2[, sample_columns] >= 1) == num_samples   # <------------- CHOOSE CUT OFF VALUE HERE
      DESeq2_clean <- DESeq2[rows_to_keep, ]

#Check for "NA" values in the gene symbols and correct it
if (any(is.na(DESeq2_clean$gene_symbol))) {
  print("Warning: NA values found in gene symbols")
    DESeq2_clean <- DESeq2_clean %>% drop_na(gene_symbol) #Remove rows with NA in the gene_symbol column
}

#Order results based on "pAdj" and "log2foldchange"
###This can be used to produce KEGG analysis for "All DEG", "Upregulated DEG" or "Downregulated DEG"
DESeq2_ordered <- DESeq2_clean[order(DESeq2_clean$pAdj, abs(DESeq2_clean$log2FoldChange)), ]
  DEG <- DESeq2_ordered[DESeq2_ordered$pAdj < 0.5 & abs(DESeq2_ordered$log2FoldChange) > 0.5, ]     # <------------- CHOOSE CUT OFF VALUE HERE
    #DEG <- DESeq2_ordered[DESeq2_ordered$pAdj < 0.5 & DESeq2_ordered$log2FoldChange > 0.5, ]       #Upregulated <------------- CHOOSE CUT OFF VALUE HERE
    #DEG <- DESeq2_ordered[DESeq2_ordered$pAdj < 0.5 & DESeq2_ordered$log2FoldChange < -0.5, ]      #Downregulated <------------- CHOOSE CUT OFF VALUE HERE

#Extract columns of interest
DESeq2_KEGG <- DEG[,c("gene_symbol","log2FoldChange","pAdj")]

#Extract invalid gene symbols from the table, a.k.a. "gene_symbols" that can't be mapped by KEGG
valid_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")             #Extract valid gene symbols from org.Mm.eg.db
  valid_gene_symbols <- DESeq2_KEGG$gene_symbol %in% valid_symbols  #Check which of your "gene_symbols" are valid
    invalid_symbols <- DESeq2_KEGG$gene_symbol[!valid_gene_symbols] #Display invalid gene symbols
      DESeq2_KEGG <- DESeq2_KEGG[valid_gene_symbols, ]              #Filter out invalid gene symbols from DESeq2_KEGG


##################################################################################################################################################

#KEGG pathways enrichment analysis using "DESeq2_KEGG" as input:

      
#######################          

#Convert gene symbols to Entrez IDs
entrez_degs <- bitr(DESeq2_KEGG$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

#Perform KEGG pathways enrichment analysis
kegg_results <- enrichKEGG(gene = entrez_degs$ENTREZID, organism = 'mmu', pvalueCutoff = 0.5)   # <------------- CHOOSE CUT OFF VALUE HERE / Use 'mmu' for mouse 
  if (is.null(kegg_results))  #Check the results
    {
      
      cat("No KEGG results found. Please check the input gene IDs and organism code.\n")
      
      } else 
    {
      
    #Create output table 
      #Customize labels to simplify "Description"
      kegg_results@result$Description <- gsub("-", "", gsub("Mus musculus \\(house mouse\\)", "", kegg_results@result$Description))  # Remove "- Mus musculus (house mouse)" from the description category
  
      #Convert results to data frame for customization
      kegg_df <- as.data.frame(kegg_results)
    
    #Save KEGG results to CSV
        write.csv(kegg_df, file = "KEGG_results.csv", row.names = FALSE)
  
  }


##################################################################################################################################################

#Select KEGG pathways of interest from the analysis using "kegg_results" as input:


#######################

#Filter the KEGG results to only include the custom pathways 
selected_kegg_results <- kegg_results
  selected_kegg_results@result <- selected_kegg_results@result[selected_kegg_results@result$ID %in% relevant_pathways$Pathway_ID, ] #This will only work if the terms "- Mus musculus house mouse" have been removed from the description in the KEGG results output
  
#Visualize the results if available
if (nrow(selected_kegg_results@result) > 0) 
  {
  
  #Create output table  
    #Convert results to a data frame
    kegg_df_selected <- as.data.frame(selected_kegg_results)
    
  #Save selected KEGG results to CSV
      write.csv(kegg_df_selected, file = "selected_KEGG_Results.csv", row.names = FALSE)
      
  } else {
    
    cat("No enriched pathways found in the selected results.\n")
} 

    
##################################################################################################################################################

#Generating plots: 
  
  
#######################

#Dotplot for selected KEGG pathways:

#Mutate the data frame to be used by ggplot
###This code can be used in a saved data frame that has been reloaded to the environment
###To replace "pvalue" with "Log10 P-value" add that information prior to plotting and change the code
kegg_df_selected <- as.data.frame(kegg_df_selected) #Replace "kegg_df_selected" to use it with a previously saved data frame
  kegg_df_selected <- kegg_df_selected %>%
    mutate(
      GeneRatio = as.numeric(gsub("/.*", "", GeneRatio)) / as.numeric(gsub(".*?/", "", GeneRatio)),
      BgRatio = as.numeric(gsub("/.*", "", BgRatio)) / as.numeric(gsub(".*?/", "", BgRatio)),
      Count = as.numeric(Count),
      pvalue = as.numeric(pvalue)
    )
  
#Set fixed size limits for both plots (adjust according to your data)
max_size <- 20  #Maximum size for the largest count
min_size <- 1   #Minimum size for the smallest count
  
#Create a dotplot
ggplot(kegg_df_selected, aes(x = reorder(Description, GeneRatio),
                            y = GeneRatio, 
                            size = Count, 
                            color = pvalue)) +
    geom_point() +
    scale_color_gradient2(
      low = "#40a6cf",      #Color for low p-values
      mid = "#40a6cf",      #Color for mid p-values
      high = "#ed3131",     #Color for high p-values
      midpoint = 0.0005,    #Specify the midpoint value for the gradient
      limits = c(0.000005, 0.05)      #Define the limits for the color scale, based on the "pvalue" output
    ) +
    
    scale_size_continuous(range = c(min_size, max_size), limits = c(1, 20)) +   #Set fixed size limits. Adjust based on the actual "Count" output
    
    scale_y_continuous(expand = expansion(add =c(0.01, 0.01))) +                #Control the range and spacing of the y-axis. This fixes the issues with dots missing from the plot
    
    theme_minimal() +               #Clean theme
    coord_flip() +                  #Flip axes for horizontal labels
    labs(
      title = "DAT;NuTRAP KEGG",                                                # <------------- CHOOSE GRAPHIC TITLE HERE
      x = "",
      y = "Gene Ratio",
      color = "P-value",    #"P-value" can be replaced by "-Log10 P-Value", but this needs to be calculated and added to the data frame prior to plotting
      size = "Gene Count"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.text.y = element_text(color = "#141313", size = 18),         #Adjust size of y-axis labels
      axis.text.x = element_text(hjust = 2, size = 14),                 #Adjust size of x-axis labels
      axis.title.x = element_text(size = 16),                           #Smaller x-axis title
      axis.title.y = element_blank(),                                   #No y-axis title
      panel.grid.major = element_line(color = "grey90", size = 0.5),    #Lighter major grid lines
      panel.grid.minor = element_blank(),                               #Remove minor grid lines
      panel.grid.major.x = element_blank(),                             #Remove vertical grid lines
      panel.grid.major.y = element_line(size = 0.5),                    #Keep light horizontal grid lines
      legend.text = element_text(size = 16),                            #Increase legend text size
      legend.title = element_text(size = 18, face = "bold"),            #Increase legend title size
      legend.key.size = unit(1.5, "cm"),                                #Increase size of legend keys (points)
      plot.margin = margin(5, 5, 5, 5)                                  #Add margins (top, right, bottom, left)
    )
  
#Save plot in PNG file with high resolution  
    ggsave("KEGG_dotplot.png", width = 9, height = 8, dpi = 600)
  