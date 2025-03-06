#########################

####This script generates "Gene Ontology" (GO) pathways analysis from a DESeq2 output matrix
####GO Pathways analysis are more complex than KEGG analysis because the code runs 3 different analysis: "Biological Process" (BP), "Cellular Component" (CC) and "Molecular Function" (MF)
####This script also does the following:
  ###Generates a filtered GO result, based on pre-selected pathways (i.e., neuronal related pathways)
  ###Plots (Dotplot or Barplot) selected GO pathways analysis
  ###Plots Cnetplots

#######################

##Libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)#Use org.Hs.eg.db for human GO analysis
  library(DOSE)
  library(pathview)
  library(enrichplot)
})


########################

#Loading files:
DESeq2 <- as.data.frame(Old_vs_Young_10_DESeq2)
relevant_pathways <- as.data.frame(Neuron_Related_GO_Pathways)


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
###This can be used to produce GO analysis for "All DEG", "Upregulated DEG" or "Downregulated DEG"
DESeq2_ordered <- DESeq2_clean[order(DESeq2_clean$pAdj, abs(DESeq2_clean$log2FoldChange)), ]
  DEG <- DESeq2_ordered[DESeq2_ordered$pAdj < 0.5 & abs(DESeq2_ordered$log2FoldChange) > 0.5, ]     # <------------- CHOOSE CUT OFF VALUE HERE
    #DEG <- DESeq2_ordered[DESeq2_ordered$pAdj < 0.5 & DESeq2_ordered$log2FoldChange > 0.5, ]       #Upregulated <------------- CHOOSE CUT OFF VALUE HERE
    #DEG <- DESeq2_ordered[DESeq2_ordered$pAdj < 0.5 & DESeq2_ordered$log2FoldChange < -0.5, ]      #Downregulated <------------- CHOOSE CUT OFF VALUE HERE

#Extract columns of interest
DESeq2_GO <- DEG[,c("gene_symbol","log2FoldChange","pAdj")]

#Extract invalid gene symbols from the table, a.k.a. "gene_symbols" that can't be mapped by GO
valid_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")             #Extract valid gene symbols from org.Mm.eg.db
  valid_gene_symbols <- DESeq2_GO$gene_symbol %in% valid_symbols  #Check which of your "gene_symbols" are valid
    invalid_symbols <- DESeq2_GO$gene_symbol[!valid_gene_symbols] #Display invalid gene symbols
      DESeq2_GO <- DESeq2_GO[valid_gene_symbols, ]              #Filter out invalid gene symbols from DESeq2_GO


##################################################################################################################################################

#GO pathways enrichment analysis using "DESeq2_GO" as input:

      
#######################          

#Convert gene symbols to Entrez IDs
entrez_degs <- bitr(DESeq2_GO$gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

#Separate GO analysis into Biological Process, Cellular Component and Molecular Function
  #Perform GO enrichment analysis for Biological Process (BP)
  go_bp <- enrichGO(gene          = entrez_degs$ENTREZID,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = "ENTREZID",
                    ont           = "BP",          # Ontology: BP
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,                               # <------------- CHOOSE CUT OFF VALUE HERE
                    qvalueCutoff  = 0.2,                                # <------------- CHOOSE CUT OFF VALUE HERE
                    readable      = TRUE)
  
    #Convert Biological Process to data frame and save to CSV
    GO_BP <- as.data.frame(go_bp)
      write.csv(GO_BP, file = "GO_BP_Results.csv", row.names = FALSE)
  
      
  #Perform GO enrichment analysis for Cellular Component (CC)
  go_cc <- enrichGO(gene          = entrez_degs$ENTREZID,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = "ENTREZID",
                    ont           = "CC",          # Ontology: CC  
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,                               # <------------- CHOOSE CUT OFF VALUE HERE
                    qvalueCutoff  = 0.2,                                # <------------- CHOOSE CUT OFF VALUE HERE
                    readable      = TRUE)
  
  
    #Convert Cellular Component to data frame and save to CSV
    GO_CC <- as.data.frame(go_cc)
      write.csv(GO_CC, file = "GO_CC_Results.csv", row.names = FALSE)
  
      
  #Perform GO enrichment analysis for Molecular Function (MF)
  go_mf <- enrichGO(gene          = entrez_degs$ENTREZID,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = "ENTREZID",
                    ont           = "MF",          # Ontology: MF
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,                               # <------------- CHOOSE CUT OFF VALUE HERE
                    qvalueCutoff  = 0.2,                                # <------------- CHOOSE CUT OFF VALUE HERE
                    readable      = TRUE)
  
    #Convert Molecular Function to data frame and save to CSV
    GO_MF <- as.data.frame(go_mf)
      write.csv(GO_MF, file = "GO_MF_Results.csv", row.names = FALSE)
 
           
#######################          
  
#Combine the results for BP, CC, and MF into a single data frame
go_combined <- rbind(
                    data.frame(go_bp, Ontology = "Biological Process"),
                    data.frame(go_cc, Ontology = "Cellular Component"),
                    data.frame(go_mf, Ontology = "Molecular Function")
                    )
   
#Save the combined data frame to to CSV
  write.csv(go_combined, file = "GO_Combined.csv", row.names = FALSE)

 
#######################  

#Plot the top 10 GO terms for each category
  
  #Filter the top 10 GO terms for each ontology category
  go_top10 <- go_combined %>%
    group_by(Ontology) %>%
      slice_max(order_by = Count, n = 5)
  
  #Filter out certain IDs (eg.: Pathways that are not neuronal related)
  go_ids_to_remove <- c("GO:0007015", "GO:0032970", "GO:0042692", "GO:0051235", "	GO:0030863", "GO:0030016","GO:0043292", "GO:0030864", "GO:0030863") #Use this to customize your "Top 10" list
    go_top10 <- go_top10[!go_top10$ID %in% go_ids_to_remove, ] 
  
  #Create a bar plot for the top 10 GO terms in each category
  ggplot(go_top10, aes(x = reorder(Description, -Count), y = Count, fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +  # Create separate plots for BP, CC, MF
    labs(title = "Top GO Enrichment Terms", x = "GO Terms", y = "Gene Count") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 20, face = "bold"),
      legend.position = "none",  # Remove legend if you want
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 18, face = "bold"),
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
    ) +
    scale_fill_manual(values = c("Biological Process" = "skyblue", 
                                 "Cellular Component" = "lightgreen", 
                                 "Molecular Function" = "lightcoral"))

#Save "Top 10"  enrichment barplot    
    ggsave("GO_enrichment_barplot.png", width = 14, height = 8, dpi = 600)
    
    
##################################################################################################################################################

#Select GO pathways of interest from the analysis using "go_combined" as input:


#######################

#Filter the GO results to only include the custom pathways 
selected_combined <- go_combined[ go_combined$ID %in% relevant_pathways$ID, ]

#Create output table  
  #Convert results to a data frame
  selected_combined <- as.data.frame(selected_combined)
    
  #Save selected GO results to CSV
      write.csv(selected_combined, file = "selected_GO_Results.csv", row.names = FALSE)
      
  
##################################################################################################################################################

#Generating plots: 
###Generic plots  
  
#######################  
      
#Plot the top 10 GO terms for each category, from selected pathways
###Same plot as above, but now with selected pathways only 
      
  #Filter the top 10 GO terms for each ontology category
  go_top10 <- selected_combined %>%
    group_by(Ontology) %>%
      slice_max(order_by = Count, n = 5)
      
  #Filter out certain IDs (eg.: Pathways that are not neuronal related)
  go_ids_to_remove <- c("GO:0007015", "GO:0032970", "GO:0042692", "GO:0051235", "	GO:0030863", "GO:0030016","GO:0043292", "GO:0030864", "GO:0030863") #Use this to customize your "Top 10" list
    go_top10 <- go_top10[!go_top10$ID %in% go_ids_to_remove, ] 
      
  #Create a bar plot for the top 10 GO terms in each category
  ggplot(go_top10, aes(x = reorder(Description, -Count), y = Count, fill = Ontology)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +  #Create separate plots for BP, CC, MF
    labs(title = "Top GO Enrichment Terms", x = "GO Terms", y = "Gene Count") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 20, face = "bold"),
      legend.position = "none",  #Remove legend if you want
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 18, face = "bold"),
      panel.grid.major = element_blank(),  #Remove major grid lines
      panel.grid.minor = element_blank()   #Remove minor grid lines
    ) +
    scale_fill_manual(values = c("Biological Process" = "skyblue", 
                                 "Cellular Component" = "lightgreen",  
                                 "Molecular Function" = "lightcoral"))
      
#Save "Top GO Terms"  bar plot in PNG file with high resolution     
    ggsave("GO_enrichment_barplot.png", width = 14, height = 8, dpi = 600)

    
#######################
    
#Generic bar plot selected GO pathways:
###This is the most simple approach
    
barplot(go_bp, showCategory = 10, title = "Top 10 GO Terms")  # <------------- CHOOSE BETWEEN "bp", "cc" or "mf"HERE

#Save plot in PNG file with high resolution  
    ggsave("GO_barplot.png", width = 9, height = 8, dpi = 600)      
    
    
#######################    
    
#Generic "Gene concept network" plot for GO pathways:    
###This is the most simple approach
    
cnetplot(go_bp,                                               # <------------- CHOOSE BETWEEN "bp", "cc" or "mf"HERE
          color.params = list(foldChange = DESeq2_GO$log2FoldChange),
          layout = "kk",          # Layout options: "circle", "kk" (Kamada-Kawai), "fr" (Fruchterman-Reingold)
          node_label = "all",          # Display labels on all nodes
          cex.params = list(gene_label = 0.9,        # Adjust gene label size
                            category_label = 1,    # Adjust category (GO term) label size
                            category_node = 1.5,          # Adjust node size for GO terms
                            gene_node = 1),              # Adjust node size for genes
          circular = FALSE              # Optional: Use a circular layout for cell talk plot
)

#Save plot in PNG file with high resolution  
    ggsave("GO_cnetplot.png", width = 9, height = 8, dpi = 600)      

        
##################################################################################################################################################
    
#Generating plots: 
###Custom plots  
    
####################### 

#Custom dot plot for selected GO pathways:
###This code can be used in a saved data frame that has been reloaded to the environment
###To replace "pvalue" with "Log10 P-value" add that information prior to plotting and change the code from "pvalue" to "Log10_P" and "-Log10 P-value"

#Start by giving your data frame a generic name "go_df"
go_df <- selected_GO_Results                             # <------------- IF NEEDED, CHOOSE THE FILE NAME HERE (DEFAULT: Combined GO after selection)
    
  #Mutate the data frame to be used by ggplot
  go_df <- go_df %>%
      mutate(
        GeneRatio = as.numeric(gsub("/.*", "", GeneRatio)) / as.numeric(gsub(".*?/", "", GeneRatio)),
        BgRatio = as.numeric(gsub("/.*", "", BgRatio)) / as.numeric(gsub(".*?/", "", BgRatio)),
        Count = as.numeric(Count),
        pvalue = as.numeric(pvalue)                 # <------------- CUSTOMIZE "pvalue" HERE, TO CHANGE TO "LOG10_P"
      )
    
    #Split the genes and convert into a list
    go_df$gene_list <- strsplit(go_df$geneID, "/")
    
    go_df <- go_df %>%
      arrange(desc(Description))  # Orders alphabetically from A to Z

#Set fixed size limits for both plots (adjust according to your data)
max_size <- 20  #Maximum size for the largest count
min_size <- 1   #Minimum size for the smallest count
  
#Create a dotplot
ggplot(go_df, aes(x = reorder(Description, GeneRatio),
                            y = GeneRatio, 
                            size = Count, 
                            color = pvalue)) +      # <------------- CUSTOMIZE "pvalue" HERE, TO CHANGE TO "LOG10_P"
    geom_point() +
    scale_color_gradient2(
      
      low = "#40a6cf",      #Color for low p-values
      mid = "#40a6cf",      #Color for mid p-values
      high = "#ed3131",     #Color for high p-values
      #midpoint = 0.0005,    #Specify the midpoint value for the gradient
      limits = c(0.0000005, 0.0005)      #Define the limits for the color scale, based on the "pvalue" output
    
    ) +
    
    scale_size_continuous(range = c(min_size, max_size), limits = c(1, 20)) +   #Set fixed size limits. Adjust based on the actual "Count" output
    
    scale_y_continuous(expand = expansion(add =c(0.008, 0.008))) +                #Control the range and spacing of the y-axis. This fixes the issues with dots missing from the plot
    
    theme_minimal() +               #Clean theme
    coord_flip() +                  #Flip axes for horizontal labels
    labs(
      
      title = "DAT;NuTRAP GO",                       # <------------- CHOOSE GRAPHIC TITLE HERE
      x = "",
      y = "Gene Ratio",
      color = "P-value",                             # <-------------"P-value" can be replaced by "-Log10 P-Value", but this needs to be calculated and added to the data frame prior to plotting
      size = "Gene Count"
    
    ) +
  
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
      axis.text.y = element_text(color = "#141313", size = 16),         #Adjust size of y-axis labels
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
    ggsave("GO_dotplot.png", width = 9, height = 8.5, dpi = 600)


#######################
    
#Gene-concept-network plot for selected GO pathways:    
###This code uses the same "go_df" input as created above
    
#Prepare data frame for cnet plot
go_df_processed <- go_df %>%   #Split and unnest genes
  mutate(genes = strsplit(geneID, "/")) %>%
  unnest(genes)
    
  #Convert gene symbols to Entrez IDs
  gene_ids <- mapIds(org.Mm.eg.db,                          #Change to appropriate organism database
                     keys = unique(go_df_processed$genes),
                     column = "ENTREZID",
                     keytype = "SYMBOL")
    
    #Create gene list for cnet plot
    go_list <- split(go_df_processed$genes, go_df_processed$Description)
    

#Create cnetplot
p <- cnetplot(go_list, 
              showCategory = 10, 
              node_label = "all",
              color.params = list(foldChange = NULL),
              ) +
  
              theme(text = element_text(size = 12),
                    plot.title = element_text(size = 14, face = "bold"),
                    legend.text = element_text(size = 16),                            #Increase legend text size
                    legend.title = element_text(size = 18, face = "bold"),            #Increase legend title size
                    legend.key.size = unit(1, "cm"),                                #Increase size of legend keys (points)
                    plot.margin = margin(5, 5, 5, 5),                                 #Add margins (top, right, bottom, left)
              ) +
  
              labs(title = "")
    
#Print the plot
    print(p)
    
#Save plot in PNG file with high resolution  
    ggsave("GO_cnetplot.png", width = 9, height = 8, dpi = 600) 
    
    

  