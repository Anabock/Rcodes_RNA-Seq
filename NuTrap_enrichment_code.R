#########################

#This code screens the enrichment DESeq2 output for each one of the gene lists and outputs tables that when copied to GraphPad Prism generate enrichment visualization graphs


#######################

#Loading files:
genes_lists <- as.data.frame(Table_S2_Enrichment_genes_lists)
My_data <- as.data.frame(enrichment_DESeq2_output)


#######################

#1.Neuronal
#search for a list of gene names in a dataset
found <- My_data[ My_data$gene_symbol %in% genes_lists$Neuronal, ]
#Extract columns of interest
found <- found[,c("gene_symbol","log2FoldChange","pvalue","pAdj")]
write.csv(found, "Neuronal_enrichment.csv", row.names = FALSE)


#2.Dopaminergic
#search for a list of gene names in a dataset
found <- My_data[ My_data$gene_symbol %in% genes_lists$Dopaminergic, ]
#Extract columns of interest
found <- found[,c("gene_symbol","log2FoldChange","pvalue","pAdj")]
write.csv(found, "Dopaminergic_enrichment.csv", row.names = FALSE)


#3.GABAergic
#search for a list of gene names in a dataset
found <- My_data[ My_data$gene_symbol %in% genes_lists$GABAergic, ]
#Extract columns of interest
found <- found[,c("gene_symbol","log2FoldChange","pvalue","pAdj")]
write.csv(found, "GABAergic_enrichment.csv", row.names = FALSE)


#4.Astrocytic
#search for a list of gene names in a dataset
found <- My_data[ My_data$gene_symbol %in% genes_lists$Astrocytic, ]
#Extract columns of interest
found <- found[,c("gene_symbol","log2FoldChange","pvalue","pAdj")]
write.csv(found, "Astrocyte_enrichment.csv", row.names = FALSE)


#5.Microglial
#search for a list of gene names in a dataset
found <- My_data[ My_data$gene_symbol %in% genes_lists$Microglial, ]
#Extract columns of interest
found <- found[,c("gene_symbol","log2FoldChange","pvalue","pAdj")]
write.csv(found, "Microglial_enrichment.csv", row.names = FALSE)


#6.Oligodendrocytes
#search for a list of gene names in a dataset
found <- My_data[ My_data$gene_symbol %in% genes_lists$Oligodendrocytic, ]
#Extract columns of interest
found <- found[,c("gene_symbol","log2FoldChange","pvalue","pAdj")]
write.csv(found, "Oligodendrocytes_enrichment.csv", row.names = FALSE)


#7.Endothelial
#search for a list of gene names in a dataset
found <- My_data[ My_data$gene_symbol %in% genes_lists$Endothelial, ]
#Extract columns of interest
found <- found[,c("gene_symbol","log2FoldChange","pvalue","pAdj")]
write.csv(found, "Endothelial_enrichment.csv", row.names = FALSE)

