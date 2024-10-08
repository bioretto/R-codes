setwd("C:/Users/mrinm/OneDrive/Documents/RfoldWD")

# Load required packages
library(GEOquery)

# Download GEO dataset
gse <- getGEO("GSE199633", GSEMatrix = TRUE)
gse <- gse[[1]]

# Extract expression data and sample information
exprs_data <- exprs(gse)
sample_info <- pData(gse)

# Check alignment of column names and row names
all(colnames(exprs_data) %in% rownames(sample_info))
all(rownames(sample_info) %in% colnames(exprs_data))

common_samples <- intersect(colnames(exprs_data), rownames(sample_info))

exprs_data <- exprs_data[, common_samples]
sample_info <- sample_info[common_samples, ]

# Ensure count data contains only integers and verify
exprs_data <- round(exprs_data)
summary(exprs_data)

# Load required packages for data manipulation
library(dplyr)

# Create receptor_status column with updated logic
sample_info <- sample_info %>%
  mutate(receptor_status = case_when(
    characteristics_ch1.4 %in% c("er status: Negative") & 
      characteristics_ch1.5 %in% c("pr status: Negative") &
      characteristics_ch1.6 %in% c("her2 status: Negative") ~ "TNBC",
    
    characteristics_ch1.4 %in% c("er status: Positive") & 
      characteristics_ch1.5 %in% c("pr status: Positive") &
      characteristics_ch1.6 %in% c("her2 status: Positive") ~ "ER_PR_HER2_Positive",
    
    characteristics_ch1.4 %in% c("er status: Positive") ~ "ER_Positive",
    characteristics_ch1.5 %in% c("pr status: Positive") ~ "PR_Positive",
    characteristics_ch1.6 %in% c("her2 status: Positive") ~ "HER2_Positive",
    
    characteristics_ch1.4 %in% c("er status: Negative") ~ "ER_Negative",
    characteristics_ch1.5 %in% c("pr status: Negative") ~ "PR_Negative",
    characteristics_ch1.6 %in% c("her2 status: Negative") ~ "HER2_Negative",
    
    TRUE ~ "Other"
  ))

# Convert receptor_status to factor
sample_info$receptor_status <- factor(sample_info$receptor_status)

# Relevel receptor_status to set TNBC as the reference
sample_info$receptor_status <- relevel(sample_info$receptor_status, ref = "TNBC")

# Convert sample_info to data frame
sample_info <- as.data.frame(sample_info)

# Set GSM identifiers as row names
sample_info$GSM <- paste0("GSM", 5979001:5979637)
rownames(sample_info) <- sample_info$GSM
sample_info$GSM <- NULL

# Get common samples
common_samples <- intersect(colnames(exprs_data), rownames(sample_info))

# Align exprs_data and sample_info by their sample identifiers
exprs_data <- exprs_data[, common_samples]
sample_info <- sample_info[common_samples, ]

# Ensure count data contains only integers
exprs_data <- round(exprs_data)

# Load required packages for differential expression analysis
library(DESeq2)

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = exprs_data,
                              colData = sample_info,
                              design = ~ receptor_status)

# Filter out low-count genes
dds <- dds[rowSums(counts(dds)) > 20, ]

vsd <- vst(dds, blind = FALSE)
#rld <- rlog(dds, blind = FALSE)
plotPCA(vsd, intgroup = "receptor_status")

# Perform differential expression analysis
dds <- DESeq(dds)
plotDispEsts(dds)

# Check levels of receptor_status
levels(sample_info$receptor_status)

# Check results names
resultsNames(dds)

# Define valid contrasts based on results names
valid_contrasts <- resultsNames(dds)[-1]  # Exclude "Intercept"

# Extract valid pairs from contrast names
valid_pairs <- lapply(valid_contrasts, function(x) {
  strsplit(sub("receptor_status_", "", x), "_vs_")[[1]]
})

# Load necessary libraries
library(dplyr)
library(DESeq2)
library(tibble)
library(GEOquery)

# Download the GPL15048 annotation file
gpl <- getGEO("GPL15048", AnnotGPL = TRUE)
annotation <- Table(gpl)

# Ensure the annotation data has the necessary columns
annotation <- annotation %>% select(ID, GeneSymbol)

# Define the path to save the results
save_path <- "C:/Users/mrinm/OneDrive/Documents/RfoldWD/results"

# Ensure the directory exists
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

# Define comparisons excluding the reference level (TNBC)
tnbc_comparisons <- setdiff(levels(sample_info$receptor_status), "TNBC")

# Initialize a list to store results
results_list <- list()

# Perform differential expression analysis for each comparison
for (group in tnbc_comparisons) {
  res <- results(dds, contrast = c("receptor_status", group, "TNBC"), independentFiltering = TRUE)
  comparison_name <- paste(group, "vs TNBC")
  
  # Convert results to a data frame and add row names as a column
  res_df <- as.data.frame(res) %>% rownames_to_column(var = "ID")
  
  # Merge the annotation data with the differential expression results
  merged_data <- merge(res_df, annotation, by = "ID", all.x = TRUE)
  
  # Filter significant genes
  sig_genes <- merged_data %>% filter(padj < 0.1)
  
  # Check if there are any significant genes
  if (nrow(sig_genes) > 0) {
    # Add a new column for regulation status
    sig_genes <- sig_genes %>% mutate(regulation = ifelse(log2FoldChange > 0, "upregulated", "downregulated"))
    
    # Save the results in the list
    results_list[[comparison_name]] <- sig_genes
    
    # Save the results to a CSV file
    write.csv(sig_genes, file = paste0(save_path, "/", comparison_name, "_significant_genes.csv"), row.names = TRUE)
    
    # Save the results to a text file
    write.table(sig_genes, file = paste0(save_path, "/", comparison_name, "_significant_genes.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
  }
}

# Check the results for each comparison
for (comparison_name in names(results_list)) {
  cat("Summary for", comparison_name, ":\n")
  print(summary(results_list[[comparison_name]]))
  cat("\n")
}

# Install and load required packages
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}
if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  install.packages("htmlwidgets")
}

library(plotly)
library(htmlwidgets)

# Function to create interactive volcano plot
create_volcano_plot <- function(res, comparison) {
  res_df <- as.data.frame(res)
  res_df$log2FoldChange <- res_df$log2FoldChange
  res_df$negLog10Pval <- -log10(res_df$pvalue)
  
  p <- plot_ly(data = res_df, x = ~log2FoldChange, y = ~negLog10Pval, text = ~rownames(res_df),
               type = 'scatter', mode = 'markers', marker = list(size = 5)) %>%
    layout(title = paste("Volcano Plot:", comparison),
           xaxis = list(title = "Log2 Fold Change"),
           yaxis = list(title = "-Log10 P-value"))
  
  return(p)
}

# Create and save volcano plots for significant comparisons
for (comparison_name in names(results_list)) {
  res <- results_list[[comparison_name]]
  if (any(res$padj < 0.1, na.rm = TRUE)) {
    p <- create_volcano_plot(res, comparison_name)
    saveWidget(p, file = paste0("C:/Users/mrinm/OneDrive/Documents/RfoldWD/plots/volcano_plot_", comparison_name, ".html"))
  }
}

# Install and load required packages
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}
if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  install.packages("htmlwidgets")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
# Install and load required packages
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}

# Function to create interactive volcano plot with gene symbols and colors
create_volcano_plot <- function(res, comparison) {
  res_df <- as.data.frame(res)
  res_df$negLog10Pval <- -log10(res_df$pvalue)
  res_df$gene <- res_df$GeneSymbol  # Use GeneSymbol for hover text
  
  # Add regulation status
  res_df$regulation <- ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated")
  
  p <- plot_ly(data = res_df, x = ~log2FoldChange, y = ~negLog10Pval, text = ~gene,
               type = 'scatter', mode = 'markers',
               marker = list(size = 5, color = ~regulation, colors = c("Upregulated" = "red", "Downregulated" = "blue"))) %>%
    layout(title = paste("Volcano Plot:", comparison),
           xaxis = list(title = "Log2 Fold Change"),
           yaxis = list(title = "-Log10 P-value"))
  
  return(p)
}

# Create and save volcano plots for significant comparisons
for (comparison_name in names(results_list)) {
  res <- results_list[[comparison_name]]
  if (any(res$padj < 0.1, na.rm = TRUE)) {
    p <- create_volcano_plot(res, comparison_name)
    saveWidget(p, file = paste0("C:/Users/mrinm/OneDrive/Documents/RfoldWD/plots/volcano_plot_", comparison_name, ".html"))
  }
}

create_ma_plot <- function(res, comparison) {
  res_df <- as.data.frame(res)
  res_df$mean <- res_df$baseMean
  res_df$gene <- res_df$GeneSymbol  # Use GeneSymbol for hover text
  
  # Add regulation status
  res_df$regulation <- ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated")
  
  p <- plot_ly(data = res_df, x = ~mean, y = ~log2FoldChange, text = ~gene,
               type = 'scatter', mode = 'markers',
               marker = list(size = 5, color = ~regulation, colors = c("Upregulated" = "red", "Downregulated" = "blue"))) %>%
    layout(title = paste("MA Plot:", comparison),
           xaxis = list(title = "Mean Expression"),
           yaxis = list(title = "Log2 Fold Change"))
  
  return(p)
}

# Create and save MA plots for significant comparisons
for (comparison_name in names(results_list)) {
  res <- results_list[[comparison_name]]
  if (any(res$padj < 0.1, na.rm = TRUE)) {
    p <- create_ma_plot(res, comparison_name)
    saveWidget(p, file = paste0("C:/Users/mrinm/OneDrive/Documents/RfoldWD/plots/ma_plot_", comparison_name, ".html"))
  }
}




# Define the directory containing the results files
results_dir <- "C:/Users/mrinm/OneDrive/Documents/RfoldWD/results"

# Define the substring of interest
substring_of_interest <- "galnac"

# List all CSV files in the results directory
results_files <- list.files(results_dir, pattern = "*.csv", full.names = TRUE)

# Initialize a list to store search results
search_results <- list()

# Loop through each file in the results directory
for (file in results_files) {
  # Read the data from the file
  res <- read.csv(file, stringsAsFactors = FALSE)
  
  # Find all genes that contain the substring in their GeneSymbol (case-insensitive)
  matching_genes <- grep(substring_of_interest, res$GeneSymbol, value = TRUE, ignore.case = TRUE)
  
  # Check if there are any matching genes
  if (length(matching_genes) > 0) {
    # Extract rows corresponding to the matching genes
    gene_data <- res[res$GeneSymbol %in% matching_genes, ]
    
    # Store the results in the search_results list with the file name as the key
    search_results[[basename(file)]] <- gene_data
  }
}

# Display the search results
if (length(search_results) > 0) {
  for (file_name in names(search_results)) {
    cat("Results for genes containing", substring_of_interest, "in file", file_name, ":\n")
    print(search_results[[file_name]])
    cat("\n")
  }
} else {
  cat("No genes containing", substring_of_interest, "found in any comparison file.\n")
}
