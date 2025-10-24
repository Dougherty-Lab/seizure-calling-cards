# Gene overlap testing using permutations

#install and load needed packages
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(readr)

# read CC gene list
df <- read.csv("genelist.csv")

# read disease gene lists
mono_ep_pg_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","Monogenic_Epilepsy_PG",col_names = FALSE)
all_ep_pg_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","All_Epilepsy_PG",col_names = FALSE)
all_ep_zhang_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","All_Epilepsy_Zhang",col_names = FALSE)
mono_ep_zhang_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","Monogenic_Epilepsy_Zhang",col_names = FALSE)
autism_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","Autism",col_names = FALSE)
id_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","Intellectual_Disability",col_names = TRUE)
pd_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","Parkinson's",col_names = FALSE)
cancer_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","Cancer",col_names = FALSE)
hypertension_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","Hypertension",col_names = FALSE)
ra_df <- read_xlsx("Supp_table_5_disease_genes.xlsx","Rheumatoid_Arthritis",col_names = FALSE)

# Convert gene lists to character vectors
mono_ep_pg <- na.omit(as.character(mono_ep_pg_df$...1))
all_ep_pg <- na.omit(as.character(all_ep_pg_df$...1))
all_ep_zhang <- na.omit(as.character(all_ep_zhang_df$...1))
mono_ep_zhang <- na.omit(as.character(mono_ep_zhang_df$...1))
autism <- na.omit(as.character(autism_df$...1))
id <- na.omit(as.character(id_df$symbol))
pd <- na.omit(as.character(pd_df$...1))
cancer <- na.omit(as.character(cancer_df$...1))
hypertension <- na.omit(as.character(hypertension_df$...1))
ra <- na.omit(as.character(ra_df$...1))

##########################################################
# Function to extract genes from the CC gene list
process_gene_list <- function(df, gene_columns) {
  gene_list <- unlist(df[, gene_columns])
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  return(unique(gene_list))
}

# # Gene columns
gene_columns <- c("Gene.Name1", "Gene.Name2", "Gene.Name3", "Gene.Name4", "Gene.Name5")
# 
# # Create a list of unique genes - these are all the CC genes
ms_cc_all <- as.list(process_gene_list(df, gene_columns))
cc_all <- unlist(ms_cc_all)

# # Filter the dataframe for rows where fisher_pval_adj < 0.05 - now we'll get significant genes
df_significant <- df[df$fisher_pval_adj < 0.05, ]
# 
# # Create a list of unique significant genes
ms_cc_sig <-  as.list(process_gene_list(df_significant, gene_columns))
cc_sig <- unlist(ms_cc_sig)

# # Now let's do this for mild and severe separately
# # Filter for severe (fisher_log2fc > 0)
df_severe <- df_significant[df_significant$fisher_log2fc > 0, ]
ms_cc_severe <- as.list(process_gene_list(df_severe, gene_columns))
cc_severe <- unlist(ms_cc_severe)
# 
# # Filter for mild (fisher_log2fc < 0)
df_mild <- df_significant[df_significant$fisher_log2fc < 0, ]
ms_cc_mild <- as.list(process_gene_list(df_mild, gene_columns))
cc_mild <- unlist(ms_cc_mild)

# Function to convert human gene IDs to mouse
convert_human_to_mouse <- function(gene_list) {
  # Initialize empty dataframe with proper column names
  output <- data.frame(human_gene = character(),
                       mouse_gene = character(),
                       stringsAsFactors = FALSE)
  
  # Track genes with no matches
  missing_genes <- c()
  
  # Read MGI database
  mouse_human_genes <- read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep="\t")
  
  for(gene in gene_list) {
    class_key <- (mouse_human_genes %>% 
                    filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    
    if(!identical(class_key, integer(0))) {
      human_genes <- (mouse_human_genes %>% 
                        filter(DB.Class.Key == class_key & 
                                 Common.Organism.Name == "mouse, laboratory"))[,"Symbol"]
      
      if(length(human_genes) > 0) {
        for(mouse_gene in human_genes) {
          output <- rbind(output, 
                          data.frame(human_gene = gene,
                                     mouse_gene = mouse_gene,
                                     stringsAsFactors = FALSE))
        }
      } else {
        missing_genes <- c(missing_genes, gene)
      }
    } else {
      missing_genes <- c(missing_genes, gene)
    }
  }
  
  return(list(
    mapping = output,
    missing_genes = missing_genes
  ))
}

# Convert from human gene IDs to mouse
mono_ep_pg_converted <- convert_human_to_mouse(mono_ep_pg)
all_ep_pg_converted <- convert_human_to_mouse(all_ep_pg)
all_ep_zhang_converted <- convert_human_to_mouse(all_ep_zhang)
mono_ep_zhang_converted <- convert_human_to_mouse(mono_ep_zhang)
autism_converted <- convert_human_to_mouse(autism)
id_converted <- convert_human_to_mouse(id)
pd_converted <- convert_human_to_mouse(pd)
cancer_converted <- convert_human_to_mouse(cancer)
hypertension_converted <- convert_human_to_mouse(hypertension)
ra_converted <- convert_human_to_mouse(ra)

# Define color dictionary
color_dict <- list(
  cc_mild = "#058a92",  # teal
  cc_sig = "white",  # white
  cc_severe = "#530592"  # purple
)

#### Permutation and overlap function
analyze_gene_overlaps <- function(total_list, sig_list, disease_list, n_permutations, color_dict = NULL) {
  # Load required packages
  require(ggplot2)
  
  # Ensure inputs are character vectors
  total_list <- as.character(total_list)
  sig_list <- as.character(sig_list)
  disease_list <- as.character(disease_list)
  
  # Input validation
  if(length(sig_list) > length(total_list)) {
    stop("sig_list cannot be larger than total_list")
  }
  
  if(!all(sig_list %in% total_list)) {
    stop("All genes in sig_list must be present in total_list")
  }
  
  # Get the actual argument names from the call
  call_args <- match.call()
  sig_list_name <- deparse(call_args$sig_list)
  disease_list_name <- deparse(call_args$disease_list)
  
  # Calculate actual overlap between sig genes and mapping
  actual_overlap <- length(intersect(sig_list, disease_list))
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Perform permutations
  permuted_overlaps <- numeric(n_permutations)
  for(i in 1:n_permutations) {
    # Randomly sample the same number of genes as in sig_list
    random_genes <- sample(total_list, length(sig_list))
    # Calculate overlap with mapping
    permuted_overlaps[i] <- length(intersect(random_genes, disease_list))
  }
  
  # Calculate empirical p-value
  p_value <- sum(permuted_overlaps >= actual_overlap) / n_permutations
  
  # Calculate mean and standard deviation of permuted overlaps
  mean_overlap <- mean(permuted_overlaps)
  sd_overlap <- sd(permuted_overlaps)
  
  # Calculate z-score of actual overlap
  z_score <- (actual_overlap - mean_overlap) / sd_overlap
  
  # Create histogram
  overlap_df <- data.frame(overlaps = permuted_overlaps)
  
  subtitle_text <- paste(
    #paste("Comparing:", sig_list_name, sprintf("(n=%d)", length(sig_list))),
    #paste("With:", disease_list_name, sprintf("(n=%d)", length(disease_list))),
    paste("z-score =", round(z_score, 2)),
    paste("p-value =", round(p_value, 4)),
    sep = "\n"
  )
  
  # Use color from dictionary if provided, otherwise default to white
  fill_color <- if (!is.null(color_dict) && sig_list_name %in% names(color_dict)) {
    color_dict[[sig_list_name]]
  } else {
    "white"
  }
  
  # Calculate x-axis range and fixed offset for label
  x_range <- max(permuted_overlaps) - min(permuted_overlaps)
  label_offset <- x_range * 0.015  # 1.5% of the x-axis range
  
  # Calculate the extended x-axis limit to accommodate the label
  max_x <- max(permuted_overlaps) + (x_range * 0.3)  # Extend right side by 30% of range
  
  p <- ggplot(overlap_df, aes(x = overlaps)) +
    geom_histogram(binwidth = 1, fill = fill_color, color = "black") +
    # Add white background line
    geom_vline(xintercept = actual_overlap, color = "white", linewidth = 3) +
    geom_vline(xintercept = actual_overlap, color = "black", linetype = "dashed", linewidth = 2) +
    annotate("text", x = actual_overlap + label_offset, y = max(table(permuted_overlaps)),
             label = paste("Actual:", actual_overlap), 
             hjust = 0, angle = 0, size = 10) +
    labs(subtitle = subtitle_text, x="", y="")+
    theme_bw() +
    theme(
      text = element_text(family = "Helvetica"),  # Set default font family
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.text = element_text(size = 32),
      axis.title = element_text(size = 30),
      plot.subtitle = element_text(size = 30),
      axis.line.x = element_line(color = "black", linewidth = 2),
      axis.line.y = element_line(color = "black", linewidth = 2),
      axis.ticks.length = unit(7, "pt"),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")  # adds more margin on right side
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1))  # No expansion at bottom, 10% at top
    ) +
    scale_x_continuous(
      limits = c(min(permuted_overlaps), max_x),  # Extended x-axis limits
      expand = expansion(mult = c(0, 0))  # No expansion on either side
    )
  
  # Return results as a list
  results <- list(
    actual_overlap = actual_overlap,
    permuted_overlaps = permuted_overlaps,
    p_value = p_value,
    mean_permuted = mean_overlap,
    sd_permuted = sd_overlap,
    z_score = z_score,
    plot = p,
    input_sizes = list(
      total_genes = length(total_list),
      severe_genes = length(sig_list),
      mapping_genes = length(disease_list)
    )
  )
  
  return(results)
}

# Example usage for running mild gene overlaps with PG comprehensinve epilepsy list using 10,000 permutations
# Can run any combination
results <- analyze_gene_overlaps(
  total_list = cc_all,
  sig_list = cc_mild,
  disease_list = all_ep_pg_converted[["mapping"]][["mouse_gene"]],
  n_permutations = 10000,
  color_dict = color_dict
)

# Display the plot
print(results$plot)
# Save the plot
ggsave('plot.pdf')
