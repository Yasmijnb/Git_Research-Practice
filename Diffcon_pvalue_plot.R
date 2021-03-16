diffcon_plot <- function(res_path){
  #
  # Script for visualizing fatty acid significance in differential connectivity
  # 
  # Inputs:
  #       res_path: txt, location of differential connectivity analysis result files
  # Output:
  #       myplot: graph per differential connectivity result
  #
  # Author: Katarina Lepasepp
  
  require(ggplot2)
  require(janitor)
  require(ggrepel)
  
  merge.multiple <- function(x, ..., by = "row.names") {
    #
    # Function for merging multiple dataframes by rownames at once
    #
    # Input: 
    #       dataframes
    # Output:
    #       dataframe: merged dataframe based on rownames
    #
    dfs <- list(...)
    for (i in seq_along(dfs)) {
      x <- merge(x, dfs[[i]], by = by)
      rownames(x) <- x$Row.names
      x$Row.names <- NULL
    }
    return(x)
  }
  
  # Read in fatty acid levels for plotting
  fa_labels <- read.table('FAclass_subset.txt', check.names = F, header = T, sep='\t', row.names = 1)
  
  # Gather files from directory
  files <- list.files(pattern = "*.rds")
  
  for (file in files){
    # Read in rds file
    data <- readRDS(file)
    
    # Get specific comparison name
    comp <- sub("\\_diffcon_subset_1502_maxperm2000.rds$", "", file)
    print(comp)
    pvalues <- data$Pval_adj
    
    # Print statements to get better overview of results
    print(comp)
    print('These are adjusted p-values')
    print(pvalues)
    print('These are p-values')
    print(data$Pval)
    
    
    #  Calculate -log10 of adjusted p-values
    pvalues <- data.frame(log10(pvalues)) *(-1)
    diffcon <- data.frame(data$Diff_Conn)
    
    # Merge dataframes by rowname
    diffcon_data <- clean_names(merge.multiple(fa_labels, pvalues, diffcon))
    
    # Make rownames pretty in the case of subset analysis
    rownames(diffcon_data) <- c("Triglycerides_VLDL", "Triglycerides_IDL", 
                                "Triglycerides_LDL", "Triglycerides_HDL", 
                                "Cholesterol_VLDL", "Cholesterol_IDL", 
                                "Cholesterol_LDL", "Cholesterol_HDL",
                                "FreeCholesterol_VLDL", "FreeCholesterol_IDL",
                                "FreeCholesterol_LDL", "FreeCholesterol_HDL",
                                "Phospholipids_VLDL", "Phospholipids_IDL",
                                "Phospholipids_LDL", "Phospholipids_HDL",
                                "ApoA1_HDL","ApoA2_HDL", "ApoB_VLDL","ApoB_IDL",        
                                "ApoB_LDL")
    
    # Plot difference in connectivity with significance
    myplot <- ggplot(diffcon_data, aes(x = data_diff_conn, y = log10_pvalues, color = type)) +
      geom_point(aes(shape = type, color = type), size=6) + 
      geom_text_repel(aes(label=ifelse(log10_pvalues>0.01, as.character(rownames(diffcon_data)), '')), color='black', size=4, nudge_y = -0.06, hjust = 0.6) +
      scale_shape_manual(name = 'Lipid Type', values = c(15, 17, 19)) + scale_color_manual(name = 'Lipid Main Fractions', values=c("#0033cc", "#ffcc00", "#808080")) +
      labs(x = "Difference in Connectivity", y = '-log10(P-value)') +
      ggtitle(comp) + 
      ylim(min(diffcon_data$log10_pvalues)-0.1, max(diffcon_data$log10_pvalues)+0.1) + xlim(min(diffcon_data$data_diff_conn), max(diffcon_data$data_diff_conn)+0.25) +
      theme_bw()
    
    print(myplot)
    
  }
  
  
}

diffcon_plot("C:/Users/Yasmijn/Documents/School/WUR/SSB-80324 - Second Thesis/Git_Research-Practice/Results/")
