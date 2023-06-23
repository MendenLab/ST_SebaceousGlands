#! /usr/bin/Rscript

# R version 4.0.3 Patched (2020-10-23 r79366)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7

# install Packages: BiocManager::install("reactome.db"), BiocManager::install("pathview"), install.packages("cowplot"), install.packages("xlsx")
# BiocManager::install("enrichplot"), BiocManager::install("ReactomePA"), install.packages("qdapTools"),  BiocManager::install("clusterProfiler"), install.packages("Hmisc")
########################################################################################

# Analysis included in that script:
# Reactome pathway enrichment analysis 

########################################################################################

# remove all variables in global environment
rm(list = ls())
library(docstring)

# libraries
# for pathway enrichemnt analysis:
library(ReactomePA)

# Gene, GO-term and Pathway database:
library(org.Hs.eg.db) # human organism

# Plot for pathway enrichemnt analysis:
library(ggplot2)
library(enrichplot) 
library(cowplot)

# variable or objects functions:
library(dplyr) 
library(Hmisc)
library(pathview)
library("xlsx")
library("qdapTools")

# Plot functions
source("/Users/christina.hillig/R_studio/NGS_Pathway_analysis/scripts/pathway_analysis/PA_plotting_functions.R")
source("/Users/christina.hillig/R_studio/NGS_Pathway_analysis/scripts/pathway_analysis/ImmunePublication_pathways.R")
source('/Users/christina.hillig/R_studio/NGS_Pathway_analysis/scripts/utils.R')
source('/Users/christina.hillig/R_studio/NGS_Pathway_analysis/scripts/pathway_analysis/Sebaceous_glands_pathways.R')


#################################### -> functions <- ##################################
load_files <- function(path_name_file)
{
  #' Read in .csv file 
  #' 
  #' Gets the full path and reads in csv file including the header
  #' @note The .csv file entries should be separated by ","
  #' @param path_name_file Path to .csv file
  #' @return Dataframe
  
  if (tail(strsplit(path_name_file, split = "[.]")[[1]], n = 1) == "xlsx") 
  {
    dge_df <- read.xlsx(path_name_file, sheetIndex = 1, header = TRUE)
  } else {
    dge_df <- read.csv(path_name_file, header = TRUE, sep = ";")
    
    if (ncol(dge_df)==1) {dge_df = read.csv(path_name_file, header = TRUE, sep = ",")}
  }
  
  return(dge_df)
}


############## initialization phase
main = function(sample, date_file, replicate_type, dge_approach, minGSSize, dge.method,
                lfc_factor, fdr_value, p_value_degs, pval_cut, correction_method) 
{
  print("-------------------->  Pathway enrichment  <---------------------")
  
  ############################### 1. Initialization phase
  sub_folders <- list.dirs(
    path = file.path(getwd(), "input", paste(dge_approach, "output", sep = "_"), dge.method,
                     date_file, sample), full.names = TRUE, recursive = FALSE)
  
  # output paths
  dir.create(file.path(getwd(), "output", "PA_Figure7_output"),  recursive = TRUE,
             showWarnings = FALSE)
  dir.create(file.path(getwd(), "output", "PA_Figure7_output", Sys.Date()),  recursive = TRUE,
             showWarnings = FALSE)
  save_dir <- file.path(getwd(), "output", "PA_Figure7_output", Sys.Date(),
                        paste(dge_approach, "output", sep = "_"))
  dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)
  
  
  ######### START ANALYSIS
  for (sub_folder in sub_folders) 
  {
    # get all DGE .csv files in subfolders
    all_filenames = list.files(path = sub_folder, pattern = c("*.csv|*.xlsx"), 
                               recursive = TRUE)
    for (dge_filename in all_filenames[!grepl("metaData*", all_filenames)])
    {
      print(dge_filename)
      # Get path to file
      dir_path = dirname(dge_filename)
      
      # ----- Get name of ref.cond and ctrl.cond 
      # identify reference condition 
      ref.cond = strsplit(tail(strsplit(dge_filename, split = "vs_", perl = TRUE)[[1]], n = 1),
                          split = "[.]")[[1]][1]
      # get control condition 
      ctrl.cond = tail(strsplit(strsplit(dge_filename, split = "_vs_", perl = TRUE)[[1]][1],
                                split = '__')[[1]],  n = 1)
      # -----
      
      name_comparison <- strsplit(strsplit(dge_filename, "DEG__")[[1]][1], split = "[.]")[[1]][1]
      results_save_dir <- file.path(
        save_dir, sample, 
        tail(strsplit(sub_folder, .Platform$file.sep)[[1]],  n = 1), 
        strsplit(strsplit(dge_filename, .Platform$file.sep)[[1]], "[.]")[[1]][1])
      dir.create(results_save_dir, recursive = TRUE, showWarnings = FALSE)
      
      # path to csv files
      dge_results_load_dir <- file.path(sub_folder,  dge_filename)

      # Load dge list with gene names and p-value (and/or log2FC)
      dge_list_results <- load_files(path_name_file = dge_results_load_dir)
      # remove duplicated rows ..
      dge_list_results <- dge_list_results[!duplicated(dge_list_results$gene_symbol), ]
      rownames(dge_list_results) <- dge_list_results$gene_symbol
      dge_list_results <- dge_list_results[, c(colnames(dge_list_results) != "X"),
                                           drop = FALSE]
      # remove NA columns
      dge_list_results <- dge_list_results[, colSums(is.na(dge_list_results)) < 
                                             nrow(dge_list_results)]
      # remove NA rows
      dge_list_results <- na.omit(dge_list_results)
      
      ########################### ---> Preparation <--- ##########################
      # I. replace gene symbol by entrezid of gene
      if ("entrezid" %nin% colnames(dge_list_results))
      {
        # works approach 
        Entrez_ID <- mapIds(org.Hs.eg.db, row.names(dge_list_results), 'ENTREZID', 'SYMBOL')
        Entrez_ID_df = qdapTools::list2df(Entrez_ID, col1 = "entrezid", col2 = "gene_symbol")
        Entrez_ID_df = Entrez_ID_df[!is.na(Entrez_ID_df$entrezid), ]
        
        gseaDat <- filter(dge_list_results, !is.na(Entrez_ID))
        dge_list_results <- filter(dge_list_results, !is.na(Entrez_ID))
        dge_list_results$entrezid <- Entrez_ID_df$entrezid
      }
      
      # I. Rank all genes based on their fold change
      ranks <- dge_list_results$log2fc 
      names(dge_list_results)[names(dge_list_results) == 'log2fc'] <- 'log2fc'
      
      names(ranks) <- dge_list_results$gene_symbol #dge_list_results$entrezid
      ranked_genes <- sort(ranks, decreasing = T)
      
      # III. Define significant DEx genes and background genes
      # 1.a) sort genes into groups belonging either to reference or test (control) condition
      # positive Log2FC 
      df.ctrl <- dge_list_results[
        dge_list_results$pval < p_value_degs & !is.na(dge_list_results$pval) & 
          dge_list_results$log2fc > lfc_factor, ]
      # sort dataframe by log2fc
      df.ctrl = df.ctrl[order(df.ctrl$log2fc, decreasing = F), ]
      de.ctrl <- df.ctrl$entrezid
      de.ctrl <- as.character(na.exclude(de.ctrl))
      ranked_genes.ctrl <- df.ctrl$log2fc
      names(ranked_genes.ctrl) <- df.ctrl$gene_symbol
      ranked_genes.ctrl <- sort(ranked_genes.ctrl, decreasing = F)
      
      # negative Log2FC 
      df.ref <- dge_list_results[
        dge_list_results$pval < p_value_degs & !is.na(dge_list_results$pval) & 
          dge_list_results$log2fc < -lfc_factor, ]
      # sort dataframe by log2fc
      df.ref = df.ref[order(df.ref$log2fc, decreasing = T), ]
      de.ref <- df.ref$entrezid
      de.ref <- as.character(na.exclude(de.ref))        
      ranked_genes.ref <- df.ref$log2fc
      names(ranked_genes.ref) <- df.ref$gene_symbol
      ranked_genes.ref <- sort(ranked_genes.ref, decreasing = T)
      
    
      # b.) Background genes are all genes from our (sup-) data set
      bg_genes <- as.character(dge_list_results$entrezid)
      
      # transform entrez_id to factor
      dge_list_results$entrezid = as.factor(dge_list_results$entrezid)
      
      #############################################################################
      #################### ---> Pathway Enrichment Analysis <--- ##################
      #############################################################################
      # 2. ReactomePA Pathway enrichment analysis of a gene set: Database REACTOME
      # Note: Used to determine occuring protein receptors in dataset
      # Given vector of genes, function returns enriched pathways with FDR control.
      # 2.1 Find enriched Pathways for reference condition
      reactome_object.ref <- enrichPathway(gene = de.ref, # a vector of entrez gene id
                                           universe = bg_genes,
                                           organism = 'human', 
                                           qvalueCutoff = fdr_value, 
                                           pvalueCutoff = pval_cut, 
                                           pAdjustMethod = correction_method, 
                                           minGSSize = minGSSize,
                                           maxGSSize = 500,
                                           readable = T)
      
      # 2.2 Find enriched Pathways for control condition
      reactome_object.ctrl <- enrichPathway(gene = de.ctrl, # a vector of entrez gene id
                                            universe = bg_genes,
                                            organism = 'human', 
                                            qvalueCutoff = fdr_value, 
                                            pvalueCutoff = pval_cut, 
                                            pAdjustMethod = correction_method, 
                                            minGSSize = minGSSize,
                                            maxGSSize = 500,
                                            readable = T)
      
      ################### ---> convert gene ID to Symbol <--- ################### 
      # Pathway Enrichment
      reactome.ctrl = setreadable_pa(paenrich_object=reactome_object.ctrl)
      reactome.ref = setreadable_pa(paenrich_object=reactome_object.ref) 
      
      
      #############################################################################
      ################### ---> Save results to csv file <--- ###################### 
      #############################################################################
      df.parameters <- data.frame(
        Parameters=c("DGE comparison", "", 'Log2FC_cut_DEGs', 'p-value_cut_DEGs', 'FDR_cut_PA',
                     'p-value_cut_PA', 'Multitesting_method', 'minGSSize'),
        Values=c(dge_filename, "", lfc_factor, p_value_degs, fdr_value, pval_cut,
                 correction_method,
                 minGSSize))
      if (!is.null(nrow(reactome.ctrl))) 
      {
        enrichobject_to_df(paenrich_object=reactome.ctrl, 
                           condition=paste0(ctrl.cond, '_Control_positiveLog2FC'), 
                           comparison=strsplit(dge_filename,  split = "[.]")[[1]][1],
                           pa_database='REACTOME', method='PA', output_path=results_save_dir,
                           df.parameters=df.parameters) 
      }

      if (!is.null(nrow(reactome.ref))) 
      {
        enrichobject_to_df(paenrich_object=reactome.ref, 
                           condition=paste0(ref.cond, '_Reference_negativeLog2FC'), 
                           comparison=strsplit(dge_filename,  split = "[.]")[[1]][1],
                           pa_database='REACTOME', method='PA', output_path=results_save_dir,
                           df.parameters=df.parameters) 
      }
      #############################################################################
      ############################### ---> PLOTS <--- ############################# 
      #############################################################################
      # Plot variables
      # select pathways or Enriched gene sets manually
      show_categories = 5
      show_dotplot_categories =  15 
      
      width_img = 16
      height_img = 8
      
      ######### ---> Save Pathway Enrichment Analysis Plots and Files <--- #########
      # If a gene is associated with two or more enriched PAs 
      # but less than those are shown than this results in a bug 
      # -> the log2fc of that gene is not correctly shown
      if (!is.null(nrow(reactome.ref)))
      {
        if (nrow(reactome.ref) > 0) # & any(show_categories %in% reactome.ref$Description)) 
        {
          # Cnetplots to visualise enriched pathways
          png(
            filename = file.path(
              results_save_dir, 
              paste0(ref.cond, "_Reference_negativeLog2FC_", 
                     strsplit(dge_filename,  split = "[.]")[[1]][1], 
                     "_REACTOME_PA_Cnetplot.png")),
            width = 20, height = 8, units = "in", bg = "white", res = 300)
          print(fig.pathways.REACTOME(reactome_res = reactome.ref, 
                                      id=reactome.ref@result$ID[1:5],
                                      entrezid_log2fc = ranked_genes.ref,
                                      showCategories = show_categories))
          dev.off()
        }
      }
      
      if (!is.null(nrow(reactome.ctrl)))
      {
        if (nrow(reactome.ctrl) > 0) # & any(show_categories %in% reactome.ctrl$Description)) 
        {
          # Cnetplots to visualise enriched pathways
          png(
            filename = file.path(
              results_save_dir, 
              paste0(ctrl.cond, "_Control_positiveLog2FC_", 
                     strsplit(dge_filename,  split = "[.]")[[1]][1], 
                     "_REACTOME_PA_Cnetplot.png")),
            width = 20, height = 8, units = "in", bg = "white", res = 300)
          print(fig.pathways.REACTOME(reactome_res = reactome.ctrl, 
                                      id=reactome.ctrl@result$ID[1:5],
                                      entrezid_log2fc = ranked_genes.ctrl,
                                      showCategories = show_categories))
          dev.off()
        }
      }
    }
  }
  sessionInfo()
}

date_file <- "2021-12-10" 
replicate_type = "biological"
sample = "Bulk_unknown"
dge_approach = "Sebaceous glands"
dges_author = 'Daniel__figure7'

lfc_factor = 1
fdr_value = 0.1
pval_cut = 0.1
p_value = 0.1
minGSSize = 10

# FDR and BH are more conservative than bonferroni
test_method =  "BH" 

# Structure of filenames
# 1. Condition
# 2. DGE comparison
# 3. design matrix
# 4. Database: KEGG, Reactome, ..
# 5. Method: PA, GSEA, ORA, ..
# 6. Plot type (optional)

main(date_file = date_file, sample = sample, dge_approach = dge_approach, 
     replicate_type = replicate_type, dge.method = dges_author,
     lfc_factor = lfc_factor, fdr_value = fdr_value, p_value_degs = p_value, pval_cut = pval_cut,
     correction_method = test_method, minGSSize=minGSSize)
