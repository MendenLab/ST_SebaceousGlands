# Installation for multienrichjam
# remotes::install_github("jmw86069/jamba")
# remotes::install_github("jmw86069/multienrichjam")
rm(list = ls())

library('xlsx')
library('clusterProfiler')
library(ggplot2)
library(cowplot)

# ------------------------------ Parameters
lfc_factor = 1
fdr_value = 0.1
pval_cut = 0.1
p_value = 0.05
minGSSize = 10
test_method =  "BH" 

# ------------------------------ Create output directory
save_dir <- file.path('/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer', 
                      "output", "SEBG_Figure_4c", Sys.Date())
dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)
input.dir.files <- '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/input/figure_5c'

# ------------------------------ Load Files
df.pa <- xlsx::read.xlsx(file.path(
  input.dir.files,'Pathway_enrichment_analysis', 'ST_cdr_annotation_condition', '20210904_SEBG_L_diseases', 
  '5_DGE_Approaches_Glands_modifiedPS_20210507', 'SEB G, L, disease__condition 1_vs_condition 2', 
  'condition 1Reference_positiveLog2FC_REACTOME_Pathway_Enrichment_Analysis_.xlsx'), sheetIndex = 1)

df.dge <- xlsx::read.xlsx(file.path(
  input.dir.files, 'DGE_analysis', 'ST_cdr_annotation_condition', '20210904_SEBG_L_diseases',
  '5_DGE_Approaches_Glands_modifiedPS_20210507', 'SEB G, L, disease__condition 1_vs_condition 2', 
  'ST_condition 1_vs_condition 2_glmGamPoi_DGE_all_genes.xlsx'), sheetIndex = 1)


if (file.exists(file.path(input.dir.files, "Enrich_obj.rds"))) 
{
  enrich.obj <- readRDS(file.path(input.dir.files, "Enrich_obj.rds"))
} else 
{
  library('multienrichjam')
  
  # ------------------------------ Convert df to enrich object
  enrich.obj <- multienrichjam::enrichDF2enrichResult(
    enrichDF = df.pa,
    pvalueCutoff = pval_cut,
    pAdjustMethod = test_method,
    keyColname = "ID",
    geneColname = "geneID",
    geneHits = "Count",
    geneRatioColname = "GeneRatio",
    geneDelim = "[,/ ]+",
    pvalueColname = 'pvalue',
    descriptionColname = "Description",
    msigdbGmtT = NULL,
    verbose = FALSE,
  )
  
  saveRDS(enrich.obj, file = file.path(input.dir.files, "Enrich_obj.rds"))
}


# ------------------------------ Get Log2FC of genes
ind_genes <- which(df.dge$gene_symbol %in% enrich.obj@gene)
sub_df.dge <- df.dge[ind_genes, ]
# Order df like in enrich object
sub_df.dge <- sub_df.dge[match(enrich.obj@gene, sub_df.dge$gene_symbol),]
entrezid_log2fc <- sub_df.dge$log2fc
names(entrezid_log2fc) <- sub_df.dge$gene_symbol



showCategories <- c('Interferon Signaling', 'Interferon gamma signaling',
                    'Antimicrobial peptides', 'Keratinization', 'SUMOylation')
showCategories = intersect(enrich.obj$Description, showCategories)

# Plot cnet MSLN
png(
  filename = file.path(
    save_dir, '5_DGE_Approaches_Glands_modifiedPS_20210507__SEB G, PSO, L__Cnetplot.png'),
  width = 25, height = 6, units = "in", bg = "white", res = 300)
p1 <- clusterProfiler::cnetplot(enrich.obj, showCategory = showCategories, 
                          categorySize = "pvalue", foldChange = entrezid_log2fc, 
                          colorEdge = TRUE) + 
  ggtitle("REACTOME: Pathway Enrichment Analysis using DB Reactome") 
# Color upper and lower border
min.value <- floor( min(p1$data$color, na.rm = TRUE) )
max.value <- ceiling( max(p1$data$color, na.rm = TRUE) )

p1 <- p1 + scale_color_gradientn(
  name = "fold change", colours = c("blue", "red"), limits= c(min.value, max.value))
cowplot::plot_grid(p1, ncol = 1, nrow = 1)

print(p1)
dev.off()

