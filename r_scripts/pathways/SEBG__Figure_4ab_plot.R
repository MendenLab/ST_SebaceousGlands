rm(list = ls())

# Installation for multienrichjam
# remotes::install_github("jmw86069/jamba")
# remotes::install_github("jmw86069/multienrichjam")

library('xlsx')
library(dplyr)

# Plotting libraries
library('clusterProfiler')
library(ggplot2)
library(cowplot)

def.get_genes_log2fc <- function(func_df.dge, func_enrich.obj) 
{
  ind_genes <- which(func_df.dge$gene_symbol %in% func_enrich.obj@gene)
  sub_df.dge <- func_df.dge[ind_genes, ]
  # Order df like in enrich object
  sub_df.dge <- sub_df.dge[match(func_enrich.obj@gene, sub_df.dge$gene_symbol),]
  # create named vector with ENTREZ IDs as names and log2FC as entries
  func_entrezid_log2fc <- sub_df.dge$log2fc
  names(func_entrezid_log2fc) <- sub_df.dge$gene_symbol
  
  return(func_entrezid_log2fc)
}

def.plot_cnet <- function(func_enrich.obj, func_categories, func_entrezidlog2fc, 
                          output.dir, func_title, width=16, height=8) 
{
  png(
    filename = file.path(output.dir, func_title),
    width = width, height = height, units = "in", bg = "white", res = 300)
  p1 <- clusterProfiler::cnetplot(
    func_enrich.obj, showCategory = func_categories, 
    categorySize = "pvalue", color.params = list(foldChange = func_entrezidlog2fc), colorEdge = TRUE) + 
    ggtitle("REACTOME: Pathway Enrichment Analysis using DB Reactome") 
  # Color upper and lower border
  min.value <- floor( min(p1$data$color, na.rm = TRUE) )
  max.value <- ceiling( max(p1$data$color, na.rm = TRUE) )
  
  p1 <- p1 + scale_color_gradientn(
    name = "fold change", colours = c("blue", "red"), limits= c(min.value, max.value))
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}


# 1. Create output directory
save_dir <- file.path('/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer', 
                      "output", "SEBG_Figure_4ab", Sys.Date())
dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)

# 2. Load .csv files
input.dir.files <- file.path('/Volumes', 'CH__data', 'Projects', 'Eyerich_AG_projects', 
                             'ST_Sebaceous_glands__Peter_Seiringer', 'input', 'figure_5c')

file.name <- 'condition 1 Reference_positiveLog2FC_REACTOME_Pathway_Enrichment_Analysis_.xlsx'

comparison.a <- xlsx::read.xlsx(file.path(
  input.dir.files, 'Pathway_enrichment_analysis', 'ST_cdr_patient_annotation_condition', 
  '20210616_LvsNL', 'SEB G, L, AE vs SEB G, NL__condition 1_vs_condition 2',
  file.name), sheetIndex = 1, header = TRUE, Comparison='SEB G, AE, L')
# keep only enriched pathways p.adjust <= 0.1
comparison.a$p.adjust <- as.numeric(comparison.a$p.adjust)
comparison.a <- comparison.a[comparison.a$p.adjust <= 0.1, ]
# sort by p.adj
comparison.a <- comparison.a[order(comparison.a$p.adjust),]


comparison.b <- xlsx::read.xlsx(file.path(
  input.dir.files, 'Pathway_enrichment_analysis', 'ST_cdr_patient_annotation_condition', 
  '20210616_LvsNL', 'SEB G, L, PSO vs SEB G, NL__condition 1_vs_condition 2', 
  file.name), sheetIndex = 1, header = TRUE,  Comparison='SEB G, Pso, L ')
# keep only enriched pathways p.adjust <= 0.1
comparison.b$p.adjust <- as.numeric(comparison.b$p.adjust)
comparison.b <- comparison.b[comparison.b$p.adjust <= 0.1, ]
# sort by p.adj
comparison.b <- comparison.b[order(comparison.b$p.adjust),]

print(head(comparison.a))
print(dim(comparison.a))
print(head(comparison.b))
print(dim(comparison.b))

# Select top 15 per comparison  <- discussed with Peter
sub_comparison.a <- comparison.a[1:15, ]
sub_comparison.b <- comparison.b[1:15, ]

# create combined dataframe join by `Comparison` column
# Way 1
# joined.df <- rbind(comparison.a, comparison.b, by='Description')
# Way 2
joined.df <- rbind(sub_comparison.a, sub_comparison.b, by='Description')
print(head(joined.df))
print(dim(joined.df))
joined.df$p.adjust <- as.numeric(joined.df$p.adjust)
joined.df <- joined.df[order(joined.df$p.adjust),]
# remove last row 
joined.df <- joined.df[-c(dim(joined.df)[1]), ]

# Plot top 15
joined.df$GeneRatio <- as.numeric(sapply(joined.df$GeneRatio, function(x) eval(parse(text=x))))

sub_joined.df <- joined.df[c(1:15), ]
# sort by GeneRatio
# sub_joined.df <- sub_joined.df[order(sub_joined.df$GeneRatio, decreasing = TRUE), ]

sub_joined.df <- sub_joined.df %>% filter(if_any(everything(), ~ !is.na(.)))

# create interval for count
sub_joined.df$Count <- as.numeric(sub_joined.df$Count)
breaks <- seq(min(sub_joined.df$Count), max(sub_joined.df$Count), 10)
sub_joined.df$Counts <- cut(sub_joined.df$Count, c(breaks, Inf), 
                            labels = breaks, include.lowest = TRUE)

# Workaround to show categoires not in aplhabetical order
sub_joined.df$Description <- factor(
  sub_joined.df$Description, levels=unique(sub_joined.df$Description))
sub_joined.df$Description <- factor(sub_joined.df$Description, levels=sub_joined.df$Description)

# https://uc-r.github.io/cleveland-dot-plots
# create grouped Barplot
png(
  filename = file.path(
    save_dir, 
    paste('Barplot', "20210616_LvsNL__SEB_G_AE_or_Pso_L.png", sep = '_')),
  width = 16, height = 8, units = "in", bg = "white", res = 300)
p1 <- ggplot(sub_joined.df, aes(Description, GeneRatio, fill = p.adjust)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + facet_wrap(~ Comparison) + theme(text = element_text(size=14)) +
  scale_fill_gradient(low="blue", high="red", trans = 'reverse')
print(p1)
dev.off()


sub_joined.df_workflow <- sub_joined.df[c(4:6), ]
png(
  filename = file.path(
    save_dir, 
    paste('Barplot', "Workflow__20210616_LvsNL__SEB_G_AE_or_Pso_L.png", sep = '_')),
  width = 8, height = 6, units = "in", bg = "white", res = 300)
p1 <- ggplot(sub_joined.df_workflow, aes(Description, GeneRatio, fill = p.adjust)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggplot2::theme_minimal() +
  coord_flip() + facet_wrap(~ Comparison) + theme(text = element_text(size=20)) +
  scale_fill_gradient(low="blue", high="red", trans = 'reverse')
print(p1)
dev.off()


# create grouped dotplot 
png(
  filename = file.path(
    save_dir, 
    paste('Dotplot', "20210616_LvsNL__SEB_G_AE_or_Pso_L.png", sep = '_')),
  width = 16, height = 8, units = "in", bg = "white", res = 300)
p2 <- ggplot(sub_joined.df, aes(Description, GeneRatio)) +
  geom_point(aes(size=Counts, color = p.adjust)) + facet_wrap(~ Comparison) +
  coord_flip()  + theme(text = element_text(size=16)) +
  scale_x_discrete(limits=rev) + 
  scale_color_gradient(low="blue", high="red", trans = 'reverse')
print(p2)
dev.off()


# ------------------------------------------------------------------------ #
# ---------------------           Cnet plot         ---------------------- #
# ------------------------------------------------------------------------ #

# ------------------------------ Parameters
lfc_factor = 1
fdr_value = 0.1
pval_cut = 0.1
p_value = 0.05
minGSSize = 10
test_method =  "BH" 

# ------------------------------ Load DGE File
dges_file_names <- 'ST_condition 1_vs_condition 2_glmGamPoi_DGE_all_genes.csv'
df.dge.comparison_a <- read.csv(
  file.path(input.dir.files, 'DGE_analysis', 'ST_cdr_patient_annotation_condition', 
            '20210616_LvsNL', 'SEB G, L, AE vs SEB G, NL__condition 1_vs_condition 2', dges_file_names))
df.dge.comparison_b <- read.csv(
  file.path(input.dir.files, 'DGE_analysis', 'ST_cdr_patient_annotation_condition', 
            '20210616_LvsNL', 'SEB G, L, PSO vs SEB G, NL__condition 1_vs_condition 2', dges_file_names))


if (file.exists(file.path(input.dir.files, "Enrich_obj_SEB G, L, AE vs SEB G, NL.rds")) & 
    file.exists(file.path(input.dir.files, "Enrich_obj_SEB G, L, PSO vs SEB G, NL.rds")))  
{
  enrich.obj.comparison_a <- readRDS(file.path(input.dir.files, "Enrich_obj_SEB G, L, AE vs SEB G, NL.rds"))
  enrich.obj.comparison_b <- readRDS(file.path(input.dir.files, "Enrich_obj_SEB G, L, PSO vs SEB G, NL.rds"))
  
  enrich.obj.comparison_a@result$p.adjust <- as.numeric(enrich.obj.comparison_a@result$p.adjust)
  enrich.obj.comparison_b@result$p.adjust <- as.numeric(enrich.obj.comparison_b@result$p.adjust)
  enrich.obj.comparison_a@result$pvalue <- as.numeric(enrich.obj.comparison_a@result$pvalue)
  enrich.obj.comparison_b@result$pvalue <- as.numeric(enrich.obj.comparison_b@result$pvalue)
  
} else 
{
  library('multienrichjam')
  
  # ------------------------------ Convert df to enrich object
  enrich.obj.comparison_a <- enrichDF2enrichResult(
    enrichDF = comparison.a,
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
  
  enrich.obj.comparison_b <- enrichDF2enrichResult(
    enrichDF = comparison.b,
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
  
  saveRDS(enrich.obj.comparison_a, file = file.path(input.dir.files, "Enrich_obj_SEB G, L, AE vs SEB G, NL.rds"))
  saveRDS(enrich.obj.comparison_b, file = file.path(input.dir.files, "Enrich_obj_SEB G, L, PSO vs SEB G, NL.rds"))
}


# ------------------------------ Get Log2FC of genes
entrezid_log2fc_com.a <- def.get_genes_log2fc(
  func_df.dge=df.dge.comparison_a, func_enrich.obj=enrich.obj.comparison_a) 
entrezid_log2fc_com.b <- def.get_genes_log2fc(
  func_df.dge=df.dge.comparison_b, func_enrich.obj=enrich.obj.comparison_b) 


# Plot cnet 
# combined categories:
# please combine „respiratory electron transport“ and „respiratory electron transport, ATP synthesis, …“
# add 'metabolism of steroids'
showCategories <- c('The citric acid (TCA) cycle and respiratory electron transport',
                    'Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.', 
                    'Fatty acid metabolism', 'Cholesterol biosynthesis',
                    'Metabolism of steroids')
def.plot_cnet(func_enrich.obj=enrich.obj.comparison_a, func_categories=showCategories, 
              func_entrezidlog2fc=entrezid_log2fc_com.a, output.dir=save_dir, 
              func_title='20210616_LvsNL__SEB G, L, AE vs SEB G, NL__condition 1_vs_condition 2_Cnetplot.png', 
              width = 20, height=6) 
# show top 5
# showCategories <- c('Keratinization',
#                     'Acyl chain remodelling of PI', 
#                     'Neutrophil degranulation', 'Antimicrobial peptides"',
#                     'Metabolism of steroids')
def.plot_cnet(func_enrich.obj=enrich.obj.comparison_b, func_categories=5, 
              func_entrezidlog2fc=entrezid_log2fc_com.b, output.dir=save_dir, 
              func_title='20210616_LvsNL__SEB G, L, Pso vs SEB G, NL__condition 1_vs_condition 2_Cnetplot.png', 
              width = 20, height=6) 
  




