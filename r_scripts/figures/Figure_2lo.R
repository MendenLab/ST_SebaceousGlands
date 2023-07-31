# installed packages:
# BiocManager::install(c("ReactomePA", 'pathview', 'enrichplot', 'org.Hs.eg.db', 'DOSE'))
rm(list = ls(all=TRUE))
set.seed(1)

# R packages
rlibs <- c("dplyr", "stringr", "tibble", "xlsx", 'hash', 'pheatmap', 'VennDiagram', 'ggvenn')
invisible(lapply(rlibs, require, character.only = TRUE))
# Bioconductor packages
bioclibs <- c("ReactomePA", "pathview",  "enrichplot", "org.Hs.eg.db", 'DOSE', 'clusterProfiler')
invisible(lapply(bioclibs, require, character.only = TRUE))

# Plots
library(ggplot2)
library(forcats)

`%nin%` = Negate(`%in%`)


# Source R-scripts
source('/Users/christina.hillig/PycharmProjects/ST_SG_publication/ST_SebaceousGlands/r_scripts/spatialDE/init.R')
source('/Users/christina.hillig/PycharmProjects/ST_SG_publication/ST_SebaceousGlands/r_scripts/spatialDE/load_data.R')
source('/Users/christina.hillig/PycharmProjects/ST_SG_publication/ST_SebaceousGlands/r_scripts/spatialDE/utils.R')
# R-script to plot the ORA results
source("/Users/christina.hillig/PycharmProjects/ST_SG_publication/ST_SebaceousGlands/r_scripts/spatialDE/ORA_GSEA_plots.R")


# Date
date = "2023-04-12"
# Data set 
dataset.type = 'patient'
# Sequencing technique
seq.technique = "spatial"  
# Comparison
comparison = 'Patterns'
create_plot = FALSE

# Used design function:
design.function = "None"

# General input directory
input.dir = file.path('/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE', 
                      paste(date, '_paper_figures', sep =""))

# Cut parameter
l2fc.factor = 1
pval.cut_degs = 0.05
fdr.value = 0.2
pval.cut = 0.05
minGSSize = 5
# Multi-test method Benjamini-Hochberg (BH)
multitest.method =  "BH" 


# Plot Parameters
show_dotplot_categories = 15
show_categories = 10

width_img = 8
height_img = 12

selected.database = 'Reactome'


## Output
output.folder = file.path('/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/spatialDE', 
                          paste(date, '_paper_figures', sep =""), 'Pathway_enrichment_analysis',
                          'Figure2def')
# General output directory
output.dir <- get_outputdir(output.folder = output.folder, dataset.type = dataset.type,
                            func_task='ORA_REACTOME', func_comparison=comparison)

##### Load files
# 1. Get all DGE .csv files in subfolders
all_filenames = list.files(path = input.dir, pattern = c("*.xlsx"), recursive = TRUE, full.names=T)
# 2. read out filename
degs.name = all_filenames[grepl("*AEH*", all_filenames) & grepl("*corrected*", all_filenames)]
# Get unique specimens
# specimens <- sapply(str_split(sapply(str_split(degs.name, pattern='_'), "[[", 10), pattern='/'), "[[", 5)
specimens <- sapply(str_split(sapply(
  str_split(degs.name, pattern='_'), "[[", 12), pattern='/'), "[[", 2)

# File to check which genes were used for pathway analysis
used_degs = c()

specimen <- "11-V19T12-012-V2"


# Create specimen specific output directory
tmp_output.dir <- file.path(output.dir, specimen)
dir.create(tmp_output.dir, recursive = TRUE)

# Get filenames
filename = degs.name[grepl(paste(specimen, '*'), degs.name)]
# 4. Load spatialDE AEH result files with colnames: 
# g	pattern	membership	specimen	SEBACEOUS GLAND 
df.dge_res = load_file(path_name_file = filename)
df.dge_res$X <- NULL
df.dge_res$hgnc_symbol <- df.dge_res$g
df.dge_res$g <- NULL

# save as .xlsx file
result <- tryCatch({
  # The code I want run
  write.xlsx(df.dge_res, file.path(output.dir, 'AEH.xlsx'), sheetName=as.character(specimen), 
             row.names = FALSE, append = TRUE)
}, warning = function(war) {
  # Is executed if warning encountered
}, error = function(err) {
  # Is executed if error encountered - if file does not exist run this first
  write.xlsx(df.dge_res, file.path(output.dir, 'AEH.xlsx'), sheetName=as.character(specimen), 
             row.names = FALSE)
})


# Get unique patterns enriched in Sebaceous Glands
patterns <- unique(df.dge_res$pattern[df.dge_res$SEBACEOUS.GLAND == 1])

df.parameters <- data.frame(
  Parameters=c("DGE comparison", "", 'specimen', 'FDR_cut_PA',
               'p-value_cut_PA', 'Multitesting_method', 'minGSSize', 'database', '',
               'Preselected PAs', 'File_preselected_PAs', 'sheetName'),
  Values=c(file.path(input.dir,  filename), "", specimen,
           fdr.value, pval.cut, multitest.method, minGSSize, 'Reactome',
           '', 'no', '', 'CH_Reactome'))

# # Pathway enrichment analysis
# The pathway enrichment analysis included the following steps:
# 1. Identify significant and background genes 
# 1. Run PA analysis
# 1. Plot Save PA analysis results  
df.dge_res = rename_genetoentrezid(df.dge_results = df.dge_res) 


################### Get significantly DEx genes and background genes
enriched.pas <- hash::hash()
# Per pattern 
for (pattern.name in patterns) 
{
  # Read out only data of enriched patterns
  df.dge_res.patterns <- df.dge_res %>% filter_at(vars(pattern), any_vars(. %in% pattern.name))  
  
  # I. Define significant pattern genes and background genes
  # I.a) sort genes into groups belonging either to reference or test (control) condition
  df.ref_degenes = get_significantgenes(df.dge_results = df.dge_res.patterns, op = `==`) 
  df.ref = df.ref_degenes[[1]]
  
  # I.b) Background genes are all genes from our (sup-) data set
  bg_genes <- as.character(df.dge_res$entrezid)
  
  
  #################### Get Genes sets from Reactome ###############
  # ORA: Hypergeometric test
  enricher.ref <- ReactomePA::enrichPathway(
    gene=df.ref$entrezid, universe=bg_genes, organism = "human",
    pvalueCutoff=pval.cut, pAdjustMethod=multitest.method,
    minGSSize = minGSSize, maxGSSize = 500,
    readable=T)
  
  
  #############################################################################
  ################### ---> Save results to csv file <--- ###################### 
  #############################################################################
  ref_group = strsplit(comparison, '_')[[1]][1]  # enriched in Sebaceous Glands
  
  # Save result per pattern enriched in Sebaceous Glands
  # TODO save in another sheet in same excel file
  save_enrichobject_as_csv(
    paenrich_object = enricher.ref, comparison=specimen,
    condition = ref_group, df.parameters=df.parameters, pattern=pattern.name,
    pa_database = selected.database, method='ORA', 
    output_path = tmp_output.dir)
  
  #############################################################################
  ###################           ---> Plots <---          ###################### 
  #############################################################################
  # 2. Barplot
  if (!is.null(nrow(enricher.ref)))
  {
    if (nrow(enricher.ref) > 1 & length(unique(enricher.ref@result$ID)) > 1) 
    {
      def.plot.enricher_barplot(
        enricher.result=enricher.ref, showCategories=15, method=paste('Pattern', pattern.name), 
        title=paste(selected.database,ref_group, specimen, pattern.name,
                    'ORA_barplot.pdf', sep = "_"),
        width=width_img, height=height_img, output.dir=tmp_output.dir, func_colours=NULL) 
    }
  }
  
  used_degs <- c(used_degs, df.dge_res.patterns$hgnc_symbol)
  enriched.pas[[as.character(pattern.name)]] <- list(enricher.ref@result$Description[
    enricher.ref@result$p.adjust < 0.05])
}

sig.pathways.pattern <- list()
names.pathways <- c()
for (pattern.name in patterns) 
{
  tmp.pathway.names <- enriched.pas[[as.character(pattern.name)]]
  # only for patterns with enriched Pathways
  if (length(tmp.pathway.names[[1]]) > 1) 
  {
    names.pathways <- c(names.pathways, as.character(pattern.name))
    names(tmp.pathway.names) <- as.character(pattern.name)
    sig.pathways.pattern <- c(sig.pathways.pattern, tmp.pathway.names)
  }
}

# Rearrange order of list objects in list
sig.pathways.pattern.v2 <- list(sig.pathways.pattern$`0`, sig.pathways.pattern$`8`, sig.pathways.pattern$`6`)
names(sig.pathways.pattern.v2) <- c('0', '8', '6')

pdf(file = file.path(tmp_output.dir, paste0("VennDiagram_enriched_pathways", specimen, ".pdf")),   
    width = 6, # The width of the plot in inches
    height = 4) # The height of the plot in inches

ggvenn.plot <- ggvenn::ggvenn(
  sig.pathways.pattern.v2, text_size = 4,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 6) 

print(ggvenn.plot)
dev.off()

# Save shared und unique pathways
n.obs <- sapply(sig.pathways.pattern, length)
df <- list2DF(lapply(sig.pathways.pattern, `length<-`, max(n.obs)))
write.table(used_degs, file.path(
  tmp_output.dir, paste0(specimen, '_ORA_REACTOME_significant_pathways.txt', sep='')), 
  append = FALSE, sep = " ", dec = ".")


v.table <- venn::venn(sig.pathways.pattern)
shared.pathways.08 <- attr(v.table, "intersections")[['0:8']]

# Save used genes per specimen enriched in Sebaceous Glands
write.table(used_degs, file.path(
  tmp_output.dir, paste0(specimen, '_ORA_REACTOME.txt', sep='')), 
  append = FALSE, sep = " ", dec = ".")


