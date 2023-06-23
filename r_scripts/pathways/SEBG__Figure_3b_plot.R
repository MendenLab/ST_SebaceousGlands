rm(list = ls())

library(xlsx)
library(dplyr)
library(tidyverse)

# Plotting libraries
library(ggplot2)


# 1. Create output directory
save_dir = file.path('/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output', 
                          "SEBG_Figure3b", Sys.Date())
dir.create(save_dir, showWarnings = FALSE,  recursive = TRUE)

# 2.1Load .csv files of 20210830_Tasksforfigures/NL, SEB G vs NL, rest of skin
input.dir <- '/Users/christina.hillig/R_studio/NGS_Pathway_analysis/output/GSEA_output/2021-09-07/Sebaceous glands_output/ST_cdr_patient_annotation_condition'

file.name <- 'condition 1Reference_positiveLog2FC_REACTOME_Pathway_Enrichment_Analysis_.xlsx'

df_sebg <- xlsx::read.xlsx(file.path(
  input.dir, '20210830_Tasksforfigures/NL, SEB G vs NL, rest of skin__condition 1_vs_condition 2/ST_condition 1_vs_condition 2_glmGamPoi_DGE_all_genes.csv',
  file.name), sheetIndex = 1, header = TRUE, Comparison='NL, SEB G')
# keep only enriched pathways p.adjust <= 0.1
df_sebg_nl <- df_sebg[df_sebg$p.adjust <= 0.1, ]
# sort by p.adj
df_sebg_nl <- df_sebg_nl[order(df_sebg_nl$p.adjust),]

# 2.2 Load .csv files of 13th Task
input.dir <- '/Users/christina.hillig/R_studio/NGS_Pathway_analysis/output/GSEA_output/2021-11-02/Sebaceous glands_output/ST_cdr_patient_annotation_condition'

file.name <- 'condition 1Reference_positiveLog2FC_REACTOME_Pathway_Enrichment_Analysis_.xlsx'

comparison.AD <- xlsx::read.xlsx(file.path(
  input.dir, '13_DGE_Approaches_Signature_AE_modifiedPS_20210507/L, SEB G__condition 1_vs_condition 2/ST_condition 1_vs_condition 2_glmGamPoi_DGE.xlsx',
  file.name), sheetIndex = 1, header = TRUE, Comparison='SEB G, L, AE')
# keep only enriched pathways p.adjust <= 0.1
comparison.a <- comparison.AD[comparison.AD$p.adjust <= 0.1, ]
# sort by p.adj
comparison.a <- comparison.a[order(comparison.a$p.adjust),]

comparison.PsO <- xlsx::read.xlsx(file.path(
  input.dir, 
  '13_DGE_Approaches_Signature_PSO_modifiedPS_20210507/L, SEB G__condition 1_vs_condition 2/ST_condition 1_vs_condition 2_glmGamPoi_DGE_all_genes.xlsx', 
  file.name), sheetIndex = 1, header = TRUE,  Comparison='SEB G, L, Pso')
# keep only enriched pathways p.adjust <= 0.1
comparison.b <- comparison.PsO[comparison.PsO$p.adjust <= 0.1, ]
# sort by p.adj
comparison.b <- comparison.b[order(comparison.b$p.adjust),]

# print(head(df_sebg_nl))
print(dim(df_sebg_nl))
# print(head(comparison.a))
print(dim(comparison.a))
# print(head(comparison.b))
print(dim(comparison.b))

# Way 1: Find common pathways
# intersection.pa <- intersect(comparison.a$Description, comparison.b$Description)
# ind.a <- match(comparison.a$Description, intersection.pa)
# ind.a <- ind.a[!is.na(ind.a)]
# comparison.a <- comparison.a[ind.a, ]
# print(dim(comparison.a))
# ind.b <- match(comparison.b$Description, intersection.pa)
# ind.b <- ind.b[!is.na(ind.b)]
# comparison.b <- comparison.b[ind.b, ]
# print(dim(comparison.b))

# Way 2: Select top 15 per comparison  <- discussed with Peter
sub_comparison.a <- comparison.a[1:15, ]
sub_comparison.b <- comparison.b[1:15, ]

# create combined dataframe join by `Comparison` column
# Way 1
# joined.df <- rbind(comparison.a, comparison.b, by='Description')
# Way 2
joined.df <- rbind(sub_comparison.a, sub_comparison.b, by='Description')
print(head(joined.df))
print(dim(joined.df))
# remove last row 
joined.df <- joined.df[-c(dim(joined.df)[1]), ]

# Select pathways from joined.df and read out same pathways in df_sebg_nl
ind.dfsebgnl <- match(joined.df$Description, df_sebg_nl$Description)
sub_df_sebg_nl <- df_sebg_nl[ind.dfsebgnl, ]
# Merge with joined.df
joined.df <- rbind(sub_df_sebg_nl, joined.df, by='Description')
# print(head(joined.df))
print(dim(joined.df))
# remove last row 
joined.df <- joined.df[-c(dim(joined.df)[1]), ]

# convert to numeric and reorder by p.adj value
joined.df$p.adjust <- as.numeric(joined.df$p.adjust)
joined.df <- joined.df[order(joined.df$p.adjust),]


# Plot top 15
# way 1
# sub_joined.df <- joined.df[1:15, ]
# sub_joined.df$GeneRatio <- as.numeric(sapply(sub_joined.df$GeneRatio, function(x) eval(parse(text=x))))
# way 2
joined.df$GeneRatio <- as.numeric(sapply(joined.df$GeneRatio, function(x) eval(parse(text=x))))

# create interval for count
joined.df$Count <- as.numeric(joined.df$Count)
breaks <- seq(min(joined.df$Count), max(joined.df$Count), 10)
joined.df$Counts <- cut(joined.df$Count, c(breaks, Inf), labels = breaks, include.lowest = TRUE)


# replace AE with AD and Pso with PsO
joined.df <- joined.df %>% mutate(Comparison = str_replace(Comparison, "AE", "AD"))
joined.df <- joined.df %>% mutate(Comparison = str_replace(Comparison, "Pso", "PSO"))

# https://uc-r.github.io/cleveland-dot-plots
# create grouped Barplot
png(
  filename = file.path(
    save_dir, 
    paste('Barplot', "20210830_Tasksforfigures_NL, SEB G vs NL, rest of skin_13_DGE_Approaches_Signature_AEorPso_modifiedPS_20210507.png", sep = '__')),
  width = 16, height = 8, units = "in", bg = "white", res = 300)

p1 <- ggplot(joined.df, aes(Description, GeneRatio, fill = p.adjust)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + facet_wrap(~ Comparison) + theme(text = element_text(size=14)) +
  scale_fill_gradient(low="red", high="blue")
print(p1)
dev.off()


# create grouped dotplot 
png(
  filename = file.path(
    save_dir, 
    paste('Dotplot', "20210830_Tasksforfigures_NL, SEB G vs NL, rest of skin_13_DGE_Approaches_Signature_AEorPso_modifiedPS_20210507.png", sep = '__')),
  width = 16, height = 8, units = "in", bg = "white", res = 300)
p2 <- ggplot(joined.df, aes(Description, GeneRatio)) +
  geom_point(aes(size=Counts, color = p.adjust)) + facet_wrap(~ Comparison) +
  coord_flip()  + theme(text = element_text(size=16)) +
  scale_color_gradient(low="red", high="blue")
print(p2)
dev.off()


df.parameters <- data.frame(
  Parameters=c("DGE comparison NL", "DGE comparison L, SEB, AE", "DGE comparision L, SEB, Pso",
               'FDR_cut_PA', 'Multitesting_method', 'database', '',
               'Preselected PAs', 'File_preselected_PAs', 'sheetName'),
  Values=c(file.path(input.dir, 
                     '20210830_Tasksforfigures', 
                     'NL, SEB G vs NL, rest of skin__condition 1_vs_condition 2', 
                     'ST_condition 1_vs_condition 2_glmGamPoi_DGE_all_genes.csv', file.name), 
           file.path(input.dir, '13_DGE_Approaches_Signature_AE_modifiedPS_20210507', 
                     'L, SEB G__condition 1_vs_condition 2', 
                     'ST_condition 1_vs_condition 2_glmGamPoi_DGE.xlsx', file.name), 
           file.path(input.dir, 
                     '13_DGE_Approaches_Signature_PSO_modifiedPS_20210507', 
                     'L, SEB G__condition 1_vs_condition 2', 
                     'ST_condition 1_vs_condition 2_glmGamPoi_DGE_all_genes.xlsx', file.name), 
           0.1, 'BH', 'Reactome', '', 'no', '', 'CH_Reactome'))

# Save df to excel file
write.xlsx(joined.df, file=file.path(
  save_dir, "20210830_Tasksforfigures_NL, SEB G vs NL, rest of skin_13_DGE_Approaches_Signature_AEorPso_modifiedPS_20210507.xlsx"), 
  sheetName="Plot data")
write.xlsx(comparison.AD, file=file.path(
  save_dir, "20210830_Tasksforfigures_NL, SEB G vs NL, rest of skin_13_DGE_Approaches_Signature_AEorPso_modifiedPS_20210507.xlsx"), 
  sheetName="SEB G, L, AE", append=TRUE)
write.xlsx(comparison.PsO, file=file.path(
  save_dir, "20210830_Tasksforfigures_NL, SEB G vs NL, rest of skin_13_DGE_Approaches_Signature_AEorPso_modifiedPS_20210507.xlsx"), 
  sheetName="SEB G, L, Pso", append=TRUE)
write.xlsx(df_sebg, file=file.path(
  save_dir, "20210830_Tasksforfigures_NL, SEB G vs NL, rest of skin_13_DGE_Approaches_Signature_AEorPso_modifiedPS_20210507.xlsx"), 
  sheetName="NL, SEB G", append=TRUE)
write.xlsx(df.parameters, file=file.path(
  save_dir, "20210830_Tasksforfigures_NL, SEB G vs NL, rest of skin_13_DGE_Approaches_Signature_AEorPso_modifiedPS_20210507.xlsx"), 
  sheetName="info", append=TRUE)

# Scale size of points
# + scale_size_continuous(range = c(3, 7))



