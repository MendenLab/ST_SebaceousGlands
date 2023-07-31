rm(list = ls(all=TRUE))
# Set seed before loading packages 
set.seed(1)

# Load packages
library('DESeq2')
library(dplyr)
library(ggplot2)
library(xlsx)

library(ggsignif)
library(ggrepel)

#' @description Creates a Volcano plot labeling significantly expressed genes and with the option to instead label custom genes
#' 
#' @param df.data : Dataframe with columns log2FoldChange, pvalue, padj
#' @param log2fc.cut : Threshold for log2FC
#' @param add.labels : (optional) Names of predefined gene names
def.get.volcanoplot <- function(df.data, log2fc.cut, add.labels=NULL) 
{
  # 1.Mark diff expressed genes
  df.data <- df.data[order(df.data$padj), ]
  # add a column of NAs
  df.data$diffexpressed <- "Not Sig"
  # if log2Foldchange > 1 and padj < 0.05, set as "Up" 
  df.data$diffexpressed[df.data$log2FoldChange > log2fc.cut & df.data$padj < 0.1] <- "Up"
  # if log2Foldchange < -1 and padj < 0.05, set as "Down"
  df.data$diffexpressed[df.data$log2FoldChange < -log2fc.cut & df.data$padj < 0.1] <- "Down"
  df.data$diffexpressed <- as.factor(df.data$diffexpressed)
  df.data$diffexpressed <- factor(df.data$diffexpressed , levels = c("Not Sig", "Up", "Down"))
  
  # 2. Add labels to genes
  df.data <- def.getlabels(func_df.data=df.data, add.labels=NULL) 
  
  # 2.1 DOWN: split into two columns if label length > 20
  mask.down <- !is.na(df.data$label) & df.data$diffexpressed == "Down"
  # mask.down <- mask.down | df.data$label %in% add.labels
  ind.temp.down <- rownames(df.data[mask.down, ])
  second.col.down <- df.data$label[mask.down]
  second.half.down <- ceiling(length(second.col.down) / 2)
  if (length(second.col.down) %% 2 == 0 & length(second.col.down) > 1) 
  {
    second.half.down <- second.half.down + 1
  }
  first.half.down <- length(second.col.down) - second.half.down
  
  # 2.2 UP: split into two columns if label length > 20
  mask.up <- !is.na(df.data$label) & df.data$diffexpressed == "Up"
  # mask.up <- mask.up | df.data$label %in% add.labels 
  ind.temp.up <-  rownames(df.data[mask.up, ])
  second.col.up <- df.data$label[mask.up]
  second.half.up <- ceiling(length(second.col.up) / 2)
  first.half.up <- length(second.col.up) - second.half.up
  if (length(second.col.up) %% 2 == 0 & length(second.col.up) > 1) 
  {
    second.half.up <- second.half.up + 1
  }
  
  
  # 3. set parameters
  # 3.1 find xaxis max limit
  xmax = max(abs(df.data$log2FoldChange[
    !is.na(df.data$log2FoldChange) & (!is.infinite(df.data$log2FoldChange))]))
  xlim_magnitude = 0
  # 3.2 Get second largest y-value
  y = -log10(df.data$padj)
  ymax = max(y)
  ymax_second = sort(y)[length(y) - 1]  
  ylim_magnitude = 0
  # 3.3 Point and label size
  size.annot <- 2.5
  point.size <- 1
  fontsize <- 1
  segment.size <- 0.1
  box.padding <- 0.1
  
  # 4. Create Volcano plot
  func.volcanoplot <- ggplot2::ggplot(
    data=df.data %>% plyr::arrange(diffexpressed), 
    ggplot2::aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed, label=label)) +
    ggplot2::geom_point(size=point.size) + 
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(
      limits = c(-xmax - 8.5 * 10**xlim_magnitude, xmax + 8.5 * 10**xlim_magnitude)) 
  
  # Split into 20 each; use different nudges
  # Anotate top diff expressed genes
  # Donw regulated genes
  func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
    data=df.data[ind.temp.down[1:first.half.down], ], 
    fill = "white",  color = 'darkblue',
    xlim = c(-xmax + 40, -xmax + 50),
    ylim = c(min(-log10(df.data$padj)), max(-log10(df.data$padj))),
    size=size.annot,  max.overlaps = Inf, 
    direction    = "y",
    hjust        = 0,
    segment.size = segment.size,
    size = fontsize, aes(fontsize=fontsize, size=fontsize),
    box.padding  = box.padding, cex.lab=1.5, show.legend  = FALSE)
  if (length(second.col.down) > 1) 
  {
    func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
      data=df.data[ind.temp.down[second.half.down:length(second.col.down)], ], 
      fill = "white", color = 'darkblue',
      xlim = c(-xmax + 70, -xmax + 80),
      ylim = c(min(-log10(df.data$padj)), max(-log10(df.data$padj))),
      size=size.annot,  max.overlaps = Inf, 
      direction    = "y",
      hjust        = 0,
      segment.size = segment.size,
      size = fontsize, aes(fontsize=fontsize, size=fontsize), 
      box.padding  = box.padding, cex.lab=1.5, show.legend  = FALSE)
  } 
  
  # Up regulated genes
  func.volcanoplot <- func.volcanoplot + ggrepel::geom_label_repel(
    data=df.data[ind.temp.up[1:first.half.up], ], 
    fill = "white",  color = 'darkred',
    xlim = c(xmax - 40,
             xmax - 50),
    size=size.annot,  max.overlaps = Inf,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.1,
    size = fontsize, aes(fontsize=fontsize, size=fontsize),
    box.padding  = box.padding, cex.lab=1.5, show.legend  = FALSE)
  if (length(second.col.up) > 1) 
  {
    func.volcanoplot <- func.volcanoplot +  ggrepel::geom_label_repel(
      data=df.data[ind.temp.up[second.half.up:length(second.col.up)], ], 
      fill = "white",  color = 'darkred',
      xlim = c(xmax - 70,
               xmax - 80),
      size=size.annot,  max.overlaps = Inf,
      direction    = "y",
      hjust        = 0,
      segment.size = segment.size,
      size = fontsize, aes(fontsize=fontsize, size=fontsize),
      box.padding  = box.padding, show.legend  = FALSE)
  }
  
  # Legend 
  func.volcanoplot <- func.volcanoplot + 
    ggplot2::labs(color = "Genes", x = bquote(log[2] ~ "FC"), y=bquote(log[10] ~ "p-adj value")) +
    ggplot2::scale_color_manual(
      values=c('Not Sig' = "grey", 'Down' = "darkblue", "Up" = "darkred"), 
      labels = c("Not Sig.", "Up", "Down")) +
    ggplot2::theme(axis.title = element_text(size = 15),axis.text = element_text(size = 12),
                   legend.text = element_text(size = 10), legend.title = element_text(size = 12)) 
  
  return(func.volcanoplot)
}


def.getlabels <- function(func_df.data, add.labels) 
{
  if (!is.null(add.labels)) 
  {
    func_df.data$label <- NA
    # Add custom labels -> Todo check in CEBPB
    ind.labels <- match(add.labels, func_df.data$hgnc_symbol)
    # remove NAs
    ind.labels <- ind.labels[!is.na(ind.labels)]
    func_df.data$label[ind.labels] <- func_df.data$hgnc_symbol[ind.labels]
    
    # add top 10 DEGs annotations
    top_genes <- 10
    mask.sig <- func_df.data$diffexpressed == "Up"
    max_ele <- top_genes + 1  # sum(mask.sig, na.rm = TRUE)
    func_df.data$label[mask.sig][1:max_ele] <- func_df.data$hgnc_symbol[mask.sig][1:max_ele]
    
    mask.sig <- func_df.data$diffexpressed == "Down"
    func_df.data$label[mask.sig][1:max_ele] <- func_df.data$hgnc_symbol[mask.sig][1:max_ele]
    
  } else 
  {
    # Create a new column "label", contains the name of DEx genes (NA in case they are not)
    func_df.data$label <- NA
    # mask.sig <- func_df.data$diffexpressed != "Not Sig"
    # func_df.data$label[mask.sig] <- func_df.data$hgnc_symbol[mask.sig]
    
    # add top 10 DEGs annotations
    top_genes <- 10
    mask.sig <- func_df.data$diffexpressed == "Up"
    max_ele <- top_genes + 1  # sum(mask.sig, na.rm = TRUE)
    func_df.data$label[mask.sig][1:max_ele] <- func_df.data$hgnc_symbol[mask.sig][1:max_ele]
    
    mask.sig <- func_df.data$diffexpressed == "Down"
    func_df.data$label[mask.sig][1:max_ele] <- func_df.data$hgnc_symbol[mask.sig][1:max_ele]
  }
  
  # sort by label 
  func_df.data <- func_df.data[order(func_df.data$label), ]
  return(func_df.data)
}

output.dir <- '/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output'
output.dir <- file.path(output.dir, 'DGE_analysis', 'NL_SG_vs_NL_RestofSkin', Sys.Date())
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# Comparison AD, L SEB vs Rest of skin
df = xlsx::read.xlsx('/Volumes/CH__data/Projects/Eyerich_AG_projects/ST_Sebaceous_glands__Peter_Seiringer/output/switch_signs_log2fc/2021-09-07/switch_signs/ST_cdr_patient_annotation_condition/20210830_Tasksforfigures/NL, SEB G vs NL, rest of skin__condition 1_vs_condition 2/ST_condition 1_vs_condition 2_glmGamPoi_DGE_all_genes.xlsx',
                     sheetIndex = 1)
df$NA. <- NULL

df <- df %>% rename("log2FoldChange" = "log2fc")
df <- df %>% rename("hgnc_symbol" = "gene_symbol")

# Plot Volcano plot
pdf(file = file.path(output.dir, paste0('NL_SG_vs_NL_RestofSkin', '_Volcanoplot', '.pdf')), 
    width = 6, height = 6)
p1 <- def.get.volcanoplot(df.data=df, log2fc.cut=1, add.labels=NULL) 
print(p1)
dev.off()

