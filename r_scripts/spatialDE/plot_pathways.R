#! /usr/bin/Rscript

library(ggplot2)
library(enrichplot) 
library(cowplot)
library("ggnewscale")
library(UpSetR)
# library(pathfindR)
library(forcats)


########################################################################################
################################ --------> plots <---------- ###########################
# Source: https://yulab-smu.github.io/clusterProfiler-book/chapter12.html
# Gene Concept Network
fig.pathways.REACTOME <- function(reactome_res, entrezid_log2fc, showCategories, title,
                                  width, height, output.dir) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(reactome_res)) 
    {
      showCategories = nrow(reactome_res)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(reactome_res$Description, showCategories)
  }

  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 = enrichplot::cnetplot(reactome_res, showCategory = showCategories, 
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + 
    ggtitle("REACTOME: Pathway Enrichment Analysis using DB Reactome") 

  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
  
  p1 
}


#' Plot pathways in a dotplot
#' 
#' @param pathway_res : enrichobject
#' @param showCategories : str, int
#'   which or how many pathways to visualise; limited by max. number of enriched pathways
fig.pathway.dotplot <- function(pathway_res, showCategories, method, title,
                                width, height, output.dir) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(pathway_res)) 
    {
      showCategories = nrow(pathway_res)
    }
  } 
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 = enrichplot::dotplot(pathway_res, showCategory = showCategories, color='p.adjust') + 
    ggtitle(paste(method, ": Pathway Enrichment Analysis", sep = " "))
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}



fig.pathway.upsetplot <- function(pathway_res, showCategories, method, title,
                                  width, height, output.dir) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(pathway_res)) 
    {
      showCategories = nrow(pathway_res)
    }
  } 
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 = enrichplot::upsetplot(pathway_res, n = showCategories) + 
    ggtitle(paste(method, ": Pathway Enrichment Analysis", sep = " "))
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
  
  # UpSetR::upset()
}


fig.pathway.barplot <- function(pathway_res, showCategories, method, title,
                                  width, height, output.dir, func_colours=NULL) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(pathway_res)) 
    {
      showCategories = nrow(pathway_res)
    }
  } 

  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 = barplot(pathway_res, showCategory = showCategories, color='p.adjust', x='GeneRatio')  +
    xlab("GeneRatio") + ylab('Pathways') + theme_classic() + 
    ggtitle(paste(method, ": Pathway Enrichment Analysis", sep = " "))
  if (!is.null(func_colours))
  {
    p1 <- p1 + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                                    guide=guide_colorbar(reverse=TRUE)) 
  }
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}


fig.pathways.cnetplot <- function(pathway_res, method, entrezid_log2fc, showCategories, title,
                                  width, height, output.dir) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(pathway_res)) 
    {
      showCategories = nrow(pathway_res)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(pathway_res$Description, showCategories)
  }
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 = enrichplot::cnetplot(pathway_res, showCategory = showCategories, 
                            categorySize = "pvalue", foldChange = entrezid_log2fc, 
                            colorEdge = TRUE) + 
    ggtitle(paste(method, ": Pathway Enrichment Analysis", sep = " ")) 
  
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}


def.plot.heatmap <- function(enricher.obj.ref, enricher.obj.ctrl, title, width, height, 
                             output.dir, ref.group, ctrl.group) 
{
  df_enrich_all.ref <- enricher.obj.ref@result
  # reverse order by p.adjusted 
  test_paths = df_enrich_all.ref[
    order(df_enrich_all.ref[,'Description'], -df_enrich_all.ref[,'p.adjust']),]
  # remove duplicated pathways?
  # test_paths = test_paths[!duplicated(test_paths$Description), ]
  test_paths_final = test_paths$p.adjust%>%as.data.frame()
  rownames(test_paths_final) = test_paths$Description
  colnames(test_paths_final) = ref.group
  
  df_enrich_all.ctrl <- enricher.obj.ctrl@result
  # reverse order by p.adjusted 
  ctrl_paths = df_enrich_all.ctrl[
    order(df_enrich_all.ctrl[,'Description'], -df_enrich_all.ctrl[,'p.adjust']),]
  # remove duplicated pathways?
  # ctrl_paths = ctrl_paths[!duplicated(ctrl_paths$Description), ]
  ctrl_paths_final = ctrl_paths$p.adjust%>%as.data.frame()
  rownames(ctrl_paths_final) = ctrl_paths$Description
  colnames(ctrl_paths_final) = ctrl.group
  
  all_paths = merge(test_paths_final, ctrl_paths_final, by = 0, all = TRUE)
  all_paths = all_paths%>%na.omit()
  rownames(all_paths) = all_paths$Row.names
  all_paths = all_paths[-1]
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
  p1 <- pheatmap::pheatmap(
    all_paths, main = "", cluster_rows = FALSE, cluster_cols = FALSE, angle_col = 0, 
    legend = TRUE, legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "p.adjust\n"))
  cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}


# Outsourced code:
# 1. Dotplot
# ggplot(reactome.ref, showCategory = 10, aes(GeneRatio, 
#                                             forcats::fct_reorder(Description, GeneRatio))) +
#   geom_segment(aes(xend=0, yend = Description)) +
#   geom_point(aes(color=p.adjust, size = Count)) +
#   scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
#                         trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) +
#   scale_size_continuous(range=c(2, 10)) + xlab("GeneRatio") + 
#   ylab('Pathways') + ggtitle("Reactome Pathways")

# 2. Barplot
# ggplot(reactome.ref, showCategory=10, 
#        aes(GeneRatio, fct_reorder(Description, GeneRatio), fill=p.adjust)) +
#   geom_col() + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
#                                     guide=guide_colorbar(reverse=TRUE)) +
#   xlab("GeneRatio") + ylab('Pathways') + ggtitle("WikiPathways")
