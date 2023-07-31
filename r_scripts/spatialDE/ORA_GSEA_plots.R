library('ggnewscale')


def.plot.enricher_barplot <- function(enricher.result, method, showCategories, title, 
                             width, height, output.dir, func_colours=NULL) 
{
  if (typeof(showCategories) != 'character') 
  {
    if (showCategories > nrow(enricher.result)) 
    {
      showCategories = nrow(enricher.result)
    }
  } else 
  {
    # get intersection Pathways which are in enrichObject and in the user provided Pathway list
    showCategories = intersect(enricher.result$Description, showCategories)
  }
  
  
  df <- enricher.result
  df@result <- df@result[seq(1,showCategories), ]
  # Convert GeneRatio to numeric
  df@result$GeneRatio.numeric <- as.numeric(
    sapply(df@result$GeneRatio, function(x) eval(parse(text=x))))
  df@result <- df@result[order(-df@result$GeneRatio.numeric), ]
  
  
  pdf(file = file.path(output.dir, title), width = width, height = height)
 
  p1 = barplot(df, showCategory = showCategories, color='p.adjust', x='GeneRatio', 
               legend.text = TRUE)  +
    xlab("GeneRatio") + ylab('') + theme_classic() + ggtitle(method)
  if (!is.null(func_colours))
  {
    p1 <- p1 + scale_fill_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                                    guide=guide_colorbar(reverse=TRUE))
  }
  p1 <- p1 + 
    theme(legend.position = c(0.85, 0.12), axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18), plot.title = element_text(size = 20),
          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) 
    # theme(legend.background = element_rect(fill = "white", colour = "black"))
  
  # cowplot::plot_grid(p1, ncol = 1, nrow = 1)
  print(p1)
  dev.off()
}
