# Preparation

#' Get input directory
#' 
#' @param input.folder
#' @param date.file
#' @param dataset.type
#' @param seq.technique
#' @param comparison
#' @param genename
#' @param design.function
get_inputdir <- function(input.folder, date.file, dataset.type, comparison, design.function) 
{
  input.dir <- file.path(input.folder, dataset.type, 'DEGs', date.file, 
                         paste0(comparison, '__', design.function))
  return(input.dir)
}


#' Create output directory
#' 
#' @param output.folder
#' @param dataset.type
#' @param func_task
#' @param genename
#' @return output path
get_outputdir <- function(output.folder, dataset.type, func_task, func_comparison) 
{
  # output paths
  output_dir <- file.path(output.folder, dataset.type, func_task, func_comparison, Sys.Date())
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE,)
  
  return(output_dir)
}

