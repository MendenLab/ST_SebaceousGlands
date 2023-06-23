library("xlsx")
library('DOSE')
library('org.Hs.eg.db')
library('qdapTools')


rename_genetoentrezid <- function(df.dge_results) 
{
  # I. replace gene symbol by entrezid of gene
  if ("entrezid" %nin% colnames(df.dge_results))
  {
    # works approach 
    Entrez_ID <- mapIds(org.Hs.eg.db, df.dge_results$hgnc_symbol, 'ENTREZID', 'SYMBOL')
    Entrez_ID_df = qdapTools::list2df(Entrez_ID, col1 = "entrezid", col2 = "hgnc_symbol")
    Entrez_ID_df = Entrez_ID_df[!is.na(Entrez_ID_df$entrezid), ]
    
    gseaDat <- filter(df.dge_results, !is.na(Entrez_ID))
    df.dge_results <- filter(df.dge_results, !is.na(Entrez_ID))
    df.dge_results$entrezid <- Entrez_ID_df$entrezid
  }
  
  return(df.dge_results)
}


do_rank_genes <- function(df.dge_results) 
{
  # II. Rank all genes based on their fold change
  ranked_genes.genesymbol <- df.dge_results$log2FoldChange
  names(ranked_genes.genesymbol) <- row.names(df.dge_results)
  ranked_genes.genesymbol <- sort(ranked_genes.genesymbol, decreasing = T)
  
  # with entrezID as names
  ranked_genes.entrezid <- df.dge_results$log2FoldChange
  names(ranked_genes.entrezid) <- df.dge_results$entrezid
  ranked_genes.entrezid <- sort(ranked_genes.entrezid, decreasing = T)
  
  return(list(ranked_genes.genesymbol, ranked_genes.entrezid))
}


get_significantgenes <- function(df.dge_results, op) 
{
  # II. Get genes belonging to enriched patterns 
  df.sig <- df.dge_results[op(df.dge_results$SEBACEOUS.GLAND,  1), ]
  degenes.sig <- df.sig$entrezid
  degenes.sig <- as.character(na.exclude(degenes.sig))       
  
  return(list(df.sig, degenes.sig))
}


setreadable_pa <- function(paenrich_object) 
{
  if (!is.null(paenrich_object) && dim(paenrich_object)[1] > 0 && dim(paenrich_object)[2] > 0)
  {
    paenrich_object <- DOSE::setReadable(paenrich_object, 'org.Hs.eg.db', 'ENTREZID')
  }else{
    paenrich_object <- paenrich_object
  }
  return(paenrich_object)
}


save_enrichobject_as_csv <- function(
  paenrich_object, condition, comparison, pattern, pa_database, method, df.parameters, output_path) 
{
  df_enrich = paenrich_object@result 
  xlsx_data =  paste(condition, comparison, pa_database, method, '.xlsx', sep = "_")
  # csv_data = paste(condition, comparison, pa_database, method, '.csv', sep = "_")
  if (nrow(df_enrich) > 1) 
  {
    result <- tryCatch({
      # The code I want run
      write.xlsx(df_enrich, file.path(output_path, xlsx_data), sheetName=as.character(pattern), 
                 row.names = FALSE, append=TRUE)
    }, warning = function(war) {
      # Is executed if warning encountered
    }, error = function(err) {
      # Is executed if error encountered - if file does not exist run this first
      write.xlsx(df_enrich, file.path(output_path, xlsx_data), sheetName=as.character(pattern), 
                 row.names = FALSE)
    })
    
    
    result <- tryCatch({
      # The code I want run
      write.xlsx(df.parameters, file.path(output_path, xlsx_data), sheetName="Infos", append=TRUE)
      # write.csv(df_enrich, file.path(output_path, csv_data), row.names = FALSE)
    }, warning = function(war) {
      # Is executed if warning encountered
    }, error = function(err) {
      # Is executed if error encountered - if file does not exist run this first
    })

  }
}


save_df_as_csv <- function(
  df, condition, comparison, pattern, pa_database, method, df.parameters, output_path) 
{
  xlsx_data =  paste(condition, comparison, pa_database, method, '.xlsx', sep = "_")
  # csv_data = paste(condition, comparison, pa_database, method, '.csv', sep = "_")
  if (nrow(df) > 1) 
  {
    result <- tryCatch({
      # The code I want run
      write.xlsx(df, file.path(output_path, xlsx_data), sheetName=as.character(pattern), 
                 row.names = FALSE, append=TRUE)
    }, warning = function(war) {
      # Is executed if warning encountered
    }, error = function(err) {
      # Is executed if error encountered - if file does not exist run this first
      write.xlsx(df, file.path(output_path, xlsx_data), sheetName=as.character(pattern), 
                 row.names = FALSE)
    })

    result <- tryCatch({
      # The code I want run
      write.xlsx(df.parameters, file.path(output_path, xlsx_data), sheetName="Infos", append=TRUE)
      # write.csv(df, file.path(output_path, csv_data), row.names = FALSE)
    }, warning = function(war) {
      # Is executed if warning encountered
    }, error = function(err) {
      # Is executed if error encountered - if file does not exist run this first
    })
  }
}

