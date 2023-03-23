hgnc_symbols=c("OLFML1","MARVELD3","MUC21","KCNK6")


## now  then match to gene names in stages
reverseMatch <-function(hgnc_symbols){
  ## function that takes a character vector of hgnc_symbols (hgnc_symbols=hgnc_symbols)
  ## and a table of hgnc compliant gene mappings (hgnc.table) # # downloaded April 2022 https://github.com/waldronlab/HGNChelper/blob/master/data/hgnc.table.rda, can download use table from dowens github
  
  ## note, the deleted uniprot IDs database is not necessary now after tweaking code.
  ## when scraping uniprot IDs, empty entries will not cause error and will just be saved into a new output df so they can be checked against the deleted IDs database later, if desired
  ## but most likely, if they are not listed on uniprot anymore, then you dont need them!
  ## this was done as loading the deleted uniprot IDs takes A LONG TIME and will sometimes die due to memory allocation failures
  
  ## and will return a list containing four df elements
  ## [[1]] a df with each uniprot ID as a row with a matched HGNC-compliant gene symbol
  ## columns are 1) hgnc_symbol 2) uniprot_id
  ## use this df for later matching uniprot IDs to gene symbols
  
  ## [[2]] a df containing all the uniprot IDs that match to more than one gene symbol
  
  ## [[3]] a df containing all the gene symbols that match to more than one uniprot ID
  
  ## NEW
  ## [[4]] a df containing all the uniprot IDs that were not found on uniprot (likely deleted /demerged IDs)
  
  
  require(tidyverse)
  require(magrittr)
  require(biomaRt)
  require(data.table)
  
  # source DO's custom functions on github
  devtools::source_url("https://github.com/d0minicO/phosphoR/blob/main/customFunctions.R?raw=TRUE")
  
  
  cat("unmachR on",length(hgnc_symbols),"Uniprot IDs\n")
  
  
  ##### matching stage 1 #####
  ## batch match using biomart
  cat("\n\n\n##### matching stage 1 #####\n")
  cat("batch match using biomart\n")
  
  
  # load the biomart object to allow retreiving uniprot IDs from ENS IDs
  cat("\nloading the mart from biomart, takes some time...\n")
  hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
  cat("Done!\n")
  
  cat("\nSubmitting batch query to biomart\n")
  # use biomart to match gene symbols ENSEMBL IDs to Uniprot
  mapping <- getBM(
    attributes =c('hgnc_symbol','uniprotswissprot'),
    filters = 'hgnc_symbol', 
    values = hgnc_symbols, 
    mart = hsmart
  )
  
  ## clean up empty rows
  mapping %<>%
    filter(uniprotswissprot!="")
  
  ## check to see if one-to-one mapping was found
  max_map =
    mapping %>%
    group_by(hgnc_symbol) %>%
    dplyr::count() %>%
    pull(n) %>%
    max()
  
  if(max_map==1){
    cat("\nOne-to-one mapping found!\n")
  } else {
    cat("\nMultiple uniprot IDs found!\n")
  }
  
  return(mapping)
  
}
