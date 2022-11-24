## now  then match to gene names in stages
unimachR <-function(ids,hgnc.table){
  ## function that takes a character vector of uniprot ids (ids=ids)
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
  
  
  cat("unmachR on",length(ids),"Uniprot IDs\n")
  
  
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
    filters = 'uniprotswissprot', 
    values = ids, 
    mart = hsmart
  )
  
  ## correct the returned names with latest hgnc
  mapping$hgnc_symbol =
    HGNChelper::checkGeneSymbols(mapping$hgnc_symbol,map=hgnc.table) %>%
    dplyr::select(Suggested.Symbol) %>%
    makeChar
  
  
  mapping %<>% 
    dplyr::rename(uniprot_id = uniprotswissprot)
  
  ## get the missed ids, and those with an NA gene name
  na_genes =
    mapping %>%
    filter(is.na(hgnc_symbol)) %>%
    dplyr::select(uniprot_id) %>%
    makeChar()
  
  missed_ids = c(ids[ids %notin% mapping$uniprot_id],na_genes)
  
  ## keep just the IDs with a valid gene name from hgnc
  mapping %<>%
    filter(!is.na(hgnc_symbol))
  
  ## quantify how many ids were matched vs not matched
  matched = length(ids[ids %in% mapping$uniprot_id])
  no_matched = length(ids[ids %notin% mapping$uniprot_id])
  total = length(ids)
  perc_matched = round(matched*100/total)
  perc_nomatched = round(no_matched*100/total)
  
  cat("\nBiomart found",matched,"ids (",perc_matched,"%)\n")
  cat("\nCould not map ",no_matched,"ids (",perc_nomatched,"%)\n")
  
  
  
  
  if(no_matched==0){
    
    mapping_out = mapping
    
  } else if (no_matched>0){
    
    
    ##### matching stage 2 #####
    ## scraping uniprot webserver
    cat("\n\n\n##### matching stage 2 #####\n")

    
    ## have XX uniprot IDs not mapped to genes, attemp to look these up on uniprot web and recover gene name
    missed_ids = ids[ids %notin% mapping$uniprot_id]
    
    # make the ids table into a data table to allow quick join to the deleted IDs
    missed_ids %<>% data.table()
    
    colnames(missed_ids) = "uniprot_id"

    ids_toScrape = makeChar(missed_ids$uniprot_id)
    cat("scraping uniprot server searching for",length(ids_toScrape),"ids\n")
    
    
    id = ids_toScrape[1]
    mapped_ids = tibble()
    del_ids = tibble()
    for(i in 1:length(ids_toScrape)){
      
      # reset temp to null
      temp=NULL
      
      id = ids_toScrape[i]
      
      # report progress
      a = signif(i*100/length(ids_toScrape),2)
      cat(id,a,"% \n")
      ## error here is likely caused by deleted ID
      acc_url = paste0("https://www.uniprot.org/uniprot/",id,".fasta")
      
      # scrape uniprot web server
      temp = fread(acc_url)
      
      ## only parse and add to table if not empty, some deleted IDs seem to be escaping
      ## or if not using the most up to date deleted IDs table
      
      if(length(temp)>0){
        
        temp = str_split(colnames(temp), "GN=")
        
        temp2 = str_split(unlist(temp)[2], " ")
        
        name = unlist(temp2)[1]
        
        temp =
          tibble(name=name,
                 id=id)
        
        mapped_ids %<>% rbind.data.frame(temp)
      } else {
        temp_df =
          tibble(
            id
          )
        del_ids %<>% rbind.data.frame(temp_df)
      }
      
    }
    
    cat("\nscraping done for",nrow(mapped_ids),"ids\n")
    
    ## remove any that have no gene name on uniprot
    mapped_ids %<>%
      filter(!is.na(name))
    
    ## combine the scraped ids and biomart matched ids
    colnames(mapped_ids) = colnames(mapping)
    mapping_out = rbind.data.frame(mapping,mapped_ids)
    
  }
    
    ## report total mapping stats
    total_mapped = length(unique(mapping_out$uniprot_id))
    total_genes = length(unique(mapping_out$hgnc_symbol))
    total_in = length(ids)
    
    perc_matched = signif(total_mapped*100/total_in,4)
    
    cat("\n\n\nunimachR mapped",total_mapped,"ids out of",total_in,"(",perc_matched,"%)\n")
    cat(total_mapped,"ids mapped to",total_genes,"HGNC symbols\n")
    
    
    ## find the uniprot IDs that map to more than one gene name
    multi_ids =
      mapping_out %>%
      group_by(uniprot_id) %>%
      summarise(count=n()) %>%
      filter(count >1) %>%
      arrange(desc(count))
    
    multi_genes =
      mapping_out %>%
      group_by(hgnc_symbol) %>%
      summarise(count=n()) %>%
      filter(count >1) %>%
      arrange(desc(count))
    
    ## compose a list as output
    mapping_out_list =
      list(mapping_out,multi_ids,multi_genes,del_ids)
    
    
    return(mapping_out_list)
    
}
