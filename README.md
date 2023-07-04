# unimachR
 R function for matching uniprot IDs to HGNC-compliant gene names


Takes a character vector of uniprot ids (ids=ids)
and a table of hgnc compliant gene mappings (hgnc.table) # # downloaded April 2022 https://github.com/waldronlab/HGNChelper/blob/master/data/hgnc.table.rda

Note, the deleted uniprot IDs database is not necessary now after tweaking code.
when scraping uniprot IDs, empty entries will not cause error and will just be saved into a new output df so they can be checked against the deleted IDs database later, if desired
but most likely, if they are not listed on uniprot anymore, then you dont need them!
this was done as loading the deleted uniprot IDs takes A LONG TIME and will sometimes die due to memory allocation failures
  
Returns a list containing four df elements
 - [[1]] a df with each uniprot ID as a row with a matched HGNC-compliant gene symbol
  columns are 1) hgnc_symbol 2) uniprot_id
  use this df for later matching uniprot IDs to gene symbols
  
 - [[2]] a df containing all the uniprot IDs that match to more than one gene symbol
  
 - [[3]] a df containing all the gene symbols that match to more than one uniprot ID
  
# NEW
 [[4]] a df containing all the uniprot IDs that were not found on uniprot (likely deleted /demerged IDs)