## script for processing deleted IDs table from uniprot
# downloaded October 2022 https://www.uniprot.org/help/deleted_accessions

library(tidyverse)
library(magrittr)
library(data.table)


## requires quite a bit of ram to load the data especially the huge delac_tr.txt
## worked on 
## OS Name	Microsoft Windows 10 Pro
## Processor	Intel(R) Core(TM) i5-10310U CPU @ 1.70GHz, 2208 Mhz, 4 Core(s), 8 Logical Processor(s)
## Installed Physical Memory (RAM)	16.0 GB


#################
#### OPTIONS ####
#################

options(max.print=100)



###################
#### Functions ####
###################

# source custom functions
devtools::source_url("https://github.com/d0minicO/phosphoR/blob/main/customFunctions.R?raw=TRUE")



################
#### inputs ####
################


base = "C:/Users/dowens/OneDrive/unimachR/unimachR/"

## deleted from UniProtKB/Swiss-Prot.
dat1 = 
  data.table::fread(paste0(base,"delac_sp.txt"),skip = 27,col.names = "id") %>%
  dplyr::slice(-(1787:1792))


## deleted from UniProtKB/TrEMBL
dat2 = 
  data.table::fread(paste0(base,"delac_tr.txt"),skip = 27,showProgress = T,col.names = "id")


# combine
dat = rbind.data.frame(dat1,dat2)

# save
saveRDS(dat,paste0(base,"Uniprot_Deleted_accessions_all.Rds"))

## quick check to make sure all entries included
nrow(dat1)+nrow(dat2)==nrow(dat)
