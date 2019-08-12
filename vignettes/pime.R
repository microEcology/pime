## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "doc/"
)

## ---- eval=FALSE---------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install("phyloseq")

## ---- eval=FALSE---------------------------------------------------------
#  seqtab = readRDS("path_to_file/sequence_table_final.rds")
#  tax= readRDS("path_to_file/tax_final.rds")
#  map <- "path_to_file/sample_data.txt"
#  ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE),
#                 tax_table(taxa))
#  sample_metadata = import_qiime_sample_data(map)
#  physeq =merge_phyloseq(ps, sample_metadata)

## ---- eval=FALSE---------------------------------------------------------
#  jsonbiomfile = "path_to_file/otu_table_fix.biom"
#  mapfile = "path_to_file/v35_map_uniquebyPSN.txt"
#  biom = import_biom(jsonbiomfile, mapfile, parseFunction=parse_taxonomy_greengenes)
#  map = import_qiime_sample_data(mapfile)
#  input = merge_phyloseq(biom,map)

## ---- eval = FALSE-------------------------------------------------------
#  #install.packages("devtools")
#  library(devtools)
#  install_github("microEcology/pime")

## ------------------------------------------------------------------------
library(pime)
data("restroom")
pime.oob.error(restroom, "Environment")

## ------------------------------------------------------------------------
per_variable_obj= pime.split.by.variable(restroom, "Environment")
per_variable_obj

## ------------------------------------------------------------------------
prevalences=pime.prevalence(per_variable_obj)
head(prevalences)

## ------------------------------------------------------------------------
set.seed(42)
best.prev=pime.best.prevalence(prevalences, "Environment")

## ------------------------------------------------------------------------
imp65=best.prev$`Importance`$`Prevalence 65`
head(knitr::kable(imp65, format="markdown"))

#To get the table with OOB error results.
#best.prev$`OOB error`

## ------------------------------------------------------------------------
prevalence.65 = prevalences$`65`
prevalence.65

## ---- eval=FALSE---------------------------------------------------------
#  randomized=pime.error.prediction(restroom, "Environment", bootstrap = 10, parallel = TRUE, max.prev = 95)
#  randomized$Plot
#  randomized$'Table results'

## ---- eval=FALSE---------------------------------------------------------
#  replicated.oob.error= pime.oob.replicate(prevalences, "Environment", bootstrap = 10, parallel = TRUE)

