#' Per variable slices
#'
#' The used input file (in phyloseq object format) is sliced according to a metadata
#' variable present in the metadata file. Code extracted from metagMisc: v.0.0.4 - DOI 10.5281/zenodo.172500
#'
#'@param physeq The input file in phyloseq object format
#'@param variable Any variable present in the metadata to be analyzed "variable to split"
#'@param drop_zeroes Remove OTUs/ASVs with zeros, default is TRUE. If FALSE, OTUs with zero counts will be kept.
#' Further analysis based on dissimilarity measurements (Bray-Curtis) might produce empty files.
#'@examples phylist=pime.split.by.variable(restroom, "Environment")
#'@importFrom phyloseq "otu_table"
#'@export
pime.split.by.variable <- function(physeq, variable, drop_zeroes = T){

  # Check the input
  if(is.null(phyloseq::sample_data(physeq, errorIfNULL = F))){
    stop("Sample data is missing in the phyloseq-object.\n")
  }

  # Extract sample meta-data
  mtd <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")

  if(!variable %in% colnames(mtd)){
    stop("Grouping variable is missing from the sample data of phyloseq-object.\n")
  }

  if(class(mtd[, variable]) %in% c("integer", "numeric") ){
    if( length( unique(mtd[, variable]) ) > 5){
      stop("Groupping variable is numeric and it has too many levels. Consider transforming it to factor.\n")
    } else {
      warning("Groupping variable is numeric and it was coerced to factor.\n")
      mtd[, variable] <- factor(mtd[, variable])
    }
  }

  if(length(table(mtd[, variable])) == 1){
    cat("Warning: there is only one group of samples in the resulting list.\n")
  }

  # Add sample IDs to the meta-data
  smp <- data.frame(
    SID = phyloseq::sample_names(physeq),
    mtd,
    stringsAsFactors = F)

  # Exatract sample names by the specified variable
  svv <- plyr::dlply(.data = smp, .variables = variable, .fun = function(z){ z$SID })

  # Extract samples by groupping variable
  res <- plyr::llply(.data = svv, .fun = function(z){ phyloseq::prune_samples(z, x = physeq) })

  # Remove taxa with zero abundance
  if(drop_zeroes == TRUE){
    res <- plyr::llply(.data = res, .fun = function(x){ phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) })
  }

  return(res)
}
