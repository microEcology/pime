#'Prevalences intervals
#'
#'
#'      This function takes the output of  \code{\link{pime.split.by.variable}}. It estimates the highest prevalence
#'possible (no empty OTU/ASV table), calculates prevalence for taxa, starting at 5% interval by increments of 5, until
#'maximum prevalence possible (no empty OTU table or dropping samples).
#'After prevalence calculation, each prevalence interval are merged.
#'
#'@param phylist Input, takes a list of phyloseq objects (one for each variable),
#'output from \code{\link{pime.split.by.variable}}
#'
#'@return A list with a phyloseq object (OTU/ASV and metadata table) per prevalence interval.
#'@keywords Prevalence
#'@examples phylist=pime.split.by.variable(restroom, "Environment")
#' prev=pime.prevalence(phylist)
#'@export
pime.prevalence = function(phylist) {
  cc=lapply(phylist, microbiome::prevalence)
  vv=sort(sapply(cc, max))
  max_cut=as.integer((vv[1]*100)-0.001)
  merged=list()
  var.phy=list()
  cutoff= as.factor(seq(5,max_cut, by=5))
  pb <- progress::progress_bar$new(total = length(cutoff), clear=F, show_after = 0)
  pb$tick(0)
  for (i in cutoff){
    for (j in names(phylist)){
      var.phy[[j]]=microbiome::core(phylist[[j]], detection = 0, prevalence= as.numeric(i)/100)
      prev.cut=do.call(phyloseq::merge_phyloseq, var.phy)
      if (any(phyloseq::sample_sums(prev.cut)==0)==FALSE){
        merged[[i]]=prev.cut
      }
    }
    pb$tick()
    Sys.sleep(1 / 100)
  }
  return(merged)
}
