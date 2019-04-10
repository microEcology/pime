#'Prevalences intervals
#'
#'
#'      This function takes the output of pime.split.by.variable(). It estimates the highest prevalence
#'possible (no empty OTU table), calculates prevalence for taxa at 5% interval by increments of 5, until
#'maximum prevalence possible (no empty OTU table or droping samples).
#'After prevalence calculation, merges phyloseq objects for each interval and variable.
#'
#'@param phylist Input, takes a list of phyloseq objects (one for each variable),
#'output from pime.split.by.variable()
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
  for (i in cutoff){
    for (j in names(phylist)){
      var.phy[[j]]=microbiome::core(phylist[[j]], detection = 0, prevalence= as.numeric(i)/100)
      prev.cut=do.call(phyloseq::merge_phyloseq, var.phy)
      if (any(phyloseq::sample_sums(prev.cut)==0)==FALSE){
        merged[[i]]=prev.cut
      }
    }
  }
  return(merged)
}
