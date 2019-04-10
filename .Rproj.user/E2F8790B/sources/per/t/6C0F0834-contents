#'Error prediction
#'
#'    For each prevalence interval it randomizes the samples labels into arbitrary groupings using n random
#' permutations (user defined). For each randomized, prevalence filtered dataset, the OOB error rate is
#' calculated to estimate whether the original differences in groups of samples occur by chance.
#' Results are in a list containing a table and a boxplot summarizing the results.
#'
#'
#'@param physeq Input in phyloseq format, original data, pime.prevalence() unfiltered
#'@param variable Variable to run the model, to be randomized
#'@param bootstrap Number to run randomizations
#'@param max.prev Max prevalence reached with pime.prevalence()
#'@param parallel Whether or not to run in parallel. Default is TRUE
#'@param ... Aditional parameters passed to ranger::ranger() on random forest classification.
#'
#'@examples phylist=pime.split.by.variable(restroom, "Environment")
#'prev=pime.prevalence(phylist)
#'pime.best.prevalence(prev, "Environment", method.dist="bray")
#'set.seed(123)
#'result=pime.error.prediction(restroom, "Environment", bootstrap=10, max.prev=90, parallel=TRUE)
#'result$Plot
#'result$'Results table'
#' @importFrom phyloseq "otu_table"
#' @importFrom phyloseq "sample_data"
#'
#'@export
pime.error.prediction = function (physeq, variable, bootstrap, max.prev, parallel=TRUE,...){
  cutoff = as.factor(seq(5, max.prev, by = 5))
  physeq2 = physeq
  if (parallel==TRUE) {
  cor=parallel::detectCores()
  cl <- parallel::makeCluster(cor-1)
  doParallel::registerDoParallel(cl)
  foreach::getDoParWorkers()
  results<- foreach::foreach (c = cutoff, .packages = c('phyloseq'), .combine = cbind) %do% {
    randon=NULL
    tab=data.frame(replicate(bootstrap,sample(sample_data(physeq)[[variable]],replace=TRUE)))
    rownames(tab) = phyloseq::sample_names(physeq)
    phyloseq::sample_data(physeq2) = tab
    df90.2 = as(sample_data(physeq2), "data.frame")
    foreach::foreach (g = colnames(df90.2),
                      .packages = c('phyloseq', 'ranger'), .combine = c) %dopar% {
                        split = list()
                        split = pime::pime.split.by.variable(physeq2, g)
                        varphy=lapply(split, function(y) microbiome::core(y,detection=0, prevalence=as.numeric(c)/100))
                        prev.cut = do.call(phyloseq::merge_phyloseq, varphy)
                        if (any(phyloseq::sample_sums(prev.cut) ==0) ==FALSE) {
                          merged = prev.cut
                        }
                        if ((phyloseq::taxa_are_rows(merged) == TRUE) == TRUE) {
                          otu = t(otu_table(merged))
                        } else {
                          otu = otu_table(merged) }
                        response = as.factor(df90.2[, g])
                        training.set <- data.frame(response, otu)
                        train.model = ranger::ranger(response ~ .,data = training.set,...)
                        randon = train.model$prediction.error
                      }
  }
  ### Stop cluster
  parallel::stopCluster(cl)
  } else {
    results<- foreach::foreach (c = cutoff, .packages = c('phyloseq'), .combine = cbind) %do% {
      randon=NULL
      tab=data.frame(replicate(bootstrap,sample(sample_data(physeq)[[variable]],replace=TRUE)))
      rownames(tab) = phyloseq::sample_names(physeq)
      sample_data(physeq2) = tab
      df90.2 = as(sample_data(physeq2), "data.frame")
      foreach::foreach (g = colnames(df90.2),
                        .packages = c('phyloseq', 'ranger'), .combine = c) %do% {
                          split = list()
                          split = pime::pime.split.by.variable(physeq2, g)
                          varphy=lapply(split, function(y) microbiome::core(y,detection=0, prevalence=as.numeric(c)/100))
                          prev.cut = do.call(phyloseq::merge_phyloseq, varphy)
                          if (any(phyloseq::sample_sums(prev.cut) ==0) ==FALSE) {
                            merged = prev.cut
                          }
                          if ((phyloseq::taxa_are_rows(merged) == TRUE) == TRUE) {
                            otu = t(otu_table(merged))
                          } else {
                            otu = otu_table(merged) }
                          response = as.factor(df90.2[, g])
                          training.set <- data.frame(response, otu)
                          train.model = ranger::ranger(response ~ .,data = training.set,...)
                          randon = train.model$prediction.error
                        }
    }
  }
  colnames(results)=paste(levels(c))
  res.plot=as.data.frame(results) %>% tidyr::gather("Prevalence", "OOBerror", 1:ncol(.)) %>%
    ggplot2::ggplot(.)+ggplot2::aes(x=forcats::fct_inorder(Prevalence), y=OOBerror)+
    ggplot2::geom_point(size=1)+ggplot2::geom_boxplot()+
    ggplot2::ylim(0,1)+ggplot2::ylab("OOB error")+ggplot2::xlab("Prevalences (%)")
  colnames(results)=paste("Prevalence",levels(c), "%", sep = "")
  return(list("Results table"=as.data.frame(results), "Plot"=res.plot))
}
