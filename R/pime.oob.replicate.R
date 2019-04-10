#' Replicated random forests
#'
#'
#'This function is part of the error detection step. It performs a random forests classification, with
#'\code{\link{ranger}}, and computes the OOB error for n replications in each prevalence
#'interval without randomizing the sample labels. Returns a results table with prediction error and a boxplot.
#'
#'@param prev.list List of phyloseq objects. Output of \code{\link{pime.prevalence}}
#'@param variable Variable to run the model
#'@param bootstrap Number to run repetitions
#'@param parallel Whether to run parallel or not. Default TRUE
#'@examples phylist=pime.split.by.variable(restroom, "Environment")
#'prevalences=pime.prevalence(phylist)
#'set.seed(42)
#'tab=pime.bootstraped.randonForest(prevalences, "Environment", bootstrap=10, parallel=TRUE)
#'tab$Plot
#'tab$'Results table'
#' @importFrom phyloseq "otu_table"
#' @importFrom phyloseq "sample_data"
#' @importFrom foreach %do% %dopar%
#'
#'@export
pime.oob.replicate = function (prev.list, variable, bootstrap, parallel=TRUE){
  if (parallel==TRUE){
    cor=parallel::detectCores()
    cl <- parallel::makeCluster(cor-1)
    doParallel::registerDoParallel(cl)
    foreach::getDoParWorkers()
    results<-foreach::foreach(v= names(prev.list), .packages = 'phyloseq', .combine = cbind)%do%{
      randon=NULL
      df90.2 = as(phyloseq::sample_data(prev.list[[v]]), "data.frame")
      foreach::foreach (g = seq_len(bootstrap), .combine = c, .packages = c('phyloseq','ranger'))%dopar% {
        if ((phyloseq::taxa_are_rows(prev.list[[v]]) == TRUE) == TRUE) {
          otu = t(otu_table(prev.list[[v]]))
        } else {
          otu = otu_table(prev.list[[v]])
        }
        response = as.factor(df90.2[, variable])
        training.set <- data.frame(response, otu)
        train.model = ranger::ranger(response ~ .,data = training.set)
        randon = train.model$prediction.error
      }
    }
    ### Stop cluster
    parallel::stopCluster(cl)
  } else {
    results<-foreach::foreach(v= names(prev.list), .packages = 'phyloseq', .combine = cbind) %do% {
      randon=NULL
      df90.2 = as(sample_data(prev.list[[v]]), "data.frame")
      foreach::foreach (g = seq_len(bootstrap), .combine = c, .packages = c('phyloseq','ranger')) %do% {
        if ((phyloseq::taxa_are_rows(prev.list[[v]]) == TRUE) == TRUE) {
          otu = t(otu_table(prev.list[[v]]))
        } else {
          otu = otu_table(prev.list[[v]])
        }
        response = as.factor(df90.2[, variable])
        training.set <- data.frame(response, otu)
        train.model = ranger::ranger(response ~ .,data = training.set)
        randon = train.model$prediction.error
      }
    }
  }
  results=as.data.frame(results)
  colnames(results)=paste(names(prev.list))
  res.plot=as.data.frame(results) %>% tidyr::gather("Prevalence", "OOBerror", 1:ncol(.)) %>%
    ggplot2::ggplot(.)+ggplot2::aes(x=forcats::fct_inorder(Prevalence), y=OOBerror)+
    ggplot2::geom_point(size=1)+ggplot2::geom_boxplot()+ggplot2::theme_bw()+
    ggplot2::ylim(0,1)+ggplot2::ylab("OOB error")+ggplot2::xlab("Prevalences (%)")
  colnames(results)=paste("Prevalence",names(prev.list), "%", sep = "")
  return(list("Results table"=as.data.frame(results), "Plot"=res.plot))
}
