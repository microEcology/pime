#'Best Prevalence
#'
#'      This function is the core of PIME. It uses \code{\link{randomForest}} to build random forests trees for samples
#'classification and variable importance computation. It performs classifications for each prevalence interval
#'returned by \code{\link{pime.prevalence}}. Variable importance is calculated, returning the Mean Decrease Accuracy (MDA),
#'Mean Decrease Impurity (MDI), overall and by sample group, and taxonomy for each feature/OTU/ASV. 
#'PIME keeps the top 30 variables with highest MDA each prevalence level.
#'
#'@param prev.list List phyloseq objects with the calculated prevalences for each interval. The output of
#' \code{\link{pime.prevalence}}
#'@param variable Any variable in the metadata to be analyzed
#'@keywords prevalence OOB classification
#'@examples
#' #Spliting data by "Environment"
#' phylist=pime.split.by.variable(restroom, "Environment")
#'
#' #Computes prevalences for each treatment "Environment" separately
#' prev=pime.prevalence(phylist)
#'
#' #Finds best prevalence interval
#' set.seed(42)
#' pime.best.prevalence(prev, "Environment")
#'
#'@return The function returns a list with a 'OOB error' table with results from random forests classification, numbers of
#'sequence and OTUs/ASVs for each prevalence level. It also returns a list 'Importance', with a table for each prevalence
#'interval, containing the OTUs/ASVs with importance (MDA and MDI) above the overall mean, with individual importance values and the
#'taxonomy classification.
#'
#'@seealso \code{\link{pime.prevalence}}
#'  \code{\link{randomForest}}
#'
#'@importFrom phyloseq "otu_table"
#'@importFrom phyloseq "sample_data"
#'@importFrom phyloseq "tax_table"
#'@importFrom rlang .data
#'@export
pime.best.prevalence<-function (prev.list, variable) {
  randon<-list()
  imp<-list()
  gs <- as(object = sample_data(prev.list[[1]]), Class = "data.frame")
  Variable <- as.factor(gs[, variable])
  pb <- progress::progress_bar$new(total = length(prev.list), clear=F, show_after = 0)
  pb$tick(0)
  for (i in prev.list){
    if ((phyloseq::taxa_are_rows(i)==TRUE)==TRUE){
      train<-t(otu_table(i))
    } else {
      train<-otu_table(i)}
    pb$tick()
    Sys.sleep(1 / 100)
    response <- Variable
    training.set <- data.frame(response, train)# Combine them into 1 data frame
    train.model <- randomForest::randomForest(response~., data = training.set, importance=TRUE)
    randon[[length(randon)+1]] <-  round(train.model$err.rate[500,1], digits=4)*100
    Importance <- train.model$importance
    imp.otu <- data.frame(rownames(Importance), Importance) 
    k=i %>% tax_table() %>% as("matrix") %>% data.frame() %>% 
      data.frame(imp.otu,.) %>% dplyr::arrange(., dplyr::desc(.data$MeanDecreaseAccuracy)) %>%
      dplyr::top_n(30, .data$MeanDecreaseAccuracy)
    names(k)[1]<-c("SequenceID")
    rownames(k)<-NULL
    imp[[length(imp)+1]]=k
   }
  #names tables from Lista as the names of the tables inside list.core
  names(randon) <- paste("Prevalence", names(prev.list))
  names(imp) <- paste("Prevalence", names(prev.list))
  #gets only the first line, all columns of every table inside perm
  Interval= paste(names(randon), "%", sep = "")
  "OOB error rate (%)" <- as.numeric(sapply(randon, cbind))
  Nseqs=sapply(prev.list, function(z) sum(phyloseq::sample_sums(z)))
  OTUs=sapply(prev.list, phyloseq::ntaxa)
  OOB.err=as.data.frame(cbind(Interval,`OOB error rate (%)`,OTUs,Nseqs))
  rownames(OOB.err) <- NULL
  print(OOB.err, row.names=FALSE)
  return(list("OOB error"=OOB.err, "Importance"=imp))
}
