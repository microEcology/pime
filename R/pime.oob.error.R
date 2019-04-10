#' Baseline noise detection
#'
#'
#'    This function predicts the likelihood of detecting any factor as source of sample variation.
#' It runs random forests on prevalence unfiltered dataset (the original dataset). The model
#' performance is indicated by the out-of-bag (OOB) error rate. Higher OOB error
#' indicates the dataset has a high relative abundance of taxa with low prevalence, which is defined as
#' noise in PIME analysis. There is no formal criteria for definition of low or high OOB error but empirical
#' tests showed that PIME can improve microbiome differences when OOB error >= 0.01.
#'
#'@param physeq The input file in phyloseq object format
#'@param variable Any variable present in the metadata to be analyzed. "variable to run the classification"
#'@param ... Aditional parameters passed to ranger::ranger() on random forest classification.
#'@keywords Classification OOBerror
#'@examples pime.oob.error(restroom, "Environment")
#'@export
pime.oob.error= function(physeq, variable,...){
  set.seed(2125)
  if ((phyloseq::taxa_are_rows(physeq))==TRUE) {
    train=t(otu_table((physeq)))
  } else {(train=otu_table(physeq))}
  #train=otu_table(physeq)
  # Make one column for our outcome/response variable
  met_data <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")
  response <- factor(met_data[, variable])
  training.set <- data.frame(response, train)
  train.model = ranger::ranger(response ~ ., data = training.set,...)
  res=train.model$prediction.error
  if (res <=0.01) {print("OOB error rate is zero. Your dataset presents large differences.
                      Prevalence evaluation might not be necessary")
  } else {return(res)}
}
