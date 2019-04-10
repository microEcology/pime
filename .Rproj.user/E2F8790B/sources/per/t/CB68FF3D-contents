#' Prevalence Interval for Microbiome Evaluation - PIME
#'
#'    This package removes the within group variation found in metataxonomic surveys (16S rRNA datasets)
#' by capturing only biological differences at high sample prevalence levels. Prevalence is defined here
#' as the percent of samples within a treatment that all contain the same taxa irrespective of relative
#' abundance. It takes a phyloseq object as input, builds hundreds of decision trees using a supervised
#' machine learning algorithm and combines them into a single model to predict the likelihood of detecting
#' any factor as source of sample variation. Higher OOB error indicates the dataset has a high relative
#' abundance of taxa with low prevalence, which is defined as noise in PIME analysis.
#'     To remove the noise, PIME applies the following steps: i) per treatment slices of the full dataset are obtained and filtered
#' using prevalence values for taxa at 5% intervals; ii) calculate Bray-Curtis (or other dissimilarity measures
#' available in vegdist function) distance matrices for each prevalence interval and compute the variance
#' partitioning using the permutational multivariate analysis of variance (PERMANOVA); iii) build decision trees
#' and determine the decline in noise in each prevalence interval. The number of taxa and the number of remaining
#' sequences for each prevalence interval are also computed.
#'
#' @name pime
#' @docType package
NULL
