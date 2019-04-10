#' Bacterial composition of men’s and women’s restrooms
#'
#' Partial dataset analyzed on the article "Differences in bacterial composition between men’s and women’s restrooms
#' and other common areas within a public building" 'DOI: 10.1007/s10482-017-0976-6'. Contains microbial composition
#' of ceiling duct samples from man's and women's restrooms in a public building from Brazil.
#'
#'
#' @format A phyloseq object with 18 samples and 3253 OTUs. Contains otu_table and sample_data slots. Sample data has 5 variables
#' SampleID BarcodeSequence   LinkerPrimerSequence  Floor Description Environment
#' \describe{
#'   \item{SampleID}{IDs of samples}
#'   \item{BarcodeSequence}{Sequences of barcodes used for sequencing}
#'   \item{LinkerPrimerSequence}{Sequence of Linker Primer used for sequencing}
#'   \item{Floor}{Which floor samples were taken}
#'   \item{Description}{Sample grouped by floor and man's or women's restroom}
#'   \item{Environment}{Sample grouped by man's or women's restroom}
#'   ...
#' }
#'
"restroom"
