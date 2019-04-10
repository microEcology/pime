
PIME: Prevalence Interval for Microbiome Evaluation
===================================================

PIME removes the within group variation found in metataxonomic surveys (16S rRNA datasets) by capturing only biological differences at high samples prevalence levels.

How to install PIME package
===========================

To install PIME first install the devtools package. Then load the library(devtools) and run install\_github using the following commands.

``` r
#install.packages("devtools")
library(devtools)
install_github("microEcology/pime")
```

PIME uses a Phyloseq object as input. A description of the phyloseq object and a tutorial on how to create this file in R using OTU tables in many different formats is detailed into the Phyloseq website <https://joey711.github.io/phyloseq/>

Step-by-step example
====================

<img src="doc/Figure1.tif" height="500px" />

The first step in PIME is to define if the microbial community presents a high relative abundance of taxa with low prevalence, which is considered as noise in PIME analysis. This is calculated by random rorests analysis. In this example we run PIME using the restroom dataset (<https://doi.org/10.1007%2Fs10482-017-0976-6>) against the metadata variable called Environment (a variable with two categories: men’s and women’s restroom).

Prediction using random forests on full dataset. Results in Out of Bag error rate. The input file was rarefied at 500 sequences for the purpose of this example (speed up the analysis). Using a rarefied dataset is recommended at this step.

``` r
library(pime)
data("restroom")
pime.oob.error(restroom, "Environment")
#> [1] 0.4444444
```

The OOB error rate &lt;=0.1, indicated the dataset present large differences, and pime might not remove much of the noise. Higher OOB error rate indicates that the next functions should be run to find the best prevalence interval for the dataset.

This function takes two parameters: The phyloseq object (restroom) and the predictor variable (Environment).

Split the dataset by predictor variable
---------------------------------------

Two parameters are required to run this function: The phyloseq object (restroom) and the predictor variable (Environment).

``` r
per_variable_obj= pime.split.by.variable(restroom, "Environment")
per_variable_obj
#> $Restroom_F
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1858 taxa and 9 samples ]
#> sample_data() Sample Data:       [ 9 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1858 taxa by 7 taxonomic ranks ]
#> 
#> $Restroom_M
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1983 taxa and 9 samples ]
#> sample_data() Sample Data:       [ 9 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1983 taxa by 7 taxonomic ranks ]
```

Calculate the highest possible prevalence intervals
---------------------------------------------------

This function calculates prevalence for different intervals by increments of 5. The input file is the output from the pime.split.by.variable (per\_variable\_obj)

``` r
prevalences=pime.prevalence(per_variable_obj)
head(prevalences)
#> $`5`
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 3253 taxa and 18 samples ]
#> sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 3253 taxa by 7 taxonomic ranks ]
#> 
#> $`10`
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 3253 taxa and 18 samples ]
#> sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 3253 taxa by 7 taxonomic ranks ]
#> 
#> $`15`
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 593 taxa and 18 samples ]
#> sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 593 taxa by 7 taxonomic ranks ]
#> 
#> $`20`
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 593 taxa and 18 samples ]
#> sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 593 taxa by 7 taxonomic ranks ]
#> 
#> $`25`
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 222 taxa and 18 samples ]
#> sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 222 taxa by 7 taxonomic ranks ]
#> 
#> $`30`
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 222 taxa and 18 samples ]
#> sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 222 taxa by 7 taxonomic ranks ]
```

Calculate the best prevalence interval for the dataset
------------------------------------------------------

This function will return a table with Out of Bag error from random forests for each prevalence interval. The number of taxa and the number of remaining sequences for each prevalence interval are also computed. The best prevalence interval value provides the clearest separation of communities while still including a majority of the taxa in the analysis. If true differences are present. It will be represented by the first interval in which the OOB error rate is zero or close to zero. The input file is the list of prevalences generated by the pime.prevalence (prevalences) and the predictor variable ("Environment"").

``` r
set.seed(42)
best.prev=pime.best.prevalence(prevalences, "Environment")
#>        Interval OOB error rate (%) OTUs Nseqs
#>   Prevalence 5%              44.44 3253  9000
#>  Prevalence 10%              55.56 3253  9000
#>  Prevalence 15%              11.11  593  5438
#>  Prevalence 20%              11.11  593  5438
#>  Prevalence 25%              16.67  222  4370
#>  Prevalence 30%               5.56  222  4370
#>  Prevalence 35%               5.56  117  3835
#>  Prevalence 40%               5.56  117  3835
#>  Prevalence 45%               5.56   77  3531
#>  Prevalence 50%               5.56   77  3531
#>  Prevalence 55%               5.56   77  3531
#>  Prevalence 60%               5.56   45  3088
#>  Prevalence 65%                  0   45  3088
#>  Prevalence 70%                  0   26  2617
#>  Prevalence 75%                  0   26  2617
#>  Prevalence 80%                  0   16  2084
#>  Prevalence 85%                  0   16  2084
#>  Prevalence 90%                  0    4  1489
#>  Prevalence 95%                  0    4  1489
```

To get the table with OTU/ASV importance of the chosen prevalence interval. Pime keeps only those OTUs/ASVs with importance higher than overal mean.

``` r
imp50=best.prev$`Importance`$`Prevalence 50`

knitr::kable(imp50, format = 'markdown')
#>    SequenceID Restroom_F Restroom_M MeanDecreaseAccuracy MeanDecreaseGini
#> 1 denovo87919 0.07763333 0.07651905           0.07043333        1.2137966
#> 2 denovo22521 0.02973333 0.02083333           0.02368175        0.4528930
#> 3 denovo33035 0.02583333 0.01526667           0.01836587        0.4317984
#> 4  denovo6450 0.02536667 0.01143333           0.01563095        0.4134463
#> 5 denovo10912 0.01236667 0.01700000           0.01281190        0.2605335
#> 6 denovo29652 0.01713333 0.00730000           0.01125794        0.2908274
#>         Rank1             Rank2                  Rank3
#> 1 k__Bacteria p__Actinobacteria      c__Actinobacteria
#> 2  Unassigned              <NA>                   <NA>
#> 3 k__Bacteria p__Proteobacteria c__Gammaproteobacteria
#> 4  Unassigned              <NA>                   <NA>
#> 5  Unassigned              <NA>                   <NA>
#> 6 k__Bacteria p__Proteobacteria c__Alphaproteobacteria
#>                  Rank4                 Rank5            Rank6 Rank7
#> 1 o__Bifidobacteriales f__Bifidobacteriaceae             <NA>  <NA>
#> 2                 <NA>                  <NA>             <NA>  <NA>
#> 3   o__Pseudomonadales      f__Moraxellaceae g__Enhydrobacter   s__
#> 4                 <NA>                  <NA>             <NA>  <NA>
#> 5                 <NA>                  <NA>             <NA>  <NA>
#> 6   o__Caulobacterales   f__Caulobacteraceae              g__   s__

#To get the table with OOB error results.
#best.prev$`OOB error`
```

Within this dataset the best prevalence interval was 50%
--------------------------------------------------------

To obtain the phyloseq object at this cutoff use the following command.

``` r
prevalence.50 = prevalences$`50`
prevalence.50
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 77 taxa and 18 samples ]
#> sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 77 taxa by 7 taxonomic ranks ]
```

Estimating prediction error
---------------------------

To estimate error in prediction, we will use pime.error.prediction() to randomly assign treatments to samples and run random forests classification on each prevalence interval. The function returns a boxplot and a table with results of each classification error. For the purposes of this example we are running only 10 randomizations for saving time but we recommend at least 100 randomizations to obtain reliable results.

``` r
randomized=pime.error.prediction(restroom, "Environment", bootstrap = 10, parallel = TRUE, max.prev = 95)
randomized$Plot
randomized$'Table results'
```

It is also possible to estimate the variation of OOB error with each prevalence interval filtering. This is done by running the random forests classification for n times, determined by the user. This function will return a boxplot figure and a table for each classification error.

``` r
replicated.oob.error= pime.oob.replicate(prevalences, "Environment", bootstrap = 10, parallel = TRUE)
```

Getting Help
============

Please contact us if you need any help: *<contact@brmicrobiome.org>*

PIME Team
=========

Luiz F. W. Roesch (Universidade Federal do Pampa - Brazil)

Priscila T. Dobbler (Universidade Federal do Pampa - Brazil)

Victor S. Pylro (Universidade Federal de Lavras - Brazil)

Bryan Kolaczkowski (University of Florida - United States of America)

Jennifer C. Drew (University of Florida - United States of America)

Eric W. Triplett (University of Florida - United States of America)

Citation
========

Roesch et al. (2018), PIME: including the concept of prevalence for uncovering differences in microbiome noised data. Frontiers, submitted.
