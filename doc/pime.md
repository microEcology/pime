---
title: "DMPC"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{DMPC}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
## How to install and use DMPC package

To install DMPC Package first install the devtools package. Then run
install_github("vpylro/DMPC")
DMPC uses a Phyloseq object as input. 

##This function will be used to predict the likelihood of KIDMED index diet categories to influence the saliva's microbiome.
The prediction is performed by randonForests on full dataset. 

```r
library(DMPC)
data("restroom")
DMPC.OOB.error(restroom, "Environment")
```

```
##       OOB 
## 0.5555556
```
If the OOB error rate is zero, the dataset present large differences, and DMPC might not be necessary. 
If the OOB error rate is greater than zero, next functions will find the best prevalence cutoff for the dataset.
This function takes two parameters: The phyloseq object (restroom) and the predictor variable (Environment).

## Split the dataset by predictor variable
Two parameters are required to run this function: The phyloseq object (restroom) and the predictor variable (Environment).

```r
per_variable_obj= DMPC.split.by.variable(restroom, "Environment")
per_variable_obj
```

```
## $Restroom_F
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1858 taxa and 9 samples ]
## sample_data() Sample Data:       [ 9 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 1858 taxa by 7 taxonomic ranks ]
## 
## $Restroom_M
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1983 taxa and 9 samples ]
## sample_data() Sample Data:       [ 9 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 1983 taxa by 7 taxonomic ranks ]
```

## Calculate the highest possible prevalence cutoffs
This function calculates prevalence for different cutoff values by increments of 5. 
The input file is the output from the DMPC.split.by.variable (per_variable_obj)

```r
prevalences=DMPC.prevalence(per_variable_obj)
prevalences
```

```
## $`5`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 3253 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 3253 taxa by 7 taxonomic ranks ]
## 
## $`10`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 3253 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 3253 taxa by 7 taxonomic ranks ]
## 
## $`15`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 593 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 593 taxa by 7 taxonomic ranks ]
## 
## $`20`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 593 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 593 taxa by 7 taxonomic ranks ]
## 
## $`25`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 222 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 222 taxa by 7 taxonomic ranks ]
## 
## $`30`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 222 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 222 taxa by 7 taxonomic ranks ]
## 
## $`35`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 117 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 117 taxa by 7 taxonomic ranks ]
## 
## $`40`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 117 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 117 taxa by 7 taxonomic ranks ]
## 
## $`45`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 77 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 77 taxa by 7 taxonomic ranks ]
## 
## $`50`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 77 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 77 taxa by 7 taxonomic ranks ]
## 
## $`55`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 77 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 77 taxa by 7 taxonomic ranks ]
## 
## $`60`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 45 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 45 taxa by 7 taxonomic ranks ]
## 
## $`65`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 45 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 45 taxa by 7 taxonomic ranks ]
## 
## $`70`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 26 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 26 taxa by 7 taxonomic ranks ]
## 
## $`75`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 26 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 26 taxa by 7 taxonomic ranks ]
## 
## $`80`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 16 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 16 taxa by 7 taxonomic ranks ]
## 
## $`85`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 16 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 16 taxa by 7 taxonomic ranks ]
## 
## $`90`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 4 taxa by 7 taxonomic ranks ]
## 
## $`95`
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 4 taxa by 7 taxonomic ranks ]
```

## Calculates the best prevalence cutoff for the dataset 
This function will return a table with PERMANOVA  and RandomForests results for each prevalence cutoff.
The best prevalence cutoff value provides the clearest separation of communities while still including a majority of the taxa in the analysis. If true differences are present.
It will be represented by the first cutoff in which the OOB error rate is zero.
The input file is the list of prevalences generated by the DMPC.prevalence (prevalences) and the predictor variable (Environment).

```r
DMPC.best.prevalence(prevalences, "Environment")
```

```
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
## [1] "Calculating..."
## [1] "Done"
```

```
##               Df SumsOfSqs   MeanSqs  F.Model         R2 Pr(>F)        OOB
## Prevalence 5   1 0.3635858 0.3635858 1.098232 0.06423071  0.333 0.44444444
## Prevalence 10  1 0.3635858 0.3635858 1.098232 0.06423071  0.314 0.55555556
## Prevalence 15  1 0.4991134 0.4991134 2.025911 0.11238882  0.006 0.22222222
## Prevalence 20  1 0.4991134 0.4991134 2.025911 0.11238882  0.015 0.11111111
## Prevalence 25  1 0.4840748 0.4840748 2.251734 0.12337094  0.011 0.11111111
## Prevalence 30  1 0.4840748 0.4840748 2.251734 0.12337094  0.011 0.11111111
## Prevalence 35  1 0.5052119 0.5052119 2.501967 0.13522707  0.005 0.05555556
## Prevalence 40  1 0.5052119 0.5052119 2.501967 0.13522707  0.010 0.05555556
## Prevalence 45  1 0.5303799 0.5303799 2.684016 0.14365307  0.006 0.05555556
## Prevalence 50  1 0.5303799 0.5303799 2.684016 0.14365307  0.007 0.05555556
## Prevalence 55  1 0.5303799 0.5303799 2.684016 0.14365307  0.008 0.05555556
## Prevalence 60  1 0.5596027 0.5596027 2.895444 0.15323503  0.009 0.00000000
## Prevalence 65  1 0.5596027 0.5596027 2.895444 0.15323503  0.006 0.00000000
## Prevalence 70  1 0.5679537 0.5679537 3.023264 0.15892456  0.006 0.00000000
## Prevalence 75  1 0.5679537 0.5679537 3.023264 0.15892456  0.012 0.00000000
## Prevalence 80  1 1.0083660 1.0083660 5.686393 0.26221019  0.001 0.00000000
## Prevalence 85  1 1.0083660 1.0083660 5.686393 0.26221019  0.001 0.00000000
## Prevalence 90  1 0.8059285 0.8059285 4.653893 0.22532765  0.006 0.00000000
## Prevalence 95  1 0.8059285 0.8059285 4.653893 0.22532765  0.005 0.00000000
##               Cutoff OTUs Nseqs
## Prevalence 5      5% 3253  9000
## Prevalence 10    10% 3253  9000
## Prevalence 15    15%  593  5438
## Prevalence 20    20%  593  5438
## Prevalence 25    25%  222  4370
## Prevalence 30    30%  222  4370
## Prevalence 35    35%  117  3835
## Prevalence 40    40%  117  3835
## Prevalence 45    45%   77  3531
## Prevalence 50    50%   77  3531
## Prevalence 55    55%   77  3531
## Prevalence 60    60%   45  3088
## Prevalence 65    65%   45  3088
## Prevalence 70    70%   26  2617
## Prevalence 75    75%   26  2617
## Prevalence 80    80%   16  2084
## Prevalence 85    85%   16  2084
## Prevalence 90    90%    4  1489
## Prevalence 95    95%    4  1489
```
## Within this dataset the best prevalence cutoff was 60%
To obtain the phyloseq object at this cutoff use the following command.


```r
prevalence.60 = prevalences$`60`
prevalence.60
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 45 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
## tax_table()   Taxonomy Table:    [ 45 taxa by 7 taxonomic ranks ]
```
## FDR
If needed, its possible to run a False Discovery Rate.

```r
DMPC.FDR(restroom, "Environment", bootstrap = 100, method.dist = "bray")
```

```
## [1] FDR value
##    5 
## 0.04 
## [1] FDR value
##    5   10 
## 0.04 0.03 
## [1] FDR value
##    5   10   15 
## 0.04 0.03 0.06 
## [1] FDR value
##    5   10   15   20 
## 0.04 0.03 0.06 0.05 
## [1] FDR value
##    5   10   15   20   25 
## 0.04 0.03 0.06 0.05 0.06 
## [1] FDR value
##    5   10   15   20   25   30 
## 0.04 0.03 0.06 0.05 0.06 0.06 
## [1] FDR value
##    5   10   15   20   25   30   35 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 
## [1] FDR value
##    5   10   15   20   25   30   35   40 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50   55 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 0.03 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50   55   60 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 0.03 0.04 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50   55   60   65 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 0.03 0.04 0.04 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50   55   60   65   70 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 0.03 0.04 0.04 0.07 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50   55   60   65   70   75 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 0.03 0.04 0.04 0.07 0.06 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50   55   60   65   70   75 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 0.03 0.04 0.04 0.07 0.06 
##   80 
## 0.04 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50   55   60   65   70   75 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 0.03 0.04 0.04 0.07 0.06 
##   80   85 
## 0.04 0.05 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50   55   60   65   70   75 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 0.03 0.04 0.04 0.07 0.06 
##   80   85   90 
## 0.04 0.05 0.04 
## [1] FDR value
##    5   10   15   20   25   30   35   40   45   50   55   60   65   70   75 
## 0.04 0.03 0.06 0.05 0.06 0.06 0.06 0.06 0.04 0.07 0.03 0.04 0.04 0.07 0.06 
##   80   85   90   95 
## 0.04 0.05 0.04 0.02
```

```
##               FDRv
## Prevalence5%  0.04
## Prevalence10% 0.03
## Prevalence15% 0.06
## Prevalence20% 0.05
## Prevalence25% 0.06
## Prevalence30% 0.06
## Prevalence35% 0.06
## Prevalence40% 0.06
## Prevalence45% 0.04
## Prevalence50% 0.07
## Prevalence55% 0.03
## Prevalence60% 0.04
## Prevalence65% 0.04
## Prevalence70% 0.07
## Prevalence75% 0.06
## Prevalence80% 0.04
## Prevalence85% 0.05
## Prevalence90% 0.04
## Prevalence95% 0.02
```
