
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

The first step in PIME is to define if the microbial community presents a high relative abundance of taxa with low prevalence, which is considered as noise in PIME analysis. This is calculated by random forests analysis. In this example we run PIME using the restroom dataset (<https://doi.org/10.1007%2Fs10482-017-0976-6>) against the metadata variable called Environment (a variable with two categories: men’s and women’s restroom).

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

This function calculates prevalence for different intervals by increments of 5. The input file is the output from the pime.split.by.variable(per\_variable\_obj)

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

In addition, it also returns the results from the random forests classification for each prevalence level. It includes SequenceID (OTU/ASV), Mean Decrease Accurracy (MDA) for each sample group, that is how much that SequenceID was important for classification of that group. The Mean Decrease Impurity (Gini Index) and taxonomy are also included.

To get the table with OTU/ASV importance of the chosen prevalence interval. Pime keeps only the top 30 OTUs/ASVs, whith highest MDA.

``` r
imp65=best.prev$`Importance`$`Prevalence 65`
knitr::kable(imp65, format="markdown") %>% kableExtra::kable_styling(full_width = F)
#> Warning in kableExtra::kable_styling(., full_width = F): Please specify
#> format in kable. kableExtra can customize either HTML or LaTeX outputs. See
#> https://haozhu233.github.io/kableExtra/ for details.
```

<table>
<colgroup>
<col width="6%" />
<col width="5%" />
<col width="5%" />
<col width="10%" />
<col width="8%" />
<col width="6%" />
<col width="9%" />
<col width="11%" />
<col width="10%" />
<col width="11%" />
<col width="10%" />
<col width="4%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">SequenceID</th>
<th align="right">Restroom_F</th>
<th align="right">Restroom_M</th>
<th align="right">MeanDecreaseAccuracy</th>
<th align="right">MeanDecreaseGini</th>
<th align="left">Rank1</th>
<th align="left">Rank2</th>
<th align="left">Rank3</th>
<th align="left">Rank4</th>
<th align="left">Rank5</th>
<th align="left">Rank6</th>
<th align="left">Rank7</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">denovo87919</td>
<td align="right">0.0757667</td>
<td align="right">0.0606667</td>
<td align="right">0.0612341</td>
<td align="right">1.0158054</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Actinobacteria</td>
<td align="left">c__Actinobacteria</td>
<td align="left">o__Bifidobacteriales</td>
<td align="left">f__Bifidobacteriaceae</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">denovo22521</td>
<td align="right">0.0356000</td>
<td align="right">0.0179000</td>
<td align="right">0.0244095</td>
<td align="right">0.5767724</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">denovo6450</td>
<td align="right">0.0330667</td>
<td align="right">0.0196000</td>
<td align="right">0.0237397</td>
<td align="right">0.5189465</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Gammaproteobacteria</td>
<td align="left">o__Pseudomonadales</td>
<td align="left">f__Moraxellaceae</td>
<td align="left">g__Enhydrobacter</td>
<td align="left">s__</td>
</tr>
<tr class="even">
<td align="left">denovo1419</td>
<td align="right">0.0150524</td>
<td align="right">0.0330333</td>
<td align="right">0.0215294</td>
<td align="right">0.4957805</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">denovo65044</td>
<td align="right">0.0141667</td>
<td align="right">0.0309000</td>
<td align="right">0.0206190</td>
<td align="right">0.5613959</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">denovo23087</td>
<td align="right">0.0203333</td>
<td align="right">0.0100524</td>
<td align="right">0.0144373</td>
<td align="right">0.3774329</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Alphaproteobacteria</td>
<td align="left">o__Caulobacterales</td>
<td align="left">f__Caulobacteraceae</td>
<td align="left">g__</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo33035</td>
<td align="right">0.0180333</td>
<td align="right">0.0089333</td>
<td align="right">0.0125206</td>
<td align="right">0.3455253</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">denovo10912</td>
<td align="right">0.0059667</td>
<td align="right">0.0199667</td>
<td align="right">0.0123397</td>
<td align="right">0.3517237</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Actinobacteria</td>
<td align="left">c__Actinobacteria</td>
<td align="left">o__Actinomycetales</td>
<td align="left">f__Micrococcaceae</td>
<td align="left">g__</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo88424</td>
<td align="right">0.0061333</td>
<td align="right">0.0151667</td>
<td align="right">0.0100698</td>
<td align="right">0.3498350</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">denovo39172</td>
<td align="right">0.0103000</td>
<td align="right">0.0067333</td>
<td align="right">0.0082071</td>
<td align="right">0.2675950</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Alphaproteobacteria</td>
<td align="left">o__Rhizobiales</td>
<td align="left">f__Rhizobiaceae</td>
<td align="left">g__Agrobacterium</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo23140</td>
<td align="right">0.0053333</td>
<td align="right">0.0085000</td>
<td align="right">0.0071500</td>
<td align="right">0.1912480</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">denovo47015</td>
<td align="right">0.0025000</td>
<td align="right">0.0155667</td>
<td align="right">0.0067405</td>
<td align="right">0.1932127</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Alphaproteobacteria</td>
<td align="left">o__Rhizobiales</td>
<td align="left">f__Methylobacteriaceae</td>
<td align="left">g__</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo29652</td>
<td align="right">0.0134000</td>
<td align="right">0.0019000</td>
<td align="right">0.0057238</td>
<td align="right">0.2779126</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">denovo3089</td>
<td align="right">0.0035667</td>
<td align="right">0.0094667</td>
<td align="right">0.0055429</td>
<td align="right">0.1930492</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="odd">
<td align="left">denovo33672</td>
<td align="right">0.0112333</td>
<td align="right">0.0031667</td>
<td align="right">0.0053452</td>
<td align="right">0.1983517</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Firmicutes</td>
<td align="left">c__Bacilli</td>
<td align="left">o__Lactobacillales</td>
<td align="left">f__Lactobacillaceae</td>
<td align="left">g__Lactobacillus</td>
<td align="left">s__iners</td>
</tr>
<tr class="even">
<td align="left">denovo61097</td>
<td align="right">0.0036667</td>
<td align="right">0.0068333</td>
<td align="right">0.0045444</td>
<td align="right">0.1276202</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Actinobacteria</td>
<td align="left">c__Actinobacteria</td>
<td align="left">o__Actinomycetales</td>
<td align="left">f__Corynebacteriaceae</td>
<td align="left">g__Corynebacterium</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo18505</td>
<td align="right">0.0016333</td>
<td align="right">0.0075333</td>
<td align="right">0.0045437</td>
<td align="right">0.1735404</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Alphaproteobacteria</td>
<td align="left">o__Sphingomonadales</td>
<td align="left">f__Sphingomonadaceae</td>
<td align="left">g__Sphingomonas</td>
<td align="left">s__</td>
</tr>
<tr class="even">
<td align="left">denovo36101</td>
<td align="right">0.0005000</td>
<td align="right">0.0077333</td>
<td align="right">0.0040000</td>
<td align="right">0.1843757</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Betaproteobacteria</td>
<td align="left">o__Burkholderiales</td>
<td align="left">f__Burkholderiaceae</td>
<td align="left">g__Burkholderia</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo43216</td>
<td align="right">0.0070000</td>
<td align="right">0.0015333</td>
<td align="right">0.0038524</td>
<td align="right">0.1996569</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Gammaproteobacteria</td>
<td align="left">o__Pseudomonadales</td>
<td align="left">f__Pseudomonadaceae</td>
<td align="left">g__Pseudomonas</td>
<td align="left">s__</td>
</tr>
<tr class="even">
<td align="left">denovo85345</td>
<td align="right">0.0050667</td>
<td align="right">0.0031333</td>
<td align="right">0.0033810</td>
<td align="right">0.1495234</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Alphaproteobacteria</td>
<td align="left">o__Sphingomonadales</td>
<td align="left">f__Sphingomonadaceae</td>
<td align="left">g__Sphingomonas</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo82372</td>
<td align="right">0.0030000</td>
<td align="right">0.0045667</td>
<td align="right">0.0033214</td>
<td align="right">0.2517545</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Actinobacteria</td>
<td align="left">c__Actinobacteria</td>
<td align="left">o__Actinomycetales</td>
<td align="left">f__Kineosporiaceae</td>
<td align="left">g__Kineococcus</td>
<td align="left">s__</td>
</tr>
<tr class="even">
<td align="left">denovo12214</td>
<td align="right">0.0054667</td>
<td align="right">0.0016667</td>
<td align="right">0.0032929</td>
<td align="right">0.1342444</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Gammaproteobacteria</td>
<td align="left">o__Pseudomonadales</td>
<td align="left">f__Pseudomonadaceae</td>
<td align="left">g__Pseudomonas</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo33246</td>
<td align="right">0.0047333</td>
<td align="right">0.0026667</td>
<td align="right">0.0030079</td>
<td align="right">0.1314324</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">denovo49990</td>
<td align="right">0.0061333</td>
<td align="right">0.0011667</td>
<td align="right">0.0029500</td>
<td align="right">0.1045538</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Actinobacteria</td>
<td align="left">c__Actinobacteria</td>
<td align="left">o__Actinomycetales</td>
<td align="left">f__Microbacteriaceae</td>
<td align="left">g__Curtobacterium</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo93818</td>
<td align="right">0.0070000</td>
<td align="right">-0.0000333</td>
<td align="right">0.0025238</td>
<td align="right">0.1543636</td>
<td align="left">Unassigned</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
<td align="left">NA</td>
</tr>
<tr class="even">
<td align="left">denovo32646</td>
<td align="right">0.0023667</td>
<td align="right">0.0015667</td>
<td align="right">0.0022738</td>
<td align="right">0.1121110</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Bacteroidetes</td>
<td align="left">c__Flavobacteriia</td>
<td align="left">o__Flavobacteriales</td>
<td align="left">f__[Weeksellaceae]</td>
<td align="left">g__Chryseobacterium</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo41355</td>
<td align="right">0.0010667</td>
<td align="right">0.0025667</td>
<td align="right">0.0017579</td>
<td align="right">0.0509914</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Cyanobacteria</td>
<td align="left">c__Chloroplast</td>
<td align="left">o__Streptophyta</td>
<td align="left">f__</td>
<td align="left">g__</td>
<td align="left">s__</td>
</tr>
<tr class="even">
<td align="left">denovo90494</td>
<td align="right">0.0053333</td>
<td align="right">-0.0003333</td>
<td align="right">0.0016000</td>
<td align="right">0.0998579</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Firmicutes</td>
<td align="left">c__Bacilli</td>
<td align="left">o__Lactobacillales</td>
<td align="left">f__Leuconostocaceae</td>
<td align="left">g__</td>
<td align="left">s__</td>
</tr>
<tr class="odd">
<td align="left">denovo39338</td>
<td align="right">0.0004333</td>
<td align="right">0.0032333</td>
<td align="right">0.0015286</td>
<td align="right">0.1188592</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Proteobacteria</td>
<td align="left">c__Alphaproteobacteria</td>
<td align="left">o__Sphingomonadales</td>
<td align="left">f__Sphingomonadaceae</td>
<td align="left">g__Sphingomonas</td>
<td align="left">s__</td>
</tr>
<tr class="even">
<td align="left">denovo61390</td>
<td align="right">0.0018000</td>
<td align="right">0.0006667</td>
<td align="right">0.0012857</td>
<td align="right">0.0689921</td>
<td align="left">k__Bacteria</td>
<td align="left">p__Actinobacteria</td>
<td align="left">c__Actinobacteria</td>
<td align="left">o__Actinomycetales</td>
<td align="left">f__Geodermatophilaceae</td>
<td align="left">g__</td>
<td align="left">s__</td>
</tr>
</tbody>
</table>

``` r
#To get the table with OOB error results.
#best.prev$`OOB error`
```

Within this dataset the best prevalence interval was 65%
--------------------------------------------------------

To obtain the phyloseq object at this cutoff use the following command.

``` r
prevalence.65 = prevalences$`65`
prevalence.65
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 18 samples ]
#> sample_data() Sample Data:       [ 18 samples by 6 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 7 taxonomic ranks ]
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
