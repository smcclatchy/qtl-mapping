---
title: "Performing a genome scan"
teaching: 30
exercises: 30
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do I perform a genome scan?
- How do I plot a genome scan?
- How do additive covariates differ from interactive covariates?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Map one trait using additive covariates.
- Map the sample trait using additive and interactive covariates.
- Plot a genome scan.

::::::::::::::::::::::::::::::::::::::::::::::::




The freely available chapter on 
[single-QTL analysis](http://www.rqtl.org/book/rqtlbook_ch04.pdf) 
from Broman and Sen's 
[A Guide to QTL Mapping with R/qtl](http://www.rqtl.org/book/) describes 
different methods for QTL analysis. Two of these methods are described here 
using data from an experiment on hypertension in the mouse 
[Sugiyama et al., Genomics 71:70-77, 2001](https://s3.amazonaws.com/academia.edu.documents/45963759/geno.2000.640120160526-29022-36mpgg.pdf?AWSAccessKeyId=AKIAIWOWYYGZ2Y53UL3A&Expires=1513786158&Signature=rtodlYwe0LDmYZFOm1ejvZjZhQ0%3D&response-content-disposition=inline%3B%20filename%3DConcordance_of_murine_quantitative_trait.pdf). 
The study employed a backcross between two mouse strains resulting in two 
possible genotypes - AA for homozygotes, and AB for heterozygotes. 

Linear regression can be employed to identify presence of QTL in a cross. To 
identify QTL using regression, we compare the fit for two models: 1) the null 
hypothesis that there are no QTL anywhere in the genome; and 2) the alternative 
hypothesis that there is a QTL near a specific position. A sloped line indicates 
that there is a difference in mean phenotype between the two genotype groups, 
and that a QTL is present. A line with no slope indicates that there is no 
difference in mean phenotype between the two groups, and that no QTL exists. 
Regression aims to find the line of best fit to the data. In the case of a 
backcross with only two genotypes, a t-test is performed at the marker to 
determine whether the difference in phenotype means is zero.

![adapted from Broman & Sen, 2009](fig/nullvalt.png)

To find the line of best fit, the residuals or errors are calculated, then 
squared for each data point.

![Squared error for a single data point](fig/althypothesis.png)
![The line of best fit minimizes the sum of squared errors](fig/nullhypothesis.png)

The line of best fit will be the one that minimizes the sum of squared 
residuals, which maximizes the likelihood of the data. 

Marker regression produces a LOD (logarithm of odds) score comparing the null 
hypothesis to the alternative. The LOD score is calculated using the sum of 
squared residuals for the null and alternative hypotheses. The LOD score is the 
difference between the log10 likelihood of the null hypothesis and the log10 
likelihood of the alternative hypothesis. It is related to the regression model 
above by identifying the line of best fit to the data. A higher LOD score 
indicates greater likelihood of the alternative hypothesis. A LOD score closer 
to zero favors the null hypothesis. 

Marker regression can identify the existence and effect of a QTL by comparing 
means between groups, however, it requires known marker genotypes and can't 
identify QTL in between typed markers. To identify QTL between typed markers, we 
use Haley-Knott regression. After 
[calculating genotype probabilities](https://smcclatchy.github.io/mapping/03-calc-genoprob/), 
we can regress the phenotypes for animals of unknown genotype on these 
conditional genotype probabilities (conditional on known marker genotypes). In 
Haley-Knott regression, phenotype values can be plotted and a regression line 
drawn through the phenotype mean for the untyped individuals.

![](fig/hk-regress.png)

As shown by the green circle in the figure, an individual of unknown genotype is 
placed between known genotypes according to the probability of its genotype 
being AA or AB. In this case, the probability of this individual having genotype
AA is 0.6, and the probability of having genotype AB is 0.4.

To perform a genome scan by Haley-Knott regression
([Haley and Knott 1992](https://www.ncbi.nlm.nih.gov/pubmed/16718932)),
use the function `scan1()`.  `scan1()` takes as input the genotype 
probabilities, a matrix of phenotypes, and then optional additive and 
interactive covariates. Another option 
is to provide a vector of weights.

## Additive Genome Scan

There are two potential covariates in the Attie data set. Let's look at the
top of the covariates in the `cross` object.


``` r
head(cross$covar)
```

``` output
             Sex pgm adipose_batch gastroc_batch hypo_batch islet_batch
Mouse3051   Male   1    12/19/2007     8/11/2008 11/26/2007  11/28/2007
Mouse3551   Male   1    12/19/2007     8/12/2008 11/27/2007  12/03/2007
Mouse3430   Male   1    12/19/2007     8/12/2008 11/27/2007  12/03/2007
Mouse3476   Male   1    12/19/2007     8/12/2008 11/27/2007  12/03/2007
Mouse3414   Male   1    12/18/2007     8/11/2008 11/26/2007  11/28/2007
Mouse3145 Female   1    12/19/2007     8/11/2008         NA  11/28/2007
          kidney_batch liver_batch
Mouse3051   07/15/2008       other
Mouse3551   07/16/2008  12/10/2007
Mouse3430   07/15/2008  12/05/2007
Mouse3476   07/16/2008  12/05/2007
Mouse3414   07/15/2008  12/05/2007
Mouse3145   07/15/2008  12/10/2007
```

`Sex` is potential covariate. It is a good idea to always include
sex in any analysis. Examples of other covariates might be age, diet, treatment,
or experimental batch. It is worth taking time to identify covariates that may
affect your results.



``` r
cross$covar$Sex <- factor(cross$covar$Sex)

addcovar <- model.matrix(~Sex, data = cross$covar)[,-1, drop = FALSE]
```

When we perform a genome scan with additive covariates, we are searching for loci
that have the same effect in both covariate groups. In this case, we are 
searching for loci that affect females and male in the same way.

<!-- DMG: We need a plot of what an additive effect would look like at one marker. -->



``` r
lod_add <- scan1(genoprobs = probs, 
                 pheno     = cross$pheno[,'log10_insulin_10wk'], 
                 addcovar  = addcovar)
```

On a multi-core machine, you can get some speed-up via the `cores` argument, as 
with `calc_genoprob()` and `calc_kinship()`.


``` r
lod_add <- scan1(genoprobs = probs, 
                 pheno     = cross$pheno[,'log10_insulin_10wk', drop = FALSE], 
                 addcovar  = addcovar,
                 cores     = 4)
```

The output of `scan1()` is a matrix of LOD scores, positions &times; phenotypes. 

Take a look at the first ten rows of the scan object. The numerical values are 
the LOD scores for the marker named at the beginning of the row. 
LOD values are shown for circulating insulin.


``` r
head(lod_add, n = 10)
```

``` output
               pheno1
rs13475697 0.04674829
rs3681603  0.04674829
rs13475703 0.04680494
rs13475710 0.12953382
rs6367205  0.13734728
rs13475716 0.13735534
rs13475717 0.13735534
rs13475719 0.13735534
rs13459050 0.13735534
rs3680898  0.13735541
```

The function `plot_scan1()` can be used to plot the LOD curves. Use the argument 
`lodcolumn` to indicate which column to plot.


``` r
plot_scan1(lod_add, 
           map  = cross$pmap,
           main = 'log(insulin): 10 weeks')
```

<img src="fig/perform-genome-scan-rendered-plot_add_lod-1.png" style="display: block; margin: auto;" />

The LOD plot for insulin shows several peaks, with the largest peak on 
chromosome 2. There are smaller peaks on other chromosomes. Which of these peaks is 
significant, and why? We'll evaluate the significance of genome scan results in 
a later episode in 
[performing a permutation test](https://smcclatchy.github.io/mapping/10-perform-perm-test/).


<!-- DMG: Suggestion: let's have the exercises have them plot the LOD score for
one chromosome, the find the marker with the maximum LOD score and try to 
find it in the physical map. -->

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1

Use the `sort()` function to sort the LOD scores for liver. Hint: run `dim(out)` 
for the row and column dimensions, or `colnames(out)` for the column names.   
Which marker has the highest LOD score? Which genotyped marker has the 
highest LOD score? What chromosome number are they on? 

:::::::::::::::::::::::: solution 

`sort(out[,1])` for column 1 or `sort(out[,"liver"])` for the column named 
"liver".  
The pseudomarker with the largest score is c16.loc29, with a LOD of 7.68. The 
genotyped marker with the largest LOD is D16Mit30 with a score of 7.28. Both are 
located on chromosome 16.

:::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 2

Use the `sort()` function to sort the LOD scores for spleen. Which pseudomarker 
has the highest LOD score? Which genotyped marker has the highest LOD score? 
What chromosome number are they on? 

:::::::::::::::::::::::: solution 

`sort(out[,2])` or `sort(out[,"spleen"])`
The pseudomarker with the largest score is c9.loc57, with a LOD of 12.1. The 
genotyped marker with the largest LOD is D9Mit182 with a score of 10.4. Both are 
located on chromosome 9.

:::::::::::::::::::::::::::::::::

## Challenge 3

Plot the LOD scores for spleen. Does the genome scan for spleen share any 
large-ish peaks with the scan for liver?

:::::::::::::::::::::::: solution 

`plot_scan1(out, map = map, lodcolumn = "spleen")`
Both liver and spleen genome scans share a peak on chromosome 8 with a LOD score 
near 4.

:::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::

## Interactive Genome Scan

Above, we mapped insulin levels using sex as an additive covariate and searched
for loci where both sexes had the same effect. But what if the two sexes have
different effects? The we would like to map in a way that allows each sex to
have different effects. We do this using an interactive genome scan.

You should always include the interactive covariates as additive covariates as
well. In this case, we only have sex as a covariate, so we can use the additive
covariate matrix. For clarity, we will make a copy and name it for the 
interactive covariates.


``` r
intcovar = addcovar
```



``` r
lod_int <- scan1(genoprobs = probs, 
                 pheno     = cross$pheno[,'log10_insulin_10wk'], 
                 addcovar  = addcovar,
                 intcovar  = intcovar)
```



``` r
plot_scan1(x = lod_int, map = cross$pmap, main = 'log10_insulin_10wk')
```

<img src="fig/perform-genome-scan-rendered-plot_int-1.png" style="display: block; margin: auto;" />

It is difficult to tell if there is a difference in LOD scores between the 
additive and interactive scans. To resolve this, we can plot both genome scans
in the same plo using the "add = TRUE" argument.


``` r
plot_scan1(x = lod_int, map = cross$pmap, main = 'log10_insulin_10wk')
plot_scan1(x = lod_add, map = cross$pmap, col = 'blue', add = TRUE)
```

<img src="fig/perform-genome-scan-rendered-plot_add_int_lod-1.png" style="display: block; margin: auto;" />

It is still difficult to tell whether any peaks differ by sex. Another way 
to view the plot is to plot the difference between the interactive and additive
scans.


``` r
plot_scan1(x = lod_int, map = cross$pmap, main = 'log10_insulin_10wk')
plot_scan1(x = lod_add, map = cross$pmap, col = 'blue', add = TRUE)
plot_scan1(x = lod_int - lod_add, map = cross$pmap, col = 'red', add = TRUE)
```

<img src="fig/perform-genome-scan-rendered-plot_add_int_diff_lod-1.png" style="display: block; margin: auto;" />

While it was important to look at the effect of sex on the trait, in this 
experiment, there do not appear to be any sex-specific peaks.

<!-- DMG: Add challenge asking students why we didn't map and plot females and
males separately. Maybe they can even try it. -->

::::::::::::::::::::::::::::::::::::: keypoints 

- A qtl2 genome scan requires genotype probabilities and a phenotype matrix.
- The output from a genome scan contains a LOD score matrix, map positions, and phenotypes.
- LOD curve plots for a genome scan can be viewed with plot_scan1().

::::::::::::::::::::::::::::::::::::::::::::::::

