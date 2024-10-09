---
title: "Special covariates for the X chromosome"
teaching: 15
exercises: 15
---

:::::::::::::::::::::::::::::::::::::: questions 

- "How do I find the chromosome X covariates for a cross?"

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Get the X covariates for a cross.

::::::::::::::::::::::::::::::::::::::::::::::::



The X chromosome must be treated differently from the autosomes in the mapping 
of quantitative trait loci (QTL). If the X chromosome is treated like an 
autosome, a sex difference in a phenotype, such as weight or height, can lead to 
spurious linkage on the X chromosome. The X chromosome varies depending on the 
sex of the animal and the direction of the cross, so accounting for these 
covariates is important under the null hypothesis of no QTL, to avoid spurious 
evidence of linkage. (See 
[Broman et al. (2006) Genetics 174:2151-2158](http://www.genetics.org/content/174/4/2151.long).) 
Specifically, 
[see the figures](https://www.genetics.org/content/174/4/2151.figures-only) 
for the behavior of the X chromosome in a backcross and in an intercross.

The particular X chromosome covariates depends on the cross, and can be obtained 
with the function `get_x_covar()`. In the iron data, sex is indicated as 0 for 
females and 1 for males. For cross direction, see Figure 2 of 
[Broman et al. (2006) Genetics 174:2151-2158](http://www.genetics.org/content/174/4/2151.long).  
For cross direction, all samples are filled in with 0 except for females from 
the *reverse* direction, who are indicated with 1.



``` r
Xcovar <- get_x_covar(iron)
head(Xcovar)
```

``` output
  sex direction
1   1         0
2   1         0
3   1         0
4   1         0
5   1         0
6   1         0
```

::::::::::::::::::::::::::::::::::::: keypoints 

- "The X chromosome requires special treatment in QTL mapping."
- "Special covariates such as sex should be included to avoid spurious evidence 
of linkage."

::::::::::::::::::::::::::::::::::::::::::::::::
