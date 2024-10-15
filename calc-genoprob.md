---
title: "Calculating Genotype Probabilities"
teaching: 30
exercises: 30
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do I calculate QTL at positions between genotyped markers?
- How do I calculate QTL genotype probabilities?
- How do I calculate allele probabilities?
- How can I speed up calculations if I have a large data set?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- To explain why the first step in QTL analysis is to calculate genotype 
probabilities.
- To calculate genotype probabilities.

::::::::::::::::::::::::::::::::::::::::::::::::

The first task in QTL analysis is to calculate conditional genotype 
probabilities, given the observed marker data, at each putative QTL position. 
For example, the first step would be to determine the probabilities for 
genotypes SS and SB at the locus indicated below.

![adapted from Broman & Sen, 2009](fig/unknown_genotype.png)

The `calc_genoprob()` function calculates QTL genotype probabilities, 
conditional on the available marker data. These are needed for most of the QTL 
mapping functions. The result is returned as a list of three-dimensional arrays 
(one per chromosome). Each 3d array of probabilities is arranged as individuals 
$\times$ genotypes $\times$ positions.

![Three-dimensional genotype probabilities array](fig/threeD_array.png){alt='Figure showing three-dimensional array of genoprobs'}

![See this page for a graphical review of data structures in R](http://venus.ifca.unican.es/Rintro/_images/dataStructuresNew.png).  

We'll use the
[Attie BL6/BTBR dataset](https://thejacksonlaboratory.box.com/shared/static/svw7ivp5hhmd7vb8fy26tc53h7r85wez.zip)
from 
[Tian et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4649649/)
(an intercross) as an example. In this study, circulating insulin levels were 
measured in an F2 cross between mouse strains C57BL/6J and BTBTR T+ <tf>. 

First, we will load in the [qtl2](https://kbroman.org/qtl2/) library, which 
provides the functions that we will use for QTL analysis.


``` r
suppressPackageStartupMessages(library(qtl2))
```



``` r
cross <- read_cross2(file = 'data/attie_b6btbr_grcm39/attie_control.json')
```

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: instructor

We need the following block for the site to build on Github. The students do
not need to see or run the next block.

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



To load your own data from your machine, you would use the file path to your 
data files. For example, if the file path to your data files is 
`/Users/myUserName/qtlProject/data`, the command to load your data would look 
like this:


``` r
myQTLdata <- read_cross2(file = "/Users/myUserName/qtlProject/data/myqtldata.json" )
```

The JSON file contains all control information for your data, including names of 
data files, cross type, column specifications for sex and cross information, and 
more. This can also be in YAML format. Alternatively, all data files can be 
zipped together for loading.


``` r
myQTLdata <- read_cross2(file = "/Users/myUserName/qtlProject/data/myqtldata.zip" )
```

Back to the BTBR data. Now look at a summary of the cross data and the names of 
each variable within the data.


``` r
summary(cross)
```

``` output
Object of class cross2 (crosstype "f2")

Total individuals             490
No. genotyped individuals     490
No. phenotyped individuals    490
No. with both geno & pheno    490

No. phenotypes                  3
No. covariates                  8
No. phenotype covariates        0

No. chromosomes                20
Total markers                2057

No. markers by chr:
  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19   X 
156 135 157 126 125 102 109  91  93 123 124 116 116  91 102  66  60  95  50  20 
```

``` r
names(cross)
```

``` output
 [1] "crosstype"  "geno"       "gmap"       "pmap"       "pheno"     
 [6] "covar"      "is_x_chr"   "is_female"  "cross_info" "alleles"   
```

Have a look at the markers listed in the genetic map, `gmap`. Markers are listed 
by chromosome and described by cM position. View only the markers on the first 
several chromosomes.


``` r
head(cross$gmap)
```

``` output
$`1`
rs13475697  rs3681603 rs13475703 rs13475710  rs6367205 rs13475716 rs13475717 
 0.1881141  0.1920975  0.4167755  0.6488793  0.6555814  0.6638576  0.6676198 
rs13475719 rs13459050  rs3680898 rs13475727 rs13475728 rs13475729 rs13475731 
 0.6711377  0.6749344  0.6775292  1.8149573  1.9596637  2.3456569  2.7186389 
rs13475737 rs13475744  rs6397513 rs13475747 rs13475748 rs13475749 rs13475750 
 3.1059517  3.8222865  4.3094607  4.3120150  4.5098582  4.8154609  4.8853505 
rs13475751 rs13475752 rs13475762 rs13475764 rs13475765 rs13475768 rs13475769 
 4.8869793  4.8902179  7.2954871  8.2102887  8.3708197  8.7178703  8.8859153 
rs13475771  rs6384194 rs13475790  rs3676270 rs13475794 rs13475801  rs4222269 
 9.1374722  9.9295192  9.9970634 10.1508878 10.3962716 11.5981956 11.9606369 
 rs6387241 rs13475822 rs13475824 rs13475826 rs13475827 rs13475834 rs13475880 
16.8770742 16.9815396 17.4434784 18.0866148 18.6276972 19.2288050 27.4056813 
rs13475883  rs6239834  rs3162895  rs6212146  rs3022802 rs13475899 rs13475900 
28.4641674 30.8427150 31.1526514 31.2751278 31.3428706 31.8493556 31.8518088 
rs13475906  rs3022803 rs13475907 rs13475909 rs13475912  rs6209698 rs13475929 
32.2967145 32.3074644 32.3683291 32.8001894 33.6026526 36.5341646 37.6881435 
rs13475931 rs13475933 rs13475934  rs4222476 rs13475939  rs8253473 rs13475941 
37.7429827 38.0416271 38.0430095 38.9647582 39.4116688 39.4192277 39.4871064 
rs13475944 rs13475947 rs13475948 rs13475950 rs13475951 rs13475954 rs13475955 
39.7672829 40.2599440 40.3380113 40.3417592 40.3439501 41.1407252 41.2887176 
rs13475963 rs13475966 rs13475967 rs13475970 rs13475960  rs6250696 rs13475973 
42.4744416 42.5667702 42.9736574 43.1427994 43.5985261 43.5992946 43.6014053 
 rs3691187 rs13475988 rs13475991 rs13476023 rs13476024  rs3684654  rs6274257 
44.6237384 45.7855528 46.0180221 47.8579278 47.8600317 48.2423958 48.9612178 
rs13476045 rs13476049  rs6319405 rs13476050 rs13476051 rs13476054 rs13476057 
49.2018340 49.3701384 49.4261039 49.4275718 49.4323558 49.4972616 49.5031830 
rs13476059 rs13476060 rs13476062 rs13476066 rs13476067  rs6259837 rs13476080 
49.5084008 49.5113545 49.6085043 49.6644819 50.1779477 50.8256056 51.0328603 
 rs6302966 rs13476085 rs13476089  rs3717360 rs13476090  rs6248251 rs13476091 
51.3659553 51.6974451 52.3869798 52.3903517 52.3936241 52.4228715 52.5787388 
 rs3088725  rs3022832  rs4222577 rs13476100  rs6263067  rs8256168  rs6327099 
53.4044231 53.4129004 53.4189013 54.3267003 54.4193890 55.1459517 55.3274320 
rs13476111 rs13476119  rs8236484  rs8270838  rs8236489 rs13476129 rs13476134 
55.9050491 56.8936305 56.9852502 57.1870637 58.0248893 58.7605079 59.5401544 
rs13476135 rs13476137 rs13476138 rs13476140 rs13476148  rs6202860 rs13476158 
59.5426193 59.6023794 60.3355828 60.3439598 61.1791787 61.9905512 61.9930265 
rs13476163 rs13476177 rs13476178 rs13476183 rs13476184  rs6194543 rs13476196 
62.0039607 62.6243588 62.6269118 63.8101331 64.0856907 66.4047817 66.7425394 
rs13476201  rs3685700  rs3022846 rs13476210 rs13476214 rs13459163  rs4222816 
67.2638714 68.7230251 68.7246243 69.1209547 70.1550813 75.5548371 75.5593190 
 rs4222820  rs3090340  rs8245949 rs13476242 rs13476251 rs13476254  rs6383012 
75.5593202 75.5637846 76.7508053 79.0157673 79.7644000 79.8248805 85.3173344 
rs13476279  rs6348421 rs13476290 rs13476300 rs13476302 rs13476304  rs3669814 
86.7653503 88.2128991 89.0565541 94.6215368 94.8227821 94.8269227 95.5413280 
rs13501301 rs13476316 
96.0784002 96.9960494 

$`2`
rs13476325 rs13476327 rs13476328 rs13476330  rs3695983 rs13476334 rs13476337 
  1.329379   1.760872   1.839732   1.950151   1.954566   2.265170   3.619681 
rs13476342  rs3696091 rs13476348  rs3681847 rs13476358 rs13476424 rs13476427 
  4.113919   5.308337   6.708720   7.382269  10.053730  21.102901  21.547272 
rs13476432 rs13476433 rs13476438 rs13476440 rs13476445 rs13476446 rs13476448 
 22.457323  22.458829  22.472305  22.502827  22.703755  22.706606  22.809631 
rs13476449 rs13476451 rs13476452 rs13476456  rs6333344 rs13476459  rs6288325 
 22.811905  22.816989  22.917366  23.496527  23.499872  23.505774  23.764498 
rs13476470 rs13476482  rs6203572 rs13476485 rs13476501 rs13476502 rs13476503 
 25.671278  26.670288  26.702994  27.001404  29.014864  29.017621  29.019099 
rs13476525 rs13476526 rs13476530  rs3660779 rs13476533 rs13476534  rs4223189 
 31.635258  31.882499  32.606343  33.182340  33.622217  33.881289  33.886227 
 rs6205317 rs13476536  rs3709716 rs13476538 rs13476540  rs6314726 rs13476543 
 34.125908  34.182229  34.777129  35.011099  35.196852  35.859566  36.045536 
rs13476544  rs6222797 rs13476546 rs13476553 rs13476554  rs8263229 rs13476565 
 36.199822  36.656504  36.783097  37.624894  37.628034  38.509801  39.813261 
rs13476566  rs3670631 rs13476583 rs13476586 rs13476592 rs13476595 rs13476645 
 40.314536  41.215814  43.089828  43.274068  45.681933  46.026997  50.383519 
rs13476655 rs13476660  rs6378047 rs13476661 rs13476666 rs13476667 rs13476669 
 50.804234  51.471634  51.474692  51.557954  52.083779  52.257580  52.543537 
rs13476672  rs3724460  rs3143279 rs13476687  rs6278009 rs13476693  rs3022892 
 52.594599  53.100014  54.071835  54.256293  54.721524  55.005598  55.109510 
rs13476700 rs13476702 rs13476703 rs13476705  rs4223406  rs3090608 rs13476739 
 55.226885  55.233738  55.400903  55.678374  55.949739  57.841429  58.879499 
rs13476747 rs13476754  rs6257970 rs13476758  rs3022901 rs13476769  rs4223486 
 59.218037  59.231937  59.742629  60.304581  61.226450  61.739777  61.741096 
rs13476774  rs4223511 rs13476778 rs13476783  rs6234650 rs13476784 rs13476786 
 61.912021  62.703328  63.294674  63.904805  64.118754  64.219897  64.224107 
rs13476788  rs3682725 rs13476801 rs13476803 rs13476816  rs6170159  rs3022909 
 64.987248  65.971246  66.931607  67.071161  69.496928  69.584633  69.689891 
rs13476819 rs13476822 rs13476823  rs3716380  rs6332517 rs13476826 rs13476827 
 69.754708  71.191302  71.405084  71.816078  71.959179  72.019045  72.087910 
rs13476830 rs13476831 rs13476832  rs4223605  rs3022932 rs13476872  rs3692409 
 72.364799  72.368511  72.371502  76.460822  76.621947  78.531760  79.528125 
 rs3726342 rs13476882  rs3673613  rs8260429  rs8275858 rs13476892 rs13476894 
 80.337464  81.259831  82.295892  82.969957  83.669130  84.046818  84.787594 
 rs3024096 rs13476907 rs13476910 rs13476918  rs3673248 rs13476928 rs13476934 
 87.473424  87.498263  87.955210  92.672743  94.072805  96.238835  98.601444 
rs13476935 rs13476936 
 99.057281  99.455950 

$`3`
 rs8246929 rs13476972 rs13476973 rs13476974  rs3693600  rs3716641 rs13476982 
  1.080189   1.347091   1.470204   1.476022   1.553792   1.553874   1.563144 
rs13476983 rs13476984  rs6398851 rs13476992 rs13476993 rs13476994 rs13476995 
  1.569943   1.573624   1.576176   2.809096   3.100703   3.300776   3.380060 
 rs6168642 rs13477004  rs6410894 rs13477005 rs13477006 rs13477007 rs13477019 
  3.552793   4.223627   4.585754   4.691027   4.791228   5.088769   7.798872 
rs13477027 rs13477030 rs13477031 rs13477043  rs3677091  rs4223883  rs6371982 
  9.282482   9.532585   9.606423  13.791018  14.074175  14.939077  15.157908 
rs13477059  rs6246699 rs13477061  rs3022959  rs3140619 rs13477063 rs13477066 
 15.918550  15.923947  15.992721  16.802271  16.802276  16.963472  16.970994 
 rs6257041 rs13477071 rs13477078 rs13477079 rs13477108 rs13477109 rs13477110 
 17.319985  17.755564  18.693494  18.744862  21.080665  21.413925  22.016411 
 rs3685081 rs13477111  rs4223969 rs13477121 rs13477134 rs13477137 rs13477138 
 22.640148  22.641512  25.049880  25.052026  27.282182  27.588349  27.591228 
rs13477140  rs4223979 rs13477142 rs13477143 rs13477144  rs6326655 rs13477153 
 27.693007  27.696543  27.696853  27.700581  28.329883  28.503475  28.562688 
rs13477154 rs13477155 rs13477156 rs13477160 rs13477168 rs13477169 rs13477174 
 28.613511  28.615742  28.694417  28.702636  29.332336  29.349436  30.386145 
rs13477176 rs13477179  rs3023594 rs13477183 rs13477187 rs13477188 rs13477189 
 30.734338  31.007845  31.105126  31.112997  31.165690  31.167829  31.170044 
 rs6351657  rs6289734 rs13477194 rs13477201  rs8273664 rs13477207 rs13477209 
 31.173493  31.174332  32.184066  32.341020  32.342736  32.914386  32.973448 
rs13477210 rs13477216 rs13477217 rs13477220 rs13477224  rs6196421 rs13477225 
 33.009762  33.484516  33.517478  34.047286  34.586747  34.589541  34.590880 
rs13477228 rs13477230 rs13477233  rs6183172 rs13477236  rs4224039 rs13477248 
 35.026131  35.324560  35.893062  36.115534  36.317250  36.816125  37.403199 
rs13477252  rs6224931 rs13477268 rs13477263 rs13477271 rs13477279 rs13477282 
 37.639920  38.791882  39.157460  39.160076  39.405423  41.031706  41.504565 
rs13477283 rs13477287 rs13477288 rs13477291  rs3159432  rs3159588  rs3162061 
 41.544820  41.692166  41.775032  42.169443  42.366378  42.367700  42.369817 
rs13477294  rs4224158 rs13477318 rs13477369 rs13477373 rs13477374 rs13477379 
 42.957412  44.537215  45.502206  50.866600  51.429778  51.504763  52.411466 
 rs3022967 rs13477381 rs13477384 rs13477385 rs13477396 rs13477397 rs13477409 
 52.430405  53.262919  53.368879  53.372015  55.570112  55.709649  59.858613 
rs13477428 rs13477429 rs13477433 rs13477434 rs13477438 rs13477439  rs3698409 
 61.510633  61.670155  61.958144  62.047669  62.546905  62.688373  62.743012 
rs13477442 rs13477444 rs13477448  rs6385968 rs13477450 rs13477451 rs13477456 
 63.177344  63.279576  63.722373  63.935680  63.936910  63.938561  65.312917 
rs13477457 rs13477458 rs13477459 rs13477460 rs13477461 rs13477467 rs13477472 
 65.314857  65.317267  65.318925  65.320875  65.322208  67.548098  69.605885 
rs13477529 rs13477473  rs3022971 rs13477478 rs13477482 rs13477485  rs3724805 
 69.608224  69.609224  70.843094  70.907058  71.149729  74.456309  74.457760 
rs13477486 rs13477488 rs13477489 rs13477492 rs13477495  rs4224567 rs13477497 
 74.739366  74.745908  74.832780  75.275285  75.481388  75.596886  75.684561 
rs13477500 rs13477511 rs13477526 
 76.429887  79.219312  80.957067 

$`4`
rs13477532 rs13477541 rs13477542 rs13477543 rs13477544 rs13477549 rs13477550 
 0.3991292  1.2560214  1.3917929  1.5928861  1.6978705  2.0515920  2.2621396 
rs13477552 rs13477553 rs13477554 rs13477556 rs13477558 rs13477559 rs13477561 
 2.4352574  2.5251773  2.7585040  3.0530491  3.2150744  3.4585919  3.4626618 
 rs6403469  rs4137687 rs13477566 rs13477579 rs13477580 rs13477585 rs13477586 
 3.7029303  3.9430172  3.9449684  5.4077682  5.4078174  5.4306980  5.5055395 
 rs3709372 rs13477589 rs13477592 rs13459074 rs13477595 rs13477596  rs6163246 
 5.6750754  5.6767158  5.9597782  6.1852266  7.1180568  7.1205590  7.1217001 
 rs6266094 rs13477605 rs13477606 rs13477611  rs4224426 rs13477616 rs13477619 
 7.1234642  8.6068427  8.6549771  8.7910377  8.9901964  9.5194951  9.7572785 
rs13477623  rs3679673 rs13477625 rs13477626 rs13477627 rs13477628 rs13477631 
11.4021807 11.4058442 11.4073390 11.4114893 11.8595138 11.9573244 12.6677290 
rs13477642 rs13477643  rs3703981 rs13477659 rs13477660 rs13477662 rs13477664 
15.2447672 15.2467497 15.6652248 18.3775143 18.5369150 18.6189365 18.6215940 
rs13477715 rs13477717  rs6239799 rs13477722 rs13477724 rs13477725 rs13477730 
27.9120152 28.1256700 29.9945260 30.0009463 30.0051069 30.1318185 30.6953939 
 rs4224515 rs13477732 rs13477733  rs6184584 rs13477739  rs6254381 rs13476245 
30.7119837 30.8926924 31.1603646 31.6555204 32.3233798 32.4243902 32.4320470 
rs13477743  rs8276754 rs13477745 rs13477746 rs13477759  rs3654162 rs13477785 
32.4370370 32.5459374 32.7020530 32.9891931 33.4688108 34.3363190 35.7680418 
 rs3717147 rs13477793 rs13477794 rs13477795  rs6323325 rs13477814 rs13477821 
35.7771234 35.7913993 35.9203487 36.0115374 37.1562808 38.6853237 39.1848617 
rs13477823 rs13477824 rs13477835  rs6361351 rs13477854 rs13477860 rs13477863 
39.1885557 39.5402592 40.6564169 40.7350605 42.5809515 43.7049055 44.1279677 
rs13477865  rs6401724 rs13477919 rs13477888 rs13477907 rs13477908 rs13477910 
44.1349533 46.3675418 47.9452479 47.9493262 49.9169893 49.9768431 49.9868206 
rs13477912 rs13477914 rs13477915 rs13477916 rs13477920  rs4224727 rs13477931 
50.0962204 50.1002715 50.1016595 50.1100592 50.1276921 51.6805005 51.8938626 
 rs6173859 rs13477938 rs13477940  rs6365760 rs13477947 rs13477954 rs13477966 
52.0518071 53.1881228 53.5004864 53.9066510 55.3233700 55.9983634 59.7633519 
rs13477976 rs13477980  rs6355453 rs13477983 rs13477985 rs13459077 rs13477991 
61.7882160 61.9716108 63.0252404 63.0303080 63.0330779 64.5217020 64.9138820 
rs13477999 rs13478021 rs13459080  rs3023026  rs3711383 rs13478027 rs13478028 
66.3759417 70.2574630 71.8991146 72.3118040 73.8896348 74.1570955 76.0540649 
rs13478031  rs3023027 rs13478039  rs3668420 rs13478048 rs13478050 rs13478054 
76.3810043 77.0004907 77.1047625 77.2968965 79.8090093 79.9367987 81.0633814 

$`5`
rs13478092  rs4225033 rs13478097  rs6245801 rs13478103 rs13478104 rs13478106 
 0.5148540  0.5601752  1.4578186  1.6459862  1.6508495  1.6521104  1.7168465 
 rs8277477 rs13478110 rs13478111 rs13481337 rs13481338 rs13481339 rs13481340 
 1.7415309  1.9911696  1.9937269  2.2408610  2.2660419  2.2898547  2.2972730 
 rs3669635 rs13481349 rs13478115 rs13478119 rs13478123 rs13478127  rs6256752 
 4.5647660  4.6978589  4.9630837  6.4001473  6.9122371  7.6593491  8.1438747 
rs13478135  rs6295831 rs13478138  rs3023036 rs13478139 rs13478141 rs13478142 
 8.3149190  8.3964728  8.3989242  8.4007825  8.4033521  8.6865375  8.6885014 
rs13478143  rs3090417 rs13478145 rs13478146 rs13478148  rs6322061  rs6309009 
 8.8977080  9.1431643 10.2406114 10.4115566 10.6617649 11.1715598 11.2373953 
rs13478154 rs13478160 rs13478161 rs13478162 rs13478167 rs13459084 rs13478170 
11.3740589 12.8940905 12.9954134 13.2099571 15.1865866 15.5507310 15.5578335 
rs13459083 rs13478174  rs6196732 rs13478175 rs13478178 rs13478179 rs13459085 
15.6270239 15.7081356 15.8280039 16.1317236 16.1364419 16.1481976 16.2543350 
rs13478181 rs13478182 rs13478184 rs13478185  rs6346598 rs13478204 rs13478205 
16.3646745 16.4752220 17.4437033 17.4455556 17.9704569 20.1679510 20.7422551 
rs13478207 rs13478208 rs13478209 rs13478221 rs13478222 rs13478241 rs13478244 
20.8151842 20.8793971 20.8804978 22.5458012 22.5484334 24.8718483 25.1324122 
rs13478245 rs13478252  rs6237983 rs13478254 rs13478255  rs6263715 rs13478257 
25.1347167 25.3480011 25.6909676 26.7835663 27.2330524 27.3672129 27.9932349 
rs13478261 rs13478262 rs13478263 rs13478264 rs13478265 rs13478268 rs13478271 
28.8161178 28.8414537 28.8728759 29.0427959 29.1322946 29.2399880 29.9682578 
rs13478276 rs13478279 rs13478280 rs13478282 rs13478283 rs13478285  rs3684754 
29.9813803 30.1772663 30.2076157 30.2106527 30.2140050 30.3329785 34.3825206 
rs13478311 rs13478320 rs13478367 rs13478368  rs6167407 rs13478386 rs13478387 
34.5022391 35.7707846 41.4115154 41.4669600 42.7614533 42.7644662 43.2203787 
rs13478391  rs6232866 rs13478393 rs13478395 rs13478396 rs13478402 rs13478403 
43.2305343 44.4209462 45.0370667 45.7058881 45.7104424 45.7465608 45.7527138 
 rs4225380  rs6370004 rs13478410 rs13478422 rs13478424 rs13478426 rs13478430 
46.0876570 46.2329778 46.2430427 47.0473590 47.0514238 47.4628441 48.7119128 
rs13478442 rs13478443 rs13478445 rs13478452 rs13478454 rs13478455 rs13478458 
49.9285133 49.9304123 49.9361657 51.7267769 51.7292199 51.7325687 51.9937088 
rs13478536 rs13478542 rs13478546 rs13478547  rs8265976  rs4225536 rs13478570 
70.3211281 73.2280807 74.3892894 74.3918300 75.2624965 75.4633566 79.8560921 
rs13478574  rs3023062  rs4138683  rs3668534 rs13478593 rs13459095 
81.1622904 83.1121830 85.1593074 86.9092920 87.4283034 87.4895523 

$`6`
rs13459096 rs13478602  rs6172481 rs13478606 rs13478610 rs13478667 rs13478694 
 0.2471448  0.2509530  0.2595298  0.2688093  0.7422945 10.1262711 11.5504738 
rs13478696 rs13478697 rs13478698 rs13478705  rs6163979 rs13478712  rs3024195 
12.4131565 12.5038661 12.5050255 13.5351920 13.5357026 13.9919190 13.9955592 
rs13478720 rs13478721 rs13478727  rs3696518 rs13478728  rs6214044 rs13478729 
16.3415722 16.3443428 19.9127889 19.9793499 20.1028494 20.1045540 20.1064628 
 rs6191358 rs13478730  rs6168553  rs3023068  rs3023069 rs13478756 rs13478757 
20.1802678 20.2900546 22.8125569 23.8816735 23.8817887 24.2350441 24.2361490 
rs13478758 rs13478759 rs13478761 rs13478762 rs13478768 rs13478802 rs13478819 
24.2385415 24.3438438 24.4218302 25.3208587 25.8918938 29.0540752 30.9696170 
rs13478827 rs13478829 rs13478830 rs13478833 rs13478834 rs13478843 rs13478845 
31.0889580 31.0940937 31.1459713 32.0501760 32.0754991 33.3794730 33.4494137 
rs13478852 rs13478856 rs13478858 rs13478863 rs13478864 rs13478865 rs13459097 
33.6822587 34.0640838 34.3991792 34.5451575 34.5462936 34.7572579 35.5171427 
rs13478871  rs4221995  rs4226063 rs13478876 rs13478897 rs13478901 rs13478904 
36.0826404 36.2882063 36.3528729 36.5675727 40.5244122 41.3426572 41.7358722 
rs13478905 rs13478906 rs13478913 rs13478915 rs13478916 rs13478917 rs13478919 
42.4408121 42.8756819 43.2898350 43.3444653 43.4229502 43.5628361 43.6888274 
rs13478927 rs13478932 rs13478933 rs13478934 rs13478935 rs13478936 rs13478941 
45.1993835 45.9039023 46.0695345 46.0719750 46.1226745 46.1662896 46.3825629 
rs13478948 rs13478950  rs6344812  rs3675615 rs13478966 rs13478971 rs13478977 
47.4986393 47.5033268 47.5109362 48.2958455 48.8980754 50.2242802 50.8859454 
rs13478978 rs13478983 rs13478985 rs13478990 rs13479058 rs13479059  rs4226339 
51.2879832 51.6367222 51.6425091 52.3587248 64.7346868 64.8288933 65.1405567 
rs13479070 rs13479071  rs3024135  rs8262456 rs13479085 rs13479086  rs3089737 
66.4917745 67.2476353 68.2936963 70.9494089 75.6833494 75.7983622 75.8754886 
rs13479087 rs13479088  rs3023105  rs3090435  rs3023103 rs13479099 rs13479091 
76.1284889 76.1348324 76.1353741 76.1370335 76.1407613 76.1430157 76.6498568 
 rs3023840 rs13479093  rs6155595 rs13479094 
76.6518167 76.9251446 76.9270209 76.9292552 
```

Next we use `calc_genoprob()` to calculate the QTL genotype probabilities.


``` r
probs <- calc_genoprob(cross      = cross, 
                       map        = cross$gmap, 
                       error_prob = 0.002)
```

The argument `error_prob` supplies an assumed genotyping error probability of 
0.002. If a value for `error_prob` is not supplied, the default probability is 
0.0001. 

Recall that the result of `calc_genoprob`, `probs`, is a list of 
three-dimensional arrays (one per chromosome). 


``` r
names(probs)
```

``` output
 [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
[16] "16" "17" "18" "19" "X" 
```

Each three-dimensional array of probabilities is arranged as individuals 
$\times$ genotypes $\times$ positions. Have a look at the names of each of the 
three dimensions for chromosome 19.


``` r
dimnames(probs$`19`)
```

``` output
[[1]]
  [1] "Mouse3051" "Mouse3551" "Mouse3430" "Mouse3476" "Mouse3414" "Mouse3145"
  [7] "Mouse3656" "Mouse3242" "Mouse3427" "Mouse3527" "Mouse3281" "Mouse3405"
 [13] "Mouse3530" "Mouse3477" "Mouse3498" "Mouse3526" "Mouse3284" "Mouse3160"
 [19] "Mouse3655" "Mouse3615" "Mouse3491" "Mouse3603" "Mouse3191" "Mouse3130"
 [25] "Mouse3330" "Mouse3199" "Mouse3614" "Mouse3577" "Mouse3081" "Mouse3204"
 [31] "Mouse3513" "Mouse3219" "Mouse3331" "Mouse3301" "Mouse3503" "Mouse3083"
 [37] "Mouse3568" "Mouse3189" "Mouse3287" "Mouse3131" "Mouse3311" "Mouse3357"
 [43] "Mouse3149" "Mouse3256" "Mouse3644" "Mouse3217" "Mouse3212" "Mouse3082"
 [49] "Mouse3156" "Mouse3535" "Mouse3481" "Mouse3123" "Mouse3359" "Mouse3555"
 [55] "Mouse3597" "Mouse3624" "Mouse3314" "Mouse3128" "Mouse3531" "Mouse3295"
 [61] "Mouse3231" "Mouse3496" "Mouse3438" "Mouse3183" "Mouse3052" "Mouse3237"
 [67] "Mouse3462" "Mouse3293" "Mouse3543" "Mouse3276" "Mouse3200" "Mouse3502"
 [73] "Mouse3171" "Mouse3364" "Mouse3524" "Mouse3334" "Mouse3355" "Mouse3254"
 [79] "Mouse3358" "Mouse3468" "Mouse3192" "Mouse3214" "Mouse3536" "Mouse3606"
 [85] "Mouse3226" "Mouse3393" "Mouse3415" "Mouse3266" "Mouse3648" "Mouse3224"
 [91] "Mouse3474" "Mouse3381" "Mouse3138" "Mouse3660" "Mouse3616" "Mouse3425"
 [97] "Mouse3554" "Mouse3196" "Mouse3528" "Mouse3312" "Mouse3045" "Mouse3585"
[103] "Mouse3471" "Mouse3308" "Mouse3628" "Mouse3429" "Mouse3324" "Mouse3124"
[109] "Mouse3291" "Mouse3452" "Mouse3373" "Mouse3367" "Mouse3579" "Mouse3647"
[115] "Mouse3169" "Mouse3335" "Mouse3122" "Mouse3635" "Mouse3154" "Mouse3484"
[121] "Mouse3652" "Mouse3612" "Mouse3668" "Mouse3233" "Mouse3175" "Mouse3306"
[127] "Mouse3046" "Mouse3663" "Mouse3165" "Mouse3519" "Mouse3592" "Mouse3127"
[133] "Mouse3184" "Mouse3650" "Mouse3599" "Mouse3494" "Mouse3605" "Mouse3505"
[139] "Mouse3573" "Mouse3561" "Mouse3489" "Mouse3480" "Mouse3186" "Mouse3421"
[145] "Mouse3607" "Mouse3346" "Mouse3375" "Mouse3633" "Mouse3589" "Mouse3094"
[151] "Mouse3611" "Mouse3307" "Mouse3133" "Mouse3152" "Mouse3518" "Mouse3209"
[157] "Mouse3056" "Mouse3320" "Mouse3365" "Mouse3313" "Mouse3441" "Mouse3339"
[163] "Mouse3352" "Mouse3159" "Mouse3619" "Mouse3238" "Mouse3203" "Mouse3137"
[169] "Mouse3509" "Mouse3289" "Mouse3054" "Mouse3432" "Mouse3487" "Mouse3179"
[175] "Mouse3572" "Mouse3285" "Mouse3466" "Mouse3252" "Mouse3517" "Mouse3546"
[181] "Mouse3185" "Mouse3665" "Mouse3537" "Mouse3096" "Mouse3600" "Mouse3349"
[187] "Mouse3098" "Mouse3275" "Mouse3667" "Mouse3342" "Mouse3333" "Mouse3300"
[193] "Mouse3244" "Mouse3478" "Mouse3560" "Mouse3501" "Mouse3315" "Mouse3440"
[199] "Mouse3669" "Mouse3486" "Mouse3632" "Mouse3319" "Mouse3453" "Mouse3172"
[205] "Mouse3121" "Mouse3590" "Mouse3215" "Mouse3447" "Mouse3618" "Mouse3340"
[211] "Mouse3047" "Mouse3666" "Mouse3516" "Mouse3225" "Mouse3167" "Mouse3207"
[217] "Mouse3631" "Mouse3444" "Mouse3168" "Mouse3298" "Mouse3602" "Mouse3309"
[223] "Mouse3416" "Mouse3260" "Mouse3146" "Mouse3374" "Mouse3144" "Mouse3485"
[229] "Mouse3610" "Mouse3348" "Mouse3500" "Mouse3613" "Mouse3253" "Mouse3384"
[235] "Mouse3664" "Mouse3206" "Mouse3426" "Mouse3332" "Mouse3210" "Mouse3283"
[241] "Mouse3670" "Mouse3120" "Mouse3274" "Mouse3461" "Mouse3202" "Mouse3472"
[247] "Mouse3437" "Mouse3434" "Mouse3593" "Mouse3055" "Mouse3234" "Mouse3422"
[253] "Mouse3571" "Mouse3236" "Mouse3049" "Mouse3350" "Mouse3249" "Mouse3326"
[259] "Mouse3134" "Mouse3143" "Mouse3493" "Mouse3361" "Mouse3636" "Mouse3436"
[265] "Mouse3510" "Mouse3117" "Mouse3601" "Mouse3303" "Mouse3497" "Mouse3544"
[271] "Mouse3463" "Mouse3118" "Mouse3354" "Mouse3162" "Mouse3464" "Mouse3181"
[277] "Mouse3188" "Mouse3356" "Mouse3521" "Mouse3591" "Mouse3241" "Mouse3467"
[283] "Mouse3469" "Mouse3262" "Mouse3643" "Mouse3548" "Mouse3372" "Mouse3542"
[289] "Mouse3563" "Mouse3583" "Mouse3584" "Mouse3208" "Mouse3661" "Mouse3659"
[295] "Mouse3195" "Mouse3459" "Mouse3653" "Mouse3649" "Mouse3382" "Mouse3180"
[301] "Mouse3386" "Mouse3084" "Mouse3205" "Mouse3299" "Mouse3515" "Mouse3540"
[307] "Mouse3255" "Mouse3177" "Mouse3523" "Mouse3366" "Mouse3567" "Mouse3557"
[313] "Mouse3114" "Mouse3623" "Mouse3419" "Mouse3580" "Mouse3271" "Mouse3385"
[319] "Mouse3492" "Mouse3119" "Mouse3232" "Mouse3598" "Mouse3150" "Mouse3310"
[325] "Mouse3164" "Mouse3587" "Mouse3050" "Mouse3627" "Mouse3506" "Mouse3413"
[331] "Mouse3435" "Mouse3151" "Mouse3112" "Mouse3630" "Mouse3646" "Mouse3223"
[337] "Mouse3187" "Mouse3263" "Mouse3637" "Mouse3662" "Mouse3508" "Mouse3550"
[343] "Mouse3125" "Mouse3545" "Mouse3570" "Mouse3641" "Mouse3136" "Mouse3626"
[349] "Mouse3166" "Mouse3269" "Mouse3529" "Mouse3218" "Mouse3625" "Mouse3448"
[355] "Mouse3378" "Mouse3227" "Mouse3651" "Mouse3182" "Mouse3304" "Mouse3617"
[361] "Mouse3141" "Mouse3552" "Mouse3479" "Mouse3658" "Mouse3539" "Mouse3190"
[367] "Mouse3093" "Mouse3097" "Mouse3126" "Mouse3170" "Mouse3229" "Mouse3520"
[373] "Mouse3582" "Mouse3351" "Mouse3129" "Mouse3153" "Mouse3450" "Mouse3113"
[379] "Mouse3586" "Mouse3549" "Mouse3538" "Mouse3201" "Mouse3556" "Mouse3247"
[385] "Mouse3455" "Mouse3176" "Mouse3344" "Mouse3343" "Mouse3439" "Mouse3629"
[391] "Mouse3286" "Mouse3216" "Mouse3588" "Mouse3488" "Mouse3221" "Mouse3142"
[397] "Mouse3428" "Mouse3111" "Mouse3353" "Mouse3211" "Mouse3569" "Mouse3280"
[403] "Mouse3325" "Mouse3368" "Mouse3553" "Mouse3245" "Mouse3228" "Mouse3135"
[409] "Mouse3622" "Mouse3095" "Mouse3369" "Mouse3609" "Mouse3410" "Mouse3302"
[415] "Mouse3594" "Mouse3483" "Mouse3197" "Mouse3336" "Mouse3507" "Mouse3305"
[421] "Mouse3532" "Mouse3250" "Mouse3194" "Mouse3449" "Mouse3178" "Mouse3198"
[427] "Mouse3620" "Mouse3596" "Mouse3638" "Mouse3222" "Mouse3147" "Mouse3163"
[433] "Mouse3273" "Mouse3473" "Mouse3578" "Mouse3465" "Mouse3279" "Mouse3558"
[439] "Mouse3443" "Mouse3490" "Mouse3460" "Mouse3248" "Mouse3243" "Mouse3431"
[445] "Mouse3564" "Mouse3347" "Mouse3565" "Mouse3525" "Mouse3574" "Mouse3329"
[451] "Mouse3140" "Mouse3257" "Mouse3328" "Mouse3193" "Mouse3132" "Mouse3220"
[457] "Mouse3235" "Mouse3499" "Mouse3246" "Mouse3270" "Mouse3608" "Mouse3442"
[463] "Mouse3157" "Mouse3642" "Mouse3566" "Mouse3139" "Mouse3282" "Mouse3053"
[469] "Mouse3454" "Mouse3363" "Mouse3213" "Mouse3654" "Mouse3514" "Mouse3341"
[475] "Mouse3401" "Mouse3388" "Mouse3604" "Mouse3161" "Mouse3451" "Mouse3634"
[481] "Mouse3482" "Mouse3559" "Mouse3645" "Mouse3264" "Mouse3155" "Mouse3251"
[487] "Mouse3297" "Mouse3541" "Mouse3158" "Mouse3294"

[[2]]
[1] "BB" "BR" "RR"

[[3]]
 [1] "rs4232073"  "rs13483548" "rs13483549" "rs13483550" "rs13483554"
 [6] "rs13483555" "rs3090321"  "rs3090137"  "rs6309315"  "rs13483577"
[11] "rs3090325"  "rs13483579" "rs13483584" "rs13483586" "rs13483587"
[16] "rs13483589" "rs13483592" "rs13483593" "rs6344448"  "rs13483594"
[21] "rs13483595" "rs3705022"  "rs13483609" "rs13483612" "rs13483648"
[26] "rs13483650" "rs13483654" "rs13483658" "rs13483660" "rs13483664"
[31] "rs13483666" "rs13483667" "rs13483670" "rs8275553"  "rs8275912" 
[36] "rs13483677" "rs13483679" "rs13483680" "rs13483681" "rs3660143" 
[41] "rs13483682" "rs13483683" "rs13483685" "rs13483686" "rs6355398" 
[46] "rs4222106"  "rs13483690" "rs13483693" "rs13483695" "rs13483699"
```

View the first three rows of genotype probabilities for the first genotyped 
marker on chromosome 19.


``` r
(probs$`19`)[1:5,,"rs4232073"]
```

``` output
                    BB           BR           RR
Mouse3051 1.317728e-11 1.235895e-07 9.999999e-01
Mouse3551 9.999840e-01 1.595361e-05 5.027172e-08
Mouse3430 1.317728e-11 1.235895e-07 9.999999e-01
Mouse3476 9.999999e-01 1.235895e-07 1.317728e-11
Mouse3414 6.179474e-08 9.999999e-01 6.179474e-08
```

We can also view the genotype probabilities using 
[plot_genoprob](https://github.com/rqtl/qtl2/blob/master/R/plot_genoprob.R). The 
arguments to this function specify:

1. `probs`: the genotype probabilities,
1. `map`: the marker map,
1. `ind`: the index of the individual to plot,
1. `chr`: the index of the chromosome to plot.


``` r
plot_genoprob(probs = probs, 
              map   = cross$pmap, 
              ind   = 1, 
              chr   = 19, 
              main  = rownames(probs[['19']])[1])
```

<img src="fig/calc-genoprob-rendered-plot_genoprob-1.png" style="display: block; margin: auto;" />

The coordinates along chromosome 19 are shown on the horizontal axis and the 
three genotypes are shown on the vertical axis. Higher genotype probabilities 
are plotted in darker shades. This mouse has a RR genotype on the proximal end 
of the chromosome and transitions to BR.

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1

<!-- DMG: Not sure about this. I'd rather plot more samples or chromosomes. -->

1). Load a second dataset from Arabidopsis recombinant inbred lines 
([Moore et al, Genetics, 2013](https://www.genetics.org/content/195/3/1077)) 
in a study of plant root response to gravity (gravitropism).  

`grav <- read_cross2(file = system.file('extdata', 'grav2.zip', package = 'qtl2'))`  
2). How many individuals were in the study? How many phenotypes? 
How many chromosomes?  
3). Insert pseudomarkers at 1cM intervals and save the results 
to an object called `gravmap`. Have a look at the first chromosome.  
4). Calculate genotype probabilities and save the results to an object 
called `gravpr`. View the genotypes for the first three markers and 
pseudomarkers on chromosome 1 for the first five individuals.   

:::::::::::::::::::::::: solution 

1). `grav <- read_cross2(file = system.file('extdata', 'grav2.zip', package = 'qtl2'))`  
2). `summary(grav)`   
3). `gravmap <- insert_pseudomarkers(map = grav$gmap, step = 1)`  
followed by `head(gravmap, n=1)`   
4). `gravpr  <- calc_genoprob(cross = grav, map = gravmap)` followed by  
`(gravpr$``1``)[1:5,,"PVV4"]`, `(gravpr$`1`)[1:5,,"c1.loc1"]`, and  
`(gravpr$`1`)[1:5,,"c1.loc2"]`

<!-- DMG: What about this? -->


``` r
m  = maxmarg(probs)
ph = guess_phase(cross, m)
plot_onegeno(ph, cross$pmap)
```

<img src="fig/calc-genoprob-rendered-unnamed-chunk-1-1.png" style="display: block; margin: auto;" />

:::::::::::::::::::::::::::::::::


## Challenge 3

Calculate genotype probabilities for a different data set from the 
[qtl2 data repository](https://github.com/rqtl/qtl2data), this one from a study 
of obesity and diabetes in a C57BL/6 (B6) × BTBR intercross.   
1). Create a new script in RStudio with File -> New File -> R Script.  
2)  Download the B6 x BTBR zip file from the 
[qtl2 data repository](https://github.com/rqtl/qtl2data) into an object 
called `b6btbr` by running this code:  
`b6btbr <- read_cross2(file = "https://raw.githubusercontent.com/rqtl/qtl2data/master/B6BTBR/b6btbr.zip")`  
3). View a summary of the `b6btbr` data. How many individuals? phenotypes? 
chromosomes? markers?  
4). View the genetic map for the `b6btbr` data.  
5). Insert pseudomarkers at 2 cM intervals. Assign the results to an object 
called `b6btbrmap`.    
6). Calculate genotype probabilities assuming a genotyping error probability 
of 0.001. Assign the results to an object called `b6btbrpr`.    
7). View the first several rows of genotype probabilities for any marker on 
chromosome 18.  

:::::::::::::::::::::::: solution 

1). Create a new script in RStudio with File -> New File -> R Script.  
2). `b6btbr <- read_cross2(file = "https://raw.githubusercontent.com/rqtl/qtl2data/master/B6BTBR/b6btbr.zip`  
3). `summary(b6btbr)` shows 544 individuals, 3 phenotypes, 20 chromosomes, 
2057 markers.  
4). `b6btbr$gmap`  
5). `b6btbrmap <- insert_pseudomarkers(map=b6btbr$gmap, step=2)`  
6). `b6btbrpr <- calc_genoprob(cross=b6btbr, map=b6btbrmap, error_prob=0.001)`  
7). `dimnames((b6btbrpr$`18`))` shows all marker names for chromosome 18. 
`head((b6btbrpr$`18`)[,,"c18.loc48"])` gives genotype probabilities for an 
example pseudomarker, while `head((b6btbrpr$`18`)[,,"rs6338896"])`  gives 
genotype probabilities for a genotyped marker.

:::::::::::::::::::::::::::::::::


:::::::::::::::::::::::::::::::::::::


**Parallel calculations (optional)** To speed up the calculations with large 
datasets on a multi-core machine, you can use the argument `cores`. With 
`cores=0`, the number of available cores will be detected via 
`parallel::detectCores()`. Otherwise, specify the number of cores as a positive 
integer.


``` r
probs <- calc_genoprob(cross = iron, map = map, error_prob = 0.002, cores = 4)
```

**Allele probabilities (optional)** The genome scan functions use genotype 
probabilities as well as a matrix of phenotypes. If you wished to perform a 
genome scan via an additive allele model, you would first convert the genotype 
probabilities to allele probabilities, using the function 
`genoprob_to_alleleprob()`.


``` r
apr <- genoprob_to_alleleprob(probs = probs)
```

The figure below shows genotype and allele probabilities for 3 samples. In the 
Diversity Outbred, there are 36 possible genotype states 
(AA, AB, AC, ..., BB, BC, BD, ..., CC, CD, CE, ..., DD, DE, DF, ..., EE,...) or 
8 + 7 + 6 + 5 + 4 + 3 + 2 + 1. The first SNP below has genotype BB. In the table 
describing alleles (8 state founder probabilities), the probability that this 
SNP has a B allele is 1. The 2nd SNP has genotype BH, so the allele table shows 
a probability of 0.5 for B and 0.5 for H. The third SNP is either BG or BH, and 
has a probability of 0.5 for each of these genotypes. The allele table shows a 
probability of 0.5 for allele B, and 0.25 for both G and H.

![](../fig/geno-to-allele-probs.png)
::::::::::::::::::::::::::::::::::::: keypoints 

- The first step in QTL analysis is to calculate genotype probabilities.
- Calculate genotype probabilities between genotyped markers with 
`calc_genoprob()`.

::::::::::::::::::::::::::::::::::::::::::::::::
