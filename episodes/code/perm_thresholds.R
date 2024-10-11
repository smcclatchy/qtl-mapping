# Script to run permutations of 10, 100, and 1000 to compare the variance
# of the threshold estimate with differing numbers of permutations.

library(tidyverse)
library(ggbeeswarm)
library(qtl2)

cross <- read_cross2(file = 'https://thejacksonlaboratory.box.com/shared/static/svw7ivp5hhmd7vb8fy26tc53h7r85wez.zip')
probs <- calc_genoprob(cross = cross, map = cross$gmap, error_prob = 0.002)
cross$covar$Sex <- factor(cross$covar$Sex)
addcovar <- model.matrix(~Sex, data = cross$covar)[,-1, drop = FALSE]

np = c(10, 100, 1000)

df = data.frame(n_perm = rep(0, 300), lod = rep(0, 300))

index = 1

for(n_perm in np) {

  print(n_perm)

  for(i in 1:100) {

    print(i)
    p = scan1perm(probs, cross$pheno[,1,drop = F], addcovar = addcovar,
                  n_perm = n_perm)
    df[index,] = c(n_perm, 	quantile(p, probs = 0.95))

    index = index + 1

  } # for(i)

} # for(n_perm)

png('C:/Users/c-dgatti/Documents/classes/JAX/qtl-mapping/episodes/fig/permutation_simulations.png',
    width = 1200, height = 800, res = 128)
p = df %>%
  mutate(n_perm = as.character(n_perm)) %>%
  ggplot(aes(n_perm, lod)) +
    geom_violin(draw_quantiles = 0.5, linewidth = 1.25) +
    geom_beeswarm(alpha = 0.1, size = 3) +
    labs(title = 'Variance of Significance Threshold Estimates',
         x = 'Number of Permutations', y = 'LOD Threshold') +
    theme(text = element_text(size = 24))
print(p)
dev.off()





