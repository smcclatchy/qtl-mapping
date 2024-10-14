##########
# Run 1000 permutations 100 times.
# Include or exclude the kinship matrices and compare the
# signifcance thresholds.

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(qtl2))
cross <- read_cross2(file = 'https://thejacksonlaboratory.box.com/shared/static/svw7ivp5hhmd7vb8fy26tc53h7r85wez.zip')
probs <- calc_genoprob(cross = cross, map = cross$gmap, error_prob = 0.002)
kinship_loco  <- calc_kinship(probs = probs, type = 'loco')
# Covariates
cross$covar$Sex <- factor(cross$covar$Sex)
addcovar <- model.matrix(~Sex, data = cross$covar)[,-1, drop = FALSE]

n_perm  = 1000
n_sim   = 100
n_cores = 5

thr = data.frame(sim = 1:n_sim,
                 wo_kinship = 0,
                 w_kinship  = 0)

set.seed(2^30-1)

for(i in 1:n_sim) {

  t1 = proc.time()[3]
  print(i)

  perm_pheno = replicate(n = n_perm,
                         expr = sample(cross$pheno[,'log10_insulin_10wk',drop = FALSE]))
  dimnames(perm_pheno) = list(rownames(cross$pheno), 1:n_perm)

  # Without kinship.
  p = scan1(genoprobs = probs, 
            pheno     = perm_pheno,
            addcovar  = addcovar,
            cores     = n_cores)
  p = apply(p, 2, max)
  thr$wo_kinship[i] = quantile(p, probs = 0.95)

  # With kinship.
  p = scan1(genoprobs = probs, 
            pheno     = perm_pheno,
            kinship   = kinship_loco,
            addcovar  = addcovar,
            cores     = n_cores)
  p = apply(p, 2, max)
  thr$w_kinship[i] = quantile(p, probs = 0.95)

  t2 = proc.time()[3]
  print(t2 - t1)

} # for(i)

write.csv(thr, file = 'C:/Users/c-dgatti/Documents/classes/JAX/qtl-mapping/episodes/data/perm_thresh_kinship.csv', 
          quote = FALSE, row.names = FALSE)

png('C:/Users/c-dgatti/Documents/classes/JAX/qtl-mapping/episodes/fig/perm_kinship.png',
    width = 1000, height = 800, res = 128)
p = thr |>
  pivot_longer(cols = ends_with('kinship'),
               names_to = 'Kinship', values_to = 'LOD') |>
  mutate(Kinship = if_else(Kinship == 'wo_kinship', 'Without', 'With')) |>
  ggplot(aes(Kinship, LOD)) +
    geom_violin(draw_quantiles  = 0.5, linewidth = 1.25) +
    geom_beeswarm(alpha = 0.5, size = 2) +
    geom_line(aes(group = sim)) +
    labs(title = 'Significance Thresholds with/without Kinship',
         y = 'LOD Threshold (0.05)') +
    theme(text = element_text(size = 20))
print(p)
dev.off()



