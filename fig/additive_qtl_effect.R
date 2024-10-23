# Make a sample plot of an additive and interactive QTL effect.
# We may have to edit this in Inkscape to get it right.
# Daniel Gatti

library(tidyverse)
library(ggbeeswarm)

base_dir = 'C:/Users/c-dgatti/Documents/classes/JAX/qtl-mapping/episodes'

# Possible genotypes.
gt = c('BB', 'BR', 'RR')

# Create a data frame with genotypes and Gaussian noise.
n = 500
data = data.frame(id    = paste0('M', 1:n),
                  sex   = sample(c('Female', 'Male'), size = n, 
                                 replace = TRUE),
                  geno  = sample(gt, size = n, replace = TRUE, 
                                 prob = c(0.25, 0.5, 0.25)),
                  pheno = rnorm(n, mean = 0, sd = 0.5))

# Add a sex effect.
sex_eff_size = 2
sex_eff      = as.numeric(factor(data$sex)) - 1
sex_eff      = sex_eff * sex_eff_size
data$pheno   = data$pheno + sex_eff

# Add an allele effect.
add_eff_size = 1
add_eff      = as.numeric(factor(data$geno)) - 2
add_eff      = add_eff * add_eff_size
data$pheno   = data$pheno + add_eff

# Get group means.
mod  = lm(pheno ~ sex + geno, data = data)
pred = unique(cbind(data$sex, data$geno, predict(mod)))
pred = data.frame(pred)
colnames(pred) = c('sex', 'geno', 'mean')
pred$mean = as.numeric(pred$mean)

# Make a plot.
p = data |>
  ggplot(aes(x = geno, y = pheno, color = sex)) +
    geom_beeswarm() +
    geom_smooth(aes(group = sex), method = 'lm') +
    geom_hline(aes(yintercept = mean), 
               filter(pred, sex == 'Female' & geno == 'BR'),
               color = 'red', linetype = 'dashed', linewidth = 1.25) +
    geom_hline(aes(yintercept = mean), 
               filter(pred, sex == 'Male' & geno == 'BR'),
               color = 'blue', linetype = 'dashed', linewidth = 1.25) +
    geom_hline(aes(yintercept = mean), 
               filter(pred, sex == 'Female' & geno == 'BB'),
               color = 'gray30', linetype = 'dashed', linewidth = 1.25) +
    geom_hline(aes(yintercept = mean), 
               filter(pred, sex == 'Male' & geno == 'RR'),
               color = 'gray30', linetype = 'dashed', linewidth = 1.25) +
    labs(title = 'Additive Sex Effect at QTL',
         x     = 'Genotype',
         y     = 'Phenotype',
         color = 'Sex') +
    theme(text = element_text(size = 24))

svg(file = file.path(base_dir, 'fig', 'additive_qtl_effect_w_sex.svg'),
    width = 9, height = 7)
print(p)
dev.off()

# Now remove the sex effect and plot again.
data2       = data
data2$pheno = data2$pheno - sex_eff

# Get group means.
mod2  = lm(pheno ~ geno, data = data2)
pred2 = unique(cbind(data2$geno, predict(mod2)))
pred2 = data.frame(pred2)
colnames(pred2) = c('geno', 'mean')
pred2$mean = as.numeric(pred2$mean)

p = data2 |>
  ggplot(aes(x = geno, y = pheno, color = sex)) +
    geom_beeswarm() +
    geom_smooth(aes(group = sex), method = 'lm') +
    geom_hline(aes(yintercept = mean), 
               filter(pred2, geno == 'BB'),
               color = 'gray30', linetype = 'dashed', linewidth = 1.25) +
    geom_hline(aes(yintercept = mean), 
               filter(pred2, geno == 'BR'),
               color = 'gray30', linetype = 'dashed', linewidth = 1.25) +
    geom_hline(aes(yintercept = mean), 
               filter(pred2, geno == 'RR'),
               color = 'gray30', linetype = 'dashed', linewidth = 1.25) +
    labs(title = 'Additive Effect at QTL',
         x     = 'Genotype',
         y     = 'Phenotype',
         color = 'Sex') +
    theme(text = element_text(size = 24))

svg(file = file.path(base_dir, 'fig', 'additive_qtl_effect_wo_sex.svg'),
    width = 9, height = 7)
print(p)
dev.off()

#############################################################################
# Interactive QTL

# Create a data frame with genotypes and Gaussian noise.
n = 500
data = data.frame(id    = paste0('M', 1:n),
                  sex   = sample(c('Female', 'Male'), size = n, 
                                 replace = TRUE),
                  geno  = sample(gt, size = n, replace = TRUE, 
                                 prob = c(0.25, 0.5, 0.25)),
                  pheno = rnorm(n, mean = 0, sd = 0.5))

# Create an allele effect.
add_eff_size = 1
add_eff      = as.numeric(factor(data$geno)) - 1
add_eff      = add_eff * add_eff_size

# Create a sex effect.
sex_eff_size = 2
sex_eff      = as.numeric(factor(data$sex)) - 1
sex_eff      = sex_eff * sex_eff_size

# Add a sex by allele effect.
data$pheno   = data$pheno + (sex_eff * add_eff)

# Get group means.
mod  = lm(pheno ~ sex + geno + sex:geno, data = data)
pred = unique(cbind(data$sex, data$geno, predict(mod)))
pred = data.frame(pred)
colnames(pred) = c('sex', 'geno', 'mean')
pred$mean = as.numeric(pred$mean)

p = data |>
  ggplot(aes(x = geno, y = pheno, color = sex)) +
    geom_beeswarm() +
    geom_smooth(aes(group = sex), method = 'lm') +
    geom_hline(aes(yintercept = mean), 
               filter(pred, sex == 'Male' & geno == 'BR'),
               color = 'blue', linetype = 'dashed', linewidth = 1.25) +
    geom_hline(aes(yintercept = mean), 
               filter(pred, sex == 'Male' & geno == 'RR'),
               color = 'blue', linetype = 'dashed', linewidth = 1.25) +
    geom_hline(aes(yintercept = mean), 
               filter(pred, sex == 'Female' & geno == 'BR'),
               color = 'red', linetype = 'dashed', linewidth = 1.25) +
    labs(title = 'Interactive Effect at QTL',
         x     = 'Genotype',
         y     = 'Phenotype',
         color = 'Sex') +
    theme(text = element_text(size = 24))

svg(file = file.path(base_dir, 'fig', 'interactive_qtl_effect_w_sex.svg'),
    width = 9, height = 7)
print(p)
dev.off()



