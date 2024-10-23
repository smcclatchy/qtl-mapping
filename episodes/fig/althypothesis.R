library(qtl2)
iron <- read_cross2(file = system.file("extdata",
                                       "iron.zip",
                                       package = "qtl2"))
iron$alleles <- c("B", "R")

map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)
g <- maxmarg(pr, map, chr=2, pos=56.8, return_char=TRUE)

mu <- 94.6
plusbeta <- 110
minusbeta <- 79.3

png(filename = "fig/althypothesis.png", width = 600, height = 600)
par(mai = c(1, 1, 1, 1.5))
plot_pxg(g, iron$pheno[,"liver"], ylab="Phenotype", 
         main = "Alternative Hypothesis",
         seg_col = "#ffa07a",
         seg_lwd = 4,
         sort = FALSE,
         jitter = 0.1)

# mu and beta marks on right y axis
segments(3.4, mu, 3.5, mu,
         col = "#2c7bb6", lwd = 3)
segments(3.4, plusbeta, 3.5, plusbeta,
         col = "#ffa07a", lwd = 3)
segments(3.4, minusbeta, 3.5, minusbeta,
         col = "#ffa07a", lwd = 3)

# fake regression line
segments(0.95, minusbeta, 3.05, plusbeta,
         col = "#ffa07a", lwd = 3)

# a single data point showing error squared
points(3, 155, lwd = 3, col = "#ba55d3")
segments(x0 = 3, y0 = plusbeta, y1 = 155,
         lwd=2, lty = 1)
segments(x0 = 2.5, y0 = plusbeta, y1 = 155,
         lwd=2, lty = 1)
segments(x0 = 2.5, y0 = plusbeta, x1 = 3,
         lwd=2, lty = 1)
segments(x0 = 2.5, y0 = 155, x1 = 3,
         lwd=2, lty = 1)
text(x = 2.5, y = 170, labels = "error\nsquared")
segments(x0 = 2.5, y0 = 160, x1 = 2.7, y1 = 140, 
         lwd=1, lty = 2)
dev.off()
