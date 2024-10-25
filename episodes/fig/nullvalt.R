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

png(filename = "fig/nullvalt.png", width = 600, height = 600)
par(mai = c(1, 1, 1, 1.5))
plot_pxg(g, iron$pheno[,"liver"], ylab="Phenotype", 
         main = "Null and Alternative\nHypotheses",
         seg_col = c("#ffa07a", "#2c7bb6"), 
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
mtext(expression(+ beta), side = 4, line = 0.3, las = 1,
      at=(plusbeta))
mtext(expression(- beta), side = 4, line = 0.3, las = 1,
      at=(minusbeta))
mtext(expression(mu), side = 4, line = 0.3, las = 1,
      at=(mu))

# fake regression lines
segments(0.95, mu, 3.05, mu,
         col = "#2c7bb6", lwd = 3)
segments(0.95, minusbeta, 3.05, plusbeta,
         col = "#ffa07a", lwd = 3)
text(x = 1.5, y = 103, 
     labels = "NULL\n",
     cex=0.8, col = "#2c7bb6")
text(x = 1.51, y = 99, 
     labels = "y = 94.6 + error",
     cex=0.8)
text(x = 1.55, y = 78, 
     labels = paste0("ALTERNATIVE\n"),
     cex=0.8, col = "#cd8162")
text(x = 1.55, y = 72, 
     labels = paste0("y = 94.6 + ",
                     expression(beta),
                     "*X\n + error"),
     cex=0.8)

# a single data point showing error
points(2, 157, lwd = 3, col = "#ba55d3")
segments(x0 = 2, y0 = 95, y1 = 155,
         lwd=2, lty = 1, col = "#ba55d3")
text(x = 1.7, y = 160, labels = "error")
segments(x0 = 1.8, y0 = 159, x1 = 2, y1 = 151, 
         lwd=1, lty = 2, col = "#2c7bb6")
dev.off()
