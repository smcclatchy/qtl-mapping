library(qtl2)
iron <- read_cross2(file = system.file("extdata",
                                       "iron.zip",
                                       package = "qtl2"))
map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)
g <- maxmarg(pr, map, chr=2, pos=56.8, return_char=TRUE)
c2eff <- scan1coef(pr[,"2"], iron$pheno[,"liver"])

png(filename = "nullvalt.png", width = 600, height = 600)
par(mai = c(1, 1, 1, 1.5))
plot_pxg(g, iron$pheno[,"liver"], ylab="Phenotype", 
         main = "Line of best fit",
         seg_col = c("#2c7bb6"), 
         seg_lwd = 4,
         jitter = 0.05)
segments(3.4, 94.6, 3.5, 94.6,
         col = "#2c7bb6", lwd = 3)
segments(3.4, 120, 3.5, 120,
         col = "#ffa07a", lwd = 3)
segments(3.4, 72, 3.5, 72,
         col = "#ffa07a", lwd = 3)
mtext(expression(+ beta), side = 4, line = 0.3, las = 1,
      at=(120))
mtext(expression(- beta), side = 4, line = 0.3, las = 1,
      at=(72))
mtext(expression(mu), side = 4, line = 0.3, las = 1,
      at=(c2eff["D2Mit17",4]))
abline(a = (c2eff["D2Mit17",4] + (2 * c2eff["D2Mit17",3])), 
       b = c2eff["D2Mit17", 1],
       col = "#ffa07a",
       lwd=1.5)
abline(a = (c2eff["D2Mit17",4]), 
       b = c2eff["D2Mit17", 2],
       col = "#2c7bb6",
       lwd=1.5)
text(x = 2.5, y = 102,
     labels = "NULL\n",
     cex=0.8, col = "#2c7bb6")
text(x = 2.5, y = 78,
     labels = paste0("ALTERNATIVE\n"),
     cex=0.8, col = "#ffa07a")
points(1, 155, lwd = 3)
segments(x0 = 1, y0 = 110, y1 = 154,
         lwd=2, lty = 1)
text(x = 1.3, y = 165, labels = "error or\nresidual", cex=0.8)
segments(x0 = 1.2, y0 = 164, x1 = 1.05, y1 = 157, 
         lwd=1, lty = 2)
dev.off()
