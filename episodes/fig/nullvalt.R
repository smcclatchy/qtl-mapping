library(qtl2)
iron <- read_cross2(file = system.file("extdata",
                                       "iron.zip",
                                       package = "qtl2"))
g <- maxmarg(pr, map, chr=2, pos=56.8, return_char=TRUE)
c2eff <- scan1coef(pr[,"2"], iron$pheno[,"liver"])

png(filename = "nullvalt.png", width = 600, height = 600)
par(mai = c(1, 1, 1, 1.5))
plot_pxg(g, iron$pheno[,"liver"], ylab="Phenotype", 
         main = "Null and Alternative\nHypotheses",
         seg_col = c("#ffa07a", "#2c7bb6"), 
         seg_lwd = 4)
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
       col = "#ffa07a")
abline(a = (c2eff["D2Mit17",4]), 
       b = c2eff["D2Mit17", 2],
       col = "#2c7bb6")
text(x = 2.5, y = 102, 
     labels = "NULL\n",
     cex=0.8, col = "#2c7bb6")
text(x = 2.5, y = 98, 
     labels = "y = 94.6 + error",
     cex=0.8)
text(x = 2.5, y = 78, 
     labels = paste0("ALTERNATIVE\n"),
     cex=0.8, col = "#ffa07a")
text(x = 2.5, y = 72, 
     labels = paste0("y = 94.6 + ",
                     expression(beta),
                     "*X\n + error"),
     cex=0.8)
points(2, 157, lwd = 3)
segments(x0 = 2, y0 = 95, y1 = 155,
         lwd=2, lty = 1)
text(x = 1.7, y = 160, labels = "ERROR", cex=0.8)
segments(x0 = 1.8, y0 = 159, x1 = 2, y1 = 151, 
         lwd=1, lty = 2, col = "#2c7bb6")
dev.off()
