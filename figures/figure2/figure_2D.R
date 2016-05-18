library(MASS)
my.heat.colors <- function(x) { rev(heat.colors(x, alpha=1)) }

# P(SC | SC)
# HomoHomo
SC = read.table("prior_SC", h=T)
proba = read.table("probas_SI_IM_AM_SC_SC", h=F)

threshold = max(SC[,4])

pNA = which(proba[, 1] < 0.8 & proba[, 2] < 0.8 & proba[, 3] < 0.8 & proba[, 4] < 0.8)
tmpHomoHomo = kde2d(SC[pNA, 4], (SC[pNA, 4]-SC[pNA, 5])/SC[pNA, 4], n=100)

dev.new(width = 6.5, height = 6)
filled.contour(tmpHomoHomo, color.palette=colorRampPalette(c('white','blue','yellow','red','darkred')), ylim=c(0,1), xlim=c(0, threshold), plot.title = title(xlab = expression(T["split"]), ylab = "Relative isolation before SC"))

dev.print(pdf, "figure2D.pdf", bg="white")
dev.off()

