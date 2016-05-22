# estimation of simulated data

library(MASS)
stat = read.table("statistics_IM.txt", h=T)
param = read.table("parameters_IM.txt", h=T)

estimatesN1 = NULL
estimatesN2 = NULL
estimatesNanc = NULL
estimatesTsplit = NULL
estimatesPropNtrlM1 = NULL
estimatesPropNtrlM2 = NULL
estimatesM1 = NULL
estimatesM2 = NULL

for(i in 1:1000){
        tmp = read.table(paste("posterior/Posterior_IM_HETERO_v1_migBinom_", i, sep=""), h=F)
	estimatesN1 = c(estimatesN1, median(tmp[, 1]))
	estimatesN2 = c(estimatesN2, median(tmp[, 2]))
	estimatesNanc = c(estimatesNanc, median(tmp[, 3]))
	estimatesTsplit = c(estimatesTsplit, median(tmp[, 4]))
	estimatesM1 = c(estimatesM1, median(tmp[, 5]))
	estimatesM2 = c(estimatesM2, median(tmp[, 6]))
	estimatesPropNtrlM1 = c(estimatesPropNtrlM1, median(tmp[, 7]))
        estimatesPropNtrlM2 = c(estimatesPropNtrlM2, median(tmp[, 8]))
}

estimatesN = c(estimatesN1, estimatesN2)
estimatesM = c(estimatesM1, estimatesM2)
estimatesAlpha = c(estimatesPropNtrlM1, estimatesPropNtrlM2) / 541 

realValuesN = c(param[, 1], param[,2])
realValuesNanc = param[, 3]
realValuesTsplit = param[, 4]
realValuesM = c(param[, 5], param[, 6])
realValuesAlpha = c(param[, 7], param[, 8]) / 541 

# estimation made from popPhylABC data
obs = read.table("estimation_popPhyl.txt",h=T)

# N1 and N2
tmp = kde2d(realValuesN*0.1, estimatesN*0.1, n = 50)

dev.new(width=10, height=5)
par(las=1)
layout(matrix(1:3, ncol=3), width=c(1/5, 1, 1))
image.scale(tmp, horiz=FALSE, col=rev(heat.colors(100)), n=10)
image(tmp, xlab = "current effective size (real values)", ylab = "current effective size (estimates)", cex.lab=1.3, col=rev(heat.colors(100)), n=10, cex.axis=1.2)
points(realValuesN*0.1, estimatesN*0.1)
abline(a=0,b=1); abline(a=0, b=0.5, lty=2); abline(a=0, b=2, lty=2, lwd=1.25)

plot(rep(log10(obs$netDiv), 2), c(obs$medianN1, obs$medianN2)*0.1, xlab="net synonymous divergence", ylab="current effective size", cex=1.5, xlim=c(log10(1e-5), log10(1)), cex.axis=1.2, xaxt="n", cex.lab=1.3)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.2)

dev.print(pdf, "figureS10.pdf", bg="white")
dev.off()

# Nanc
tmp = kde2d(realValuesNanc*0.1, estimatesNanc*0.1, n = 50)

dev.new(width=10, height=5)

par(las=1)
layout(matrix(1:3, ncol=3), width=c(1/5, 1, 1))
image.scale(tmp, horiz=FALSE, col=rev(heat.colors(100)), n=10)
image(tmp, xlab = "ancestral effective size (real values)", ylab = "ancestral effective size (estimates)", cex.lab=1.3, col=rev(heat.colors(100)), n=10, cex.axis=1.2)
points(realValuesNanc*0.1, estimatesNanc*0.1)
abline(a=0,b=1); abline(a=0, b=0.5, lty=2); abline(a=0, b=2, lty=2, lwd=1.25)

plot(log10(obs$netDiv), obs$medianNanc*0.1, xlab="net synonymous divergence", ylab="Nanc", cex=1.5, xlim=c(log10(1e-5), log10(1)), cex.axis=1.2, xaxt="n", cex.lab=1.3)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.2)

dev.print(pdf, "figureS11.pdf", bg="white")
dev.off()

# Tsplit 
tmp = kde2d(realValuesTsplit*0.4, estimatesTsplit*0.4, n = 50)

dev.new(width=10, height=5)

par(las=1)
layout(matrix(1:3, ncol=3), width=c(1/5, 1, 1))
image.scale(tmp, horiz=FALSE, col=rev(heat.colors(100)), n=10)
image(tmp, xlab = "Tsplit (real values)", ylab = "Tsplit (estimates)", cex.lab=1.3, col=rev(heat.colors(100)), n=10, cex.axis=1.2)
points(realValuesTsplit*0.4, estimatesTsplit*0.4)
abline(a=0,b=1); abline(a=0, b=0.5, lty=2); abline(a=0, b=2, lty=2, lwd=1.25)

plot(log10(obs$netDiv), obs$medianT*0.4, xlab="net synonymous divergence", ylab="Tsplit", cex=1.5, xlim=c(log10(1e-5), log10(1)), cex.axis=1.2, xaxt="n", cex.lab=1.3)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.2)

dev.print(pdf, "figureS12.pdf", bg="white")
dev.off()


# alpha
tmp = kde2d(realValuesAlpha, estimatesAlpha, n = 50)

dev.new(width=10, height=5)

par(las=1)
layout(matrix(1:3, ncol=3), width=c(1/5, 1, 1))
image.scale(tmp, horiz=FALSE, col=rev(heat.colors(100)), n=10)
image(tmp, xlab = "proportion of introgressed genome (real values)", ylab = "proportion of introgressed genome (estimates)", cex.lab=1.3, col=rev(heat.colors(100)), xlim=c(0,1), ylim=c(0,1), n=10, cex.axis=1.2)
points(realValuesAlpha, estimatesAlpha)
abline(a=0,b=1); abline(a=0, b=0.5, lty=2); abline(a=0, b=2, lty=2)

plot(log10(obs$netDiv), (obs$median_alphaM1+obs$median_alphaM2)/(2*obs$nLocus),ylim=c(0,1), xlab="net synonymous divergence", ylab="proportion of introgressed genome", cex=1.5, xlim=c(log10(1e-5), log10(1)), cex.axis=1.2, xaxt="n", cex.lab=1.3)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.2)

dev.print(pdf, "figureS8.pdf", bg="white")
dev.off()

# M1 and M2
tmp = kde2d(realValuesM/4, estimatesM/4, n = 50)

dev.new(width=10, height=5)

par(las=1)
layout(matrix(1:3, ncol=3), width=c(1/5, 1, 1))
image.scale(tmp, horiz=FALSE, col=rev(heat.colors(100)), n=10)
image(tmp, xlab = "Nm (real values)", ylab = "Nm (estimates)", cex.lab=1.3, col=rev(heat.colors(100)), n=10, cex.axis=1.2,)
points(realValuesM/4, estimatesM/4)
abline(a=0,b=1); abline(a=0, b=0.5, lty=2); abline(a=0, b=2, lty=2, lwd=1.25)

plot(rep(log10(obs$netDiv), 2), c(obs$medianM1, obs$medianM2)/4, xlab="net synonymous divergence", ylab="Nm", cex=1.5, xlim=c(log10(1e-5), log10(1)), cex.axis=1.2, xaxt="n", cex.lab=1.3)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.2)

dev.print(pdf, "figureS9.pdf", bg="white")
dev.off()
   

