# datasets are simulated under the SI model
# with theta = 5, T = seq(0.01, 1, 0.01) and no migration by definition.
# for a given value of T, 100 of datasets were simulated. These datasets
# were analyzed to estimate the parameters of the IM model

x = read.table("summary_IM1.txt", h=T)

comp = NULL
for(i in as.numeric(names(table(x$simulated_tsplit)))){
	tmp = NULL
	tmp = x[x$simulated_tsplit==i, ]
	nTruePos = length(which(tmp$best_model=="SI" & tmp$robustness>=0.95))
	nFalPos = length(which(tmp$best_model=="IM" & tmp$robustness>=0.95))
	nAmbig = length(which(tmp$robustness<0.95))
	comp = rbind(comp, c(i, nTruePos, nFalPos, nAmbig))
}
colnames(comp) = c("tsplit", "nTruePositives", "nFalsePositives", "nAmbigusCases")
comp = as.data.frame(comp)

# model comparison
par(las=1)
plot(comp$tsplit, comp$nFalsePositives/100, pch=16, cex=1.2, col = "red", xlab=expression(T["split"]), ylab = "rate", ylim = c(0,1), cex.lab = 1.5, cex.axis = 1.25)
points(comp$tsplit, comp$nTruePositives/100, pch=16, cex=1.2, col = "green")
points(comp$tsplit, comp$nAmbigusCases/100, pch=16, cex=1.2, col = "grey")

dev.print(pdf, "model_comparison.pdf", bg="white")
dev.off()


# estimation of parameters
pdf("parameters_estimation.pdf", bg="white", paper="a4",height = 11, width = 7 )
par(mfrow=c(4,1), mar = c(5, 5, 2, 1.25), las = 1)

boxplot(log10(x$netdivAB_avg)~x$simulated_tsplit, xlab = expression(T["split"]), ylab = "log10(simulated Da)", cex.lab = 1.5, cex.axis = 1.25)

boxplot(x$estimated_theta~x$simulated_tsplit, ylim = c(0, 10), xlab = expression(T["split"]), ylab = "estimated theta", cex.lab = 1.5, cex.axis = 1.25)
abline(h=5, col = "red", lty=2)

boxplot(x$estimated_tsplit~x$simulated_tsplit, ylim = c(0, 10), xlab = expression(T["split"]), ylab = expression(paste("estimated ", T["split"], paste="")), cex.lab = 1.5, cex.axis = 1.25)
lines(seq(0.01, 1, 0.01), col="red", lty=2)

boxplot(x$estimated_Nm~x$simulated_tsplit, ylim = c(0, 2.5), xlab = expression(T["split"]), ylab = "estimated N.m", cex.lab = 1.5, cex.axis = 1.25)
abline(h=0, col = "red", lty=2)
#dev.print(pdf, "parameters_estimation.pdf", bg="white")
dev.off()

# plot div and pIM
plot(log10(x$netdivAB_avg), x$pSI, col="white", xlab = "log10(simulated Da)", ylab = expression(P["isolation"]), cex.axis = 1.25, cex.lab = 1.25)
falseP = which(x$pIM>0.5 & x$robustness>=0.95)
trueP = which(x$pSI>0.5 & x$robustness>=0.95)
ambig = which(x$robustness<0.95)

points(log10(x$netdivAB_avg)[ambig], x$pSI[ambig], col="grey")
points(log10(x$netdivAB_avg)[falseP], x$pSI[falseP], col="red")
points(log10(x$netdivAB_avg)[trueP], x$pSI[trueP], col="green")

plot((x$FST_avg), x$pSI, col="white", xlab = "Fst", ylab = expression(P["isolation"]), cex.axis = 1.25, cex.lab = 1.25)
points((x$FST_avg)[ambig], x$pSI[ambig], col="grey")
points((x$FST_avg)[falseP], x$pSI[falseP], col="red")
points((x$FST_avg)[trueP], x$pSI[trueP], col="green")


