#! /usr/bin/env Rscript
png(filename="figure2.png", width = 35, height = 25, units = "cm", bg = "white", res = 900)
#pdf(file="fig_tmp.pdf", width = 13.7795, height = 9.84252, bg = "white")
facteur = 9 # to correct for the wrong mutation rate per locus used before launching simulations.
modelsMig = c("summary_3_IM_MHetero_NHetero_robustness.txt", "summary_3_IM_MHetero_NHomo_robustness.txt", "summary_3_IM_MHomo_NHetero_robustness.txt", "summary_3_IM_MHomo_NHomo_robustness.txt", "summary_3_PAN_NHetero_robustness.txt", "summary_3_PAN_NHomo_robustness.txt", "summary_3_SC_MHetero_NHetero_robustness.txt", "summary_3_SC_MHetero_NHomo_robustness.txt", "summary_3_SC_MHomo_NHetero_robustness.txt", "summary_3_SC_MHomo_NHomo_robustness.txt")
modelsNoMig = c("summary_3_AM_MHetero_NHetero_robustness.txt", "summary_3_AM_MHetero_NHomo_robustness.txt", "summary_3_AM_MHomo_NHetero_robustness.txt", "summary_3_AM_MHomo_NHomo_robustness.txt", "summary_3_SI_NHetero_robustness.txt", "summary_3_SI_NHomo_robustness.txt")

par(mfrow=c(2,2))

# migration models
cnt = 0
tmp_mig = NULL
for(i in modelsMig){
	cnt = cnt + 1
	x = read.table(paste("../../../revisions/analyse2/", i, sep=""), h=T)
	i = strsplit(i, ".", fixed=T)[[1]][1]
	i = paste(strsplit(i, "_", fixed=T)[[1]][-1], collapse= " ")
	tmp_mig = rbind(tmp_mig, cbind(x$netdivAB_avg*facteur, x$pMigration, x$robustness)) # divergence, pMigration, robustness
	if(cnt == 1){
		#plot(log10(x$netdivAB_avg), x$pMigration, ylim = c(0,1), ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25, main="IM, SC, PAN")
		plot(0, 0, ylim = c(0,1), ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25, main="IM, SC, PAN")
		par(las=3)
		mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
		par(las=1)
		axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
		axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
	}else{
	#	points(log10(x$netdivAB_avg), x$pMigration, cex.axis=1.25, col="white")
	}
}
trueP = which(tmp_mig[,2] > 0.5 & tmp_mig[,3] >= 0.95)
falseP = which(tmp_mig[,2] < 0.5 & tmp_mig[,3] >= 0.95)
ambig = which(tmp_mig[,3] < 0.95)
points(log10(tmp_mig[,1][trueP]), tmp_mig[,2][trueP], pch = 16, col = "green")
points(log10(tmp_mig[,1][falseP]), tmp_mig[,2][falseP], pch = 16, col = "red")
points(log10(tmp_mig[,1][ambig]), tmp_mig[,2][ambig], pch = 16, col = "grey")


# no-migration models
cnt = 0
tmp_noMig = NULL
for(i in modelsNoMig){
	cnt = cnt + 1
	x = read.table(paste("../../../revisions/analyse2/", i, sep=""), h=T)
	tmp_noMig = rbind(tmp_noMig, cbind(x$netdivAB_avg*facteur, x$pMigration, x$robustness))
	
	if(cnt == 1){
		#plot(log10(x$netdivAB_avg), x$pMigration, ylim = c(0,1), ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25, main= "AM, SI")
		plot(0, 0, ylim = c(0,1), ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25, main= "AM, SI")
		par(las=3)
		mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
		par(las=1)
		axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
		axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
	}else{
	#	points(log10(x$netdivAB_avg), x$pMigration, cex.axis=1.25, col="white")
	}
}
trueP = which(tmp_noMig[,2] < 0.5 & tmp_noMig[,3] >= 0.95)
falseP = which(tmp_noMig[,2] > 0.5 & tmp_noMig[,3] >= 0.95)
ambig = which(tmp_noMig[,3] < 0.95)
points(log10(tmp_noMig[,1][trueP]), tmp_noMig[,2][trueP], pch = 16, col = "green")
points(log10(tmp_noMig[,1][falseP]), tmp_noMig[,2][falseP], pch = 16, col = "red")
points(log10(tmp_noMig[,1][ambig]), tmp_noMig[,2][ambig], pch = 16, col = "grey")


colnames(tmp_noMig) = c("divergence", "pMigration", "robustness")
colnames(tmp_mig) = c("divergence", "pMigration", "robustness")
# true positive rate, false positive rate, rate ambiguity
# mig
trueP = NULL
falseP = NULL
ambig = NULL
bins = c(1e-4, 1e-3, 1e-2, 1e-1, 1)
for(i in 1:length(bins)){
	if(i==1){
		tmp = tmp_mig[which(tmp_mig[,1] < bins[i]),]
	}else{
		tmp = tmp_mig[which(tmp_mig[,1] > bins[i-1] & tmp_mig[,1] < bins[i]),]
	}
	trueP = c(trueP, length(which(tmp[,2] > 0.5 & tmp[,3] >= 0.95))/nrow(tmp))
	falseP = c(falseP, length(which(tmp[,2] < 0.5 & tmp[,3] >= 0.95))/nrow(tmp))
	ambig = c(ambig, length(which(tmp[,3] < 0.95))/nrow(tmp))
}

#dev.new(width=13, height=7)
par(las=1)
res = rbind(trueP, ambig, falseP)
names = c("Da<10⁻⁴", "10⁻⁴<Da<10⁻³", "10⁻³<Da<10⁻²", "10⁻²<Da<10⁻¹", "0.1<Da")
barplot(res, beside = T, col=c("green", "grey", "red"), names.arg = names, ylab = "rate", cex.lab = 1.25, cex.axis = 1.25, cex = 0.75, ylim = c(0, 1))
abline(h=0.05, lwd = 2, col="red", lty=2)

# true positive rate, false positive rate, rate ambiguity
# no-mig
trueP = NULL
falseP = NULL
ambig = NULL
bins = c(1e-3, 1e-2, 1e-1, 1)
for(i in 1:length(bins)){
	if(i==1){
		tmp = tmp_noMig[which(tmp_noMig[,1] < bins[i]),]
	}else{
		tmp = tmp_noMig[which(tmp_noMig[,1] > bins[i-1] & tmp_noMig[,1] < bins[i]),]
	}
	trueP = c(trueP, length(which(tmp[,2] < 0.5 & tmp[,3] >= 0.95))/nrow(tmp))
	falseP = c(falseP, length(which(tmp[,2] > 0.5 & tmp[,3] >= 0.95))/nrow(tmp))
	ambig = c(ambig, length(which(tmp[,3] < 0.95))/nrow(tmp))
}

#dev.new(width=13, height=7)
par(las=1)
res = rbind(trueP, ambig, falseP)
names = c("Da<10⁻³", "10⁻³<Da<10⁻²", "10⁻²<Da<10⁻¹", "0.1<Da")
barplot(res, beside = T, col=c("green", "grey", "red"), names.arg = names, ylab = "rate", cex.lab = 1.25, cex.axis = 1.25, cex = 0.75, ylim = c(0, 1))
abline(h=0.05, lwd = 2, col="red", lty=2)
dev.off()

