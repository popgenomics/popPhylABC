facteur = 10 # to correct for the wrong mutation rate per locus used before launching simulations.
modelsMig = c("summary_IM_MHetero_NHetero.txt", "summary_IM_MHetero_NHomo.txt", "summary_IM_MHomo_NHetero.txt", "summary_IM_MHomo_NHomo.txt", "summary_PAN_NHetero.txt", "summary_PAN_NHomo.txt", "summary_SC_MHetero_NHetero.txt", "summary_SC_MHetero_NHomo.txt", "summary_SC_MHomo_NHetero.txt", "summary_SC_MHomo_NHomo.txt")
modelsNoMig = c("summary_AM_MHetero_NHetero.txt", "summary_AM_MHetero_NHomo.txt", "summary_AM_MHomo_NHetero.txt", "summary_AM_MHomo_NHomo.txt", "summary_SI_NHetero.txt", "summary_SI_NHomo.txt")

par(mfrow=c(1,2))

# migration models
cnt = 0
for(i in modelsMig){
	cnt = cnt + 1
	x = read.table(i, h=T)
	i = strsplit(i, ".", fixed=T)[[1]][1]
	i = paste(strsplit(i, "_", fixed=T)[[1]][-1], collapse= " ")

	SI = grep("pSI", colnames(x))
	AM = grep("pAM", colnames(x))
	IM = grep("pIM", colnames(x))
	SC = grep("pSC", colnames(x))
	PAN = grep("pPAN", colnames(x))

	if(cnt == 1){
		plot(log10(x$netdivAB_avg*facteur), apply(x[, c(IM, SC, PAN)], MARGIN=1, FUN="sum"), ylim = c(0,1), ylab ="" , xlab = "net synonymous divergence", col = "black", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25, main="IM, SC, PAN")
		par(las=3)
		mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
		par(las=1)
		axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
		axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
	}else{
		points(log10(x$netdivAB_avg*facteur), apply(x[, c(IM, SC, PAN)], MARGIN=1, FUN="sum"), cex.axis=1.25)
	}
}


# no-migration models
cnt = 0
for(i in modelsNoMig){
	cnt = cnt + 1
	x = read.table(i, h=T)
	i = strsplit(i, ".", fixed=T)[[1]][1]
	i = paste(strsplit(i, "_", fixed=T)[[1]][-1], collapse= " ")

	SI = grep("pSI", colnames(x))
	AM = grep("pAM", colnames(x))
	IM = grep("pIM", colnames(x))
	SC = grep("pSC", colnames(x))
	PAN = grep("pPAN", colnames(x))

	if(cnt == 1){
		plot(log10(x$netdivAB_avg*facteur), apply(x[, c(IM, SC, PAN)], MARGIN=1, FUN="sum"), ylim = c(0,1), ylab ="" , xlab = "net synonymous divergence", col = "black", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25, main= "AM, SI")
		par(las=3)
		mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
		par(las=1)
		axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
		axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
	}else{
		points(log10(x$netdivAB_avg*facteur), apply(x[, c(IM, SC, PAN)], MARGIN=1, FUN="sum"), cex.axis=1.25)
	}
}

