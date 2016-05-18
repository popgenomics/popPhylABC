library(shape)
HH = read.table("migHetero_NHetero.txt", h=T)
st = read.table("species_statues.txt",h=T)

paste_tmp = function(a, b){
	return(paste(as.character(a), as.character(b),sep="_"))
}

name_st = NULL
for(i in 1:nrow(st)){
	name_st = c(name_st, paste_tmp(st[i,1], st[i,2]))
}
st = cbind(st, name_st)
st = st[,c(4,3)]

name_HH = NULL
for(i in 1:nrow(st)){
	name_HH = c(name_HH, paste_tmp(HH[i,1], HH[i,2]))
}
HH = cbind(HH, name_HH)

HH = merge(HH, st, by.x=ncol(HH), by.y=1)[,-1] 

bilanHH = HH$pIMHomov1 + HH$pIMHomov2 + HH$pSCHomov1 + HH$pSCHomov2 + HH$pIMHeterov1 + HH$pIMHeterov2 + HH$pSCHeterov1 + HH$pSCHeterov2 + HH$pPanmixiav1 + HH$pPanmixiav2


tmp = rep(NA, length(bilanHH))
tmp[which(bilanHH>=0.8)] = "Migration"
tmp[which(bilanHH<.2)] = "NoMigration"
bilanHH = tmp

pmig_HH = HH$pIMHomov1 + HH$pIMHomov2 + HH$pSCHomov1 + HH$pSCHomov2 + HH$pIMHeterov1 + HH$pIMHeterov2 + HH$pSCHeterov1 + HH$pSCHeterov2 + HH$pPanmixiav1 + HH$pPanmixiav2

# plot MHetero_NHetero
dev.new(width=9, height=5)
par(las=2, mar=c(5,4,4,12))
plot(log10(HH$netdivAB_avg), pmig_HH, ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)

div_NA = log10(sort(HH[is.na(bilanHH),]$netdivAB_avg))
div_Mig = log10(sort(HH[which(bilanHH=="Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(HH[which(bilanHH=="NoMigration"),]$netdivAB_avg))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(0.2, 0.8), lty=2, lwd=1.25)

par(las=1)
# formes des points en fonction du statu d'espces
formes = rep(0, nrow(HH))
formes[which(HH$status == 0)] = 22 
formes[which(HH$status == 0.5)] = 22
formes[which(HH$status == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(HH))
heteroM = apply(HH[,grep("Hetero", colnames(HH))], FUN="sum", MARGIN=1)

heteroN = apply(HH[,grep("v2", colnames(HH))], FUN="sum", MARGIN=1)

#homoM_homoN: species with homogeneity for both introgression and Ne
couleurs[which(pmig_HH>= 0.8 & heteroM >= 0.8)] = "orange"
couleurs[which(pmig_HH>= 0.8 & heteroM <= 0.8)] = "turquoise"
couleurs[which(pmig_HH<= 0.2)] = "red"
couleurs[which(pmig_HH> 0.2 & pmig_HH<0.8)] = grey(0.25)

points(log10(HH$netdivAB_avg), pmig_HH, pch=formes, col="black", bg=couleurs, cex=1.8)

dev.print(pdf, "figure4.pdf", bg="white")
dev.off()

