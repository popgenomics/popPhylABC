library(shape)
HH = read.table("migHetero_NHetero.txt", h=T)
Hh = read.table("migHetero_NHomo.txt", h=T)
hH = read.table("migHomo_NHetero.txt", h=T)
hh = read.table("migHomo_NHomo.txt", h=T)
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
name_Hh = NULL
name_hH = NULL
name_hh = NULL
for(i in 1:nrow(HH)){
	name_HH = c(name_HH, paste_tmp(HH[i,1], HH[i,2]))
	name_Hh = c(name_Hh, paste_tmp(Hh[i,1], Hh[i,2]))
	name_hH = c(name_hH, paste_tmp(hH[i,1], hH[i,2]))
	name_hh = c(name_hh, paste_tmp(hh[i,1], hh[i,2]))
}
HH = cbind(HH, name_HH)
HH = merge(HH, st, by.x=ncol(HH), by.y=1)[,-1] 
Hh = cbind(Hh, name_Hh)
Hh = merge(Hh, st, by.x=ncol(Hh), by.y=1)[,-1] 
hH = cbind(hH, name_hH)
hH = merge(hH, st, by.x=ncol(hH), by.y=1)[,-1] 
hh = cbind(hh, name_hh)
hh = merge(hh, st, by.x=ncol(hh), by.y=1)[,-1] 

bilanhh = hh$pSCHomov1 + hh$pIMHomov1 + hh$pPanmixiav1 
bilanHh = Hh$pIMHomov1 + Hh$pIMHeterov1 + Hh$pSCHomov1 + Hh$pSCHeterov1 + Hh$pPanmixiav1
bilanhH = hH$pIMHomov1 + hH$pIMHomov2 + hH$pSCHomov1 + hH$pSCHomov2 + hH$pPanmixiav1 + hH$pPanmixiav2
bilanHH = HH$pIMHomov1 + HH$pIMHomov2 + HH$pSCHomov1 + HH$pSCHomov2 + HH$pIMHeterov1 + HH$pIMHeterov2 + HH$pSCHeterov1 + HH$pSCHeterov2 + HH$pPanmixiav1 + HH$pPanmixiav2

tmp = rep(NA, length(bilanhh))
tmp[which(bilanhh>=0.8)] = "Migration"
tmp[which(bilanhh<.2)] = "NoMigration"
bilanhh = tmp

tmp = rep(NA, length(bilanHh))
tmp[which(bilanHh>=0.8)] = "Migration"
tmp[which(bilanHh<.2)] = "NoMigration"
bilanHh = tmp

tmp = rep(NA, length(bilanhH))
tmp[which(bilanhH>=0.8)] = "Migration"
tmp[which(bilanhH<.2)] = "NoMigration"
bilanhH = tmp

tmp = rep(NA, length(bilanHH))
tmp[which(bilanHH>=0.8)] = "Migration"
tmp[which(bilanHH<.2)] = "NoMigration"
bilanHH = tmp

pmig_hh = hh$pSCHomov1 + hh$pIMHomov1 + hh$pPanmixiav1
pmig_Hh = Hh$pIMHomov1 + Hh$pIMHeterov1 + Hh$pSCHomov1 + Hh$pSCHeterov1 + Hh$pPanmixiav1
pmig_hH = hH$pIMHomov1 + hH$pIMHomov2 + hH$pSCHomov1 + hH$pSCHomov2 + hH$pPanmixiav1 + hH$pPanmixiav2
pmig_HH = HH$pIMHomov1 + HH$pIMHomov2 + HH$pSCHomov1 + HH$pSCHomov2 + HH$pIMHeterov1 + HH$pIMHeterov2 + HH$pSCHeterov1 + HH$pSCHeterov2 + HH$pPanmixiav1 + HH$pPanmixiav2

# plot MHetero_NHetero
dev.new(width=10, height=8)
par(las=2, mar=c(5,4,2,2), mfrow=c(2,2))

# plot HH
plot(log10(HH$netdivAB_avg), pmig_HH, ylab ="" , xlab = "synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
mtext(side = 3, text = "hetero M + N", cex = 1.5)

#rect(-10, 0.2, 10, 0.8, col = rgb(1, 0, 0, 0.15))
div_NA = log10(sort(HH[is.na(bilanHH),]$netdivAB_avg))
div_Mig = log10(sort(HH[which(bilanHH=="Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(HH[which(bilanHH=="NoMigration"),]$netdivAB_avg))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
#xA2 = div_NA[2]
#xB2 = div_NA[length(div_NA)-1]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
#rect(xA2, -1, xB2, 2, col = rgb(0, 0, 0, 0.15), border=NA)
abline(h=c(0.2, 0.8), lty=2, lwd=1.25)

#points(log10(HH$netdivAB_avg), pmig_HH, col = rgb(0,0,0, 0.5))
par(las=1)
#mtext("MigHetero_NHetero", side = 3, line = 0.25)

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

# plot Hh
plot(log10(Hh$netdivAB_avg), pmig_Hh, ylab ="" , xlab = "synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)


#rect(-10, 0.2, 10, 0.8, col = rgb(1, 0, 0, 0.15))
div_NA = log10(sort(Hh[is.na(bilanHh),]$netdivAB_avg))
div_Mig = log10(sort(Hh[which(bilanHh=="Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(Hh[which(bilanHh=="NoMigration"),]$netdivAB_avg))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
#xA2 = div_NA[2]
#xB2 = div_NA[length(div_NA)-1]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
#rect(xA2, -1, xB2, 2, col = rgb(0, 0, 0, 0.15), border=NA)
abline(h=c(0.2, 0.8), lty=2, lwd=1.25)
#points(log10(Hh$netdivAB_avg), pmig_Hh, col = rgb(0,0,0, 0.5))
par(las=1)
#mtext("MigHetero_NHetero", side = 3, line = 0.25)

# formes des points en fonction du statu d'espces
formes = rep(0, nrow(Hh))
formes[which(Hh$status == 0)] = 22 
formes[which(Hh$status == 0.5)] = 22
formes[which(Hh$status == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(Hh))
heteroM = apply(Hh[,grep("Hetero", colnames(Hh))], FUN="sum", MARGIN=1)

heteroN = apply(Hh[,grep("v2", colnames(Hh))], FUN="sum", MARGIN=1)

#homoM_homoN: species with homogeneity for both introgression and Ne
couleurs[which(pmig_Hh>= 0.8 & heteroM >= 0.8)] = "orange"
couleurs[which(pmig_Hh>= 0.8 & heteroM <= 0.8)] = "turquoise"
couleurs[which(pmig_Hh<= 0.2)] = "red"
couleurs[which(pmig_Hh> 0.2 & pmig_Hh<0.8)] = grey(0.25)

points(log10(Hh$netdivAB_avg), pmig_Hh, pch=formes, col="black", bg=couleurs, cex=1.8)
mtext(side = 3, text = "hetero M + homo N", cex = 1.5)
for(i in 1:length(bilanHH)){
	testA = as.character(bilanHH[i])
	testB = as.character(bilanHh[i])
	if(is.na(testA)==T){testA="NA"}
	if(is.na(testB)==T){testB="NA"}
	if(testA != testB & testB!="NA"){
		Arrows(log10(HH$netdivAB_avg)[i], pmig_HH[i], log10(Hh$netdivAB_avg)[i], pmig_Hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(0,1,0, 0.75), lwd=1.75)
	}
	if(testA != testB & testB=="NA"){
		Arrows(log10(HH$netdivAB_avg)[i], pmig_HH[i], log10(Hh$netdivAB_avg)[i], pmig_Hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(1,0,0, 0.75), lwd=1.75)
	}
}


# plot hH
plot(log10(hH$netdivAB_avg), pmig_hH, ylab ="" , xlab = "synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)

#rect(-10, 0.2, 10, 0.8, col = rgb(1, 0, 0, 0.15))
div_NA = log10(sort(hH[is.na(bilanhH),]$netdivAB_avg))
div_Mig = log10(sort(hH[which(bilanhH=="Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(hH[which(bilanhH=="NoMigration"),]$netdivAB_avg))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
#xA2 = div_NA[2]
#xB2 = div_NA[length(div_NA)-1]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
#rect(xA2, -1, xB2, 2, col = rgb(0, 0, 0, 0.15), border=NA)
abline(h=c(0.2, 0.8), lty=2, lwd=1.25)
#points(log10(hH$netdivAB_avg), pmig_hH, col = rgb(0,0,0, 0.5))
par(las=1)
#mtext("MigHetero_NHetero", side = 3, line = 0.25)

#points(log10(hH$netdivAB_avg), pmig_hH, col = rgb(0,0,0, 0.5))
par(las=1)
#mtext("MigHetero_NHetero", side = 3, line = 0.25)

# formes des points en fonction du statu d'espces
formes = rep(0, nrow(hH))
formes[which(hH$status == 0)] = 22 
formes[which(hH$status == 0.5)] = 22
formes[which(hH$status == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(hH))
heteroM = apply(hH[,grep("Hetero", colnames(hH))], FUN="sum", MARGIN=1)

heteroN = apply(hH[,grep("v2", colnames(hH))], FUN="sum", MARGIN=1)

#homoM_homoN: species with homogeneity for both introgression and Ne
couleurs[which(pmig_hH>= 0.8 & heteroM >= 0.8)] = "orange"
couleurs[which(pmig_hH>= 0.8 & heteroM <= 0.8)] = "turquoise"
couleurs[which(pmig_hH<= 0.2)] = "red"
couleurs[which(pmig_hH> 0.2 & pmig_hH<0.8)] = grey(0.25)

points(log10(hH$netdivAB_avg), pmig_hH, pch=formes, col="black", bg=couleurs, cex=1.8)

mtext(side = 3, text = "homo M + hetero N", cex = 1.5)

for(i in 1:length(bilanHH)){
	testA = as.character(bilanHH[i])
	testB = as.character(bilanhH[i])
	if(is.na(testA)==T){testA="NA"}
	if(is.na(testB)==T){testB="NA"}
	if(testA != testB & testB!="NA"){
		Arrows(log10(HH$netdivAB_avg)[i], pmig_HH[i], log10(hH$netdivAB_avg)[i], pmig_hH[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(0,1,0, 0.75), lwd=1.75)
	}
	if(testA != testB & testB=="NA"){
		Arrows(log10(HH$netdivAB_avg)[i], pmig_HH[i], log10(hH$netdivAB_avg)[i], pmig_hH[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(1,0,0, 0.75), lwd=1.75)
	}
}


# plot hh
plot(log10(hh$netdivAB_avg), pmig_hh, ylab ="" , xlab = "synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)


#rect(-10, 0.2, 10, 0.8, col = rgb(1, 0, 0, 0.15))
div_NA = log10(sort(hh[is.na(bilanhh),]$netdivAB_avg))
div_Mig = log10(sort(hh[which(bilanhh=="Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(hh[which(bilanhh=="NoMigration"),]$netdivAB_avg))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
#xA2 = div_NA[2]
#xB2 = div_NA[length(div_NA)-1]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
#rect(xA2, -1, xB2, 2, col = rgb(0, 0, 0, 0.15), border=NA)
abline(h=c(0.2, 0.8), lty=2, lwd=1.25)
#points(log10(hh$netdivAB_avg), pmig_hh, col = rgb(0,0,0, 0.5))
par(las=1)
#mtext("MigHetero_NHetero", side = 3, line = 0.25)

# formes des points en fonction du statu d'espces
formes = rep(0, nrow(hh))
formes[which(hh$status == 0)] = 22 
formes[which(hh$status == 0.5)] = 22
formes[which(hh$status == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(hh))
heteroM = apply(hh[,grep("Hetero", colnames(hh))], FUN="sum", MARGIN=1)

heteroN = apply(hh[,grep("v2", colnames(hh))], FUN="sum", MARGIN=1)

#homoM_homoN: species with homogeneity for both introgression and Ne
couleurs[which(pmig_hh>= 0.8 & heteroM >= 0.8)] = "orange"
couleurs[which(pmig_hh>= 0.8 & heteroM <= 0.8)] = "turquoise"
couleurs[which(pmig_hh<= 0.2)] = "red"
couleurs[which(pmig_hh> 0.2 & pmig_hh<0.8)] = grey(0.25)

points(log10(hh$netdivAB_avg), pmig_hh, pch=formes, col="black", bg=couleurs, cex=1.8)

mtext(side = 3, text = "homo M + N", cex = 1.5)
for(i in 1:length(bilanHH)){
	testA = as.character(bilanHH[i])
	testB = as.character(bilanhh[i])
	if(is.na(testA)==T){testA="NA"}
	if(is.na(testB)==T){testB="NA"}
	if(testA != testB & testB!="NA"){
		Arrows(log10(HH$netdivAB_avg)[i], pmig_HH[i], log10(hh$netdivAB_avg)[i], pmig_hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(0,1,0, 0.75), lwd=1.75)
	}
	if(testA != testB & testB=="NA"){
		Arrows(log10(HH$netdivAB_avg)[i], pmig_HH[i], log10(hh$netdivAB_avg)[i], pmig_hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(1,0,0, 0.75), lwd=1.75)
	}
}


dev.print(pdf, "figureS5.pdf", bg="white")
dev.off()

