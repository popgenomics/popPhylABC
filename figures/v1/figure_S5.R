library(shape)
seuil1 = 0.6419199
seuil2 = 0.1304469

x = read.csv("../tables/table_S2.csv", h=T)

bilanHH = x$Pongoing_migration_Mhetero_Nhetero
bilanHh = x$Pongoing_migration_Mhetero_Nhomo
bilanhH = x$Pongoing_migration_Mhomo_Nhetero
bilanhh = x$Pongoing_migration_Mhomo_Nhomo

pmig_HH = bilanHH 
pmig_Hh = bilanHh 
pmig_hH = bilanhH 
pmig_hh = bilanhh 

tmp = rep(NA, length(bilanHH))
tmp[which(bilanHH>=seuil1)] = "Migration"
tmp[which(bilanHH<=seuil2)] = "NoMigration"
bilanHH = tmp

tmp = rep(NA, length(bilanHh))
tmp[which(bilanHh>=seuil1)] = "Migration"
tmp[which(bilanHh<=seuil2)] = "NoMigration"
bilanHh = tmp

tmp = rep(NA, length(bilanhH))
tmp[which(bilanhH>=seuil1)] = "Migration"
tmp[which(bilanhH<=seuil2)] = "NoMigration"
bilanhH = tmp

tmp = rep(NA, length(bilanhh))
tmp[which(bilanhh>=seuil1)] = "Migration"
tmp[which(bilanhh<=seuil2)] = "NoMigration"
bilanhh = tmp

# plot MHetero_NHetero
dev.new(width=10, height=8)
par(las=2, mar=c(5,4.5,2,2), mfrow=c(2,2))

# plot HH
plot(log10(x$netdivAB_avg), pmig_HH, ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
mtext(side = 3, text = "hetero M + Ne", cex = 1.5)

div_NA = log10(sort(x[is.na(bilanHH),]$netdivAB_avg))
div_Mig = log10(sort(x[which(bilanHH=="Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(x[which(bilanHH=="NoMigration"),]$netdivAB_avg))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
#xA2 = div_NA[1]
#xB2 = div_NA[length(div_NA)]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
#rect(xA2, -1, xB2, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(seuil1, seuil2), lty=2, lwd=1.25)

#points(log10(HH$netdivAB_avg), pmig_HH, col = rgb(0,0,0, 0.5))
par(las=1)
#mtext("MigHetero_NHetero", side = 3, line = 0.25)

# formes des points en fonction du statu d'espces
formes = rep(0, nrow(x))
formes[which(x$species_status_NG == 0)] = 22 
formes[which(x$species_status_NG == 0.5)] = 22
formes[which(x$species_status_NG == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(x))
#heteroM = apply(x[, c(39, 40, 43, 44, 47, 48)], FUN="sum", MARGIN=1)
pattern=c("Mhetero_Nhetero", "Hetero")
selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(x))))
heteroM = apply(x[, selectedCol], FUN="sum", MARGIN=1)

couleurs[which(pmig_HH>= seuil1 & heteroM >= seuil1)] = "purple"
couleurs[which(pmig_HH>= seuil1 & heteroM <= seuil1)] = "turquoise"
couleurs[which(pmig_HH<= seuil2)] = "red"
couleurs[which(pmig_HH> seuil2 & pmig_HH<seuil1)] = grey(0.25)

points(log10(x$netdivAB_avg), pmig_HH, pch=formes, col="black", bg=couleurs, cex=1.8)


# HETERO M + HOMO N
# plot Hh
plot(log10(x$netdivAB_avg), pmig_Hh, ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)

div_NA = log10(sort(x[is.na(bilanHh),]$netdivAB_avg))
div_Mig = log10(sort(x[which(bilanHh=="Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(x[which(bilanHh=="NoMigration"),]$netdivAB_avg))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
#xA2 = div_NA[1]
#xB2 = div_NA[length(div_NA)]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
#rect(xA2, -1, xB2, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(seuil1, seuil2), lty=2, lwd=1.25)

par(las=1)

# formes des points en fonction du statu d'espces
formes = rep(0, nrow(x))
formes[which(x$species_status_NG == 0)] = 22 
formes[which(x$species_status_NG == 0.5)] = 22
formes[which(x$species_status_NG == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(x))
pattern=c("Mhetero_Nhomo", "Hetero")
selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(x))))
heteroM = apply(x[, selectedCol], FUN="sum", MARGIN=1)

#homoM_homoN: species with homogeneity for both introgression and Ne
couleurs[which(pmig_Hh>= seuil1 & heteroM >= seuil1)] = "purple"
couleurs[which(pmig_Hh>= seuil1 & heteroM <= seuil1)] = "turquoise"
couleurs[which(pmig_Hh<= seuil2)] = "red"
couleurs[which(pmig_Hh> seuil2 & pmig_Hh<seuil1)] = grey(0.25)

points(log10(x$netdivAB_avg), pmig_Hh, pch=formes, col="black", bg=couleurs, cex=1.8)
mtext(side = 3, text = "hetero M + homo Ne", cex = 1.5)
for(i in 1:length(bilanHH)){
	testA = as.character(bilanHH[i])
	testB = as.character(bilanHh[i])
	if(is.na(testA)==T){testA="NA"}
	if(is.na(testB)==T){testB="NA"}
	if(testA != testB & testB!="NA"){
		Arrows(log10(x$netdivAB_avg)[i], pmig_HH[i], log10(x$netdivAB_avg)[i], pmig_Hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(0,1,0, 0.75), lwd = 2.75)
	}
	if(testA != testB & testB=="NA"){
		Arrows(log10(x$netdivAB_avg)[i], pmig_HH[i], log10(x$netdivAB_avg)[i], pmig_Hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(1,0,0, 0.75), lwd = 2.75)
	}
}

# HOMO M + HETERO N 
# plot hH 
plot(log10(x$netdivAB_avg), pmig_hH, ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)

div_NA = log10(sort(x[is.na(bilanhH),]$netdivAB_avg))
div_Mig = log10(sort(x[which(bilanhH == "Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(x[which(bilanhH == "NoMigration"),]$netdivAB_avg))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
#xA2 = div_NA[1]
#xB2 = div_NA[length(div_NA)]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
#rect(xA2, -1, xB2, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(seuil1, seuil2), lty=2, lwd=1.25)

par(las=1)

# formes des points en fonction du statu d'espces
formes = rep(0, nrow(x))
formes[which(x$species_status_NG == 0)] = 22 
formes[which(x$species_status_NG == 0.5)] = 22
formes[which(x$species_status_NG == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(x))
heteroM = rep(0, nrow(x)) 

couleurs[which(pmig_hH>= seuil1 & heteroM >= seuil1)] = "purple"
couleurs[which(pmig_hH>= seuil1 & heteroM <= seuil1)] = "turquoise"
couleurs[which(pmig_hH<= seuil2)] = "red"
couleurs[which(pmig_hH> seuil2 & pmig_hH<seuil1)] = grey(0.25)

points(log10(x$netdivAB_avg), pmig_hH, pch=formes, col="black", bg=couleurs, cex=1.8)
mtext(side = 3, text = "homo M + hetero Ne", cex = 1.5)
for(i in 1:length(bilanHH)){
	testA = as.character(bilanHH[i])
	testB = as.character(bilanhH[i])
	if(is.na(testA)==T){testA="NA"}
	if(is.na(testB)==T){testB="NA"}
	if(testA != testB & testB!="NA"){
		Arrows(log10(x$netdivAB_avg)[i], pmig_HH[i], log10(x$netdivAB_avg)[i], pmig_hH[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(0,1,0, 0.75), lwd = 2.75)
	}
	if(testA != testB & testB=="NA"){
		Arrows(log10(x$netdivAB_avg)[i], pmig_HH[i], log10(x$netdivAB_avg)[i], pmig_hH[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(1,0,0, 0.75), lwd = 2.75)
	}
}

# HOMO M + HOMO N 
# plot hh 
plot(log10(x$netdivAB_avg), pmig_hh, ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)

div_NA = log10(sort(x[is.na(bilanhh),]$netdivAB_avg))
div_Mig = log10(sort(x[which(bilanhh == "Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(x[which(bilanhh == "NoMigration"),]$netdivAB_avg))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
#xA2 = div_NA[1]
#xB2 = div_NA[length(div_NA)]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
#rect(xA2, -1, xB2, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(seuil1, seuil2), lty=2, lwd=1.25)

par(las=1)

# formes des points en fonction du statu d'espces
formes = rep(0, nrow(x))
formes[which(x$species_status_NG == 0)] = 22 
formes[which(x$species_status_NG == 0.5)] = 22
formes[which(x$species_status_NG == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(x))
heteroM = rep(0, nrow(x)) 

#homoM_homoN: species with homogeneity for both introgression and Ne
couleurs[which(pmig_hh>= seuil1 & heteroM >= seuil1)] = "purple"
couleurs[which(pmig_hh>= seuil1 & heteroM <= seuil1)] = "turquoise"
couleurs[which(pmig_hh<= seuil2)] = "red"
couleurs[which(pmig_hh> seuil2 & pmig_hh<seuil1)] = grey(0.25)

points(log10(x$netdivAB_avg), pmig_hh, pch=formes, col="black", bg=couleurs, cex=1.8)
mtext(side = 3, text = "homo M + Ne", cex = 1.5)
for(i in 1:length(bilanHH)){
	testA = as.character(bilanHH[i])
	testB = as.character(bilanhh[i])
	if(is.na(testA)==T){testA="NA"}
	if(is.na(testB)==T){testB="NA"}
	if(testA != testB & testB!="NA"){
		Arrows(log10(x$netdivAB_avg)[i], pmig_HH[i], log10(x$netdivAB_avg)[i], pmig_hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(0,1,0, 0.75), lwd = 2.75)
	}
	if(testA != testB & testB=="NA"){
		Arrows(log10(x$netdivAB_avg)[i], pmig_HH[i], log10(x$netdivAB_avg)[i], pmig_hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(1,0,0, 0.75), lwd = 2.75)
	}
}

dev.print(pdf, "figureS6new.pdf", bg="white")
dev.off()

