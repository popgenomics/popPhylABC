library(shape)
seuil1 = 0.6419199
seuil2 = 0.1304469

pouet = function(x){
	return(paste(x[1], x[2], sep="_"))
}

x=read.table("../../../../supp_online/alpha.estimates",h=T)
y=read.table("../../../../supp_online/estimation_parameters.tmp",h=T)
z=read.csv("../../../tables/table_S2.csv",h=T)

names = apply(y[,1:2], MARGIN=1, FUN="pouet")
y=cbind(names, y)

names = apply(z[,1:2], MARGIN=1, FUN="pouet")
z=cbind(names, z)

tot = merge(z, y, by.x=1, by.y=1)

# 
bilanHH = tot$Pongoing_migration_Mhetero_Nhetero
bilanHh = tot$Pongoing_migration_Mhetero_Nhomo
bilanhH = tot$Pongoing_migration_Mhomo_Nhetero
bilanhh = tot$Pongoing_migration_Mhomo_Nhomo

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

# heteroM heteroN
pattern=c("Mhetero_Nhetero", "Hetero")
selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(tot))))
heteroM = apply(tot[, selectedCol], FUN="sum", MARGIN=1)

dev.new(width=10, height=8)
par(las=2, mar=c(5,4.5,2,2), mfrow=c(2,2))

# plot HH
plot((tot$medianT), pmig_HH, ylab ="" , xlab = "estimated Tsplit", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(0, 30), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
mtext(side = 3, text = "hetero M + Ne", cex = 1.5)

div_NA = (sort(tot[is.na(bilanHH),]$medianT))
div_Mig = (sort(tot[which(bilanHH=="Migration"),]$medianT))
div_NoMig = (sort(tot[which(bilanHH=="NoMigration"),]$medianT))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
xA2 = div_NA[1]
xB2 = div_NA[length(div_NA)]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(seuil2, seuil1), lty=2, lwd=1.25)

par(las=1)
# formes des points en fonction du statu d'espces
formes = rep(0, nrow(tot))
formes[which(tot$species_status_NG == 0)] = 22 
formes[which(tot$species_status_NG == 0.5)] = 22
formes[which(tot$species_status_NG == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(tot))
pattern=c("Mhetero_Nhetero", "Hetero")
selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(tot))))
heteroM = apply(tot[, selectedCol], FUN="sum", MARGIN=1)


couleurs[which(pmig_HH>= seuil1 & heteroM >= seuil1)] = "purple"
couleurs[which(pmig_HH>= seuil1 & heteroM <= seuil1)] = "turquoise"
couleurs[which(pmig_HH<= seuil2)] = "red"
couleurs[which(pmig_HH> seuil2 & pmig_HH<seuil1)] = grey(0.25)

points(tot$medianT, pmig_HH, pch=formes, col="black", bg=couleurs, cex=1.8)



# plot Hh
plot((tot$medianT), pmig_Hh, ylab ="" , xlab = "estimated Tsplit", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(0, 30), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
mtext(side = 3, text = "hetero M", cex = 1.5)

div_NA = (sort(tot[is.na(bilanHh),]$medianT))
div_Mig = (sort(tot[which(bilanHh=="Migration"),]$medianT))
div_NoMig = (sort(tot[which(bilanHh=="NoMigration"),]$medianT))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
xA2 = div_NA[1]
xB2 = div_NA[length(div_NA)]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(seuil2, seuil1), lty=2, lwd=1.25)

par(las=1)
# formes des points en fonction du statu d'espces
formes = rep(0, nrow(tot))
formes[which(tot$species_status_NG == 0)] = 22 
formes[which(tot$species_status_NG == 0.5)] = 22
formes[which(tot$species_status_NG == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(tot))
pattern=c("Mhetero_Nhomo", "Hetero")
selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(tot))))
heteroM = apply(tot[, selectedCol], FUN="sum", MARGIN=1)


couleurs[which(pmig_Hh>= seuil1 & heteroM >= seuil1)] = "purple"
couleurs[which(pmig_Hh>= seuil1 & heteroM <= seuil1)] = "turquoise"
couleurs[which(pmig_Hh<= seuil2)] = "red"
couleurs[which(pmig_Hh> seuil2 & pmig_Hh<seuil1)] = grey(0.25)

points(tot$medianT, pmig_Hh, pch=formes, col="black", bg=couleurs, cex=1.8)

#
for(i in 1:length(bilanHH)){
	testA = as.character(bilanHH[i])
	testB = as.character(bilanHh[i])
	if(is.na(testA)==T){testA="NA"}
	if(is.na(testB)==T){testB="NA"}
	if(testA != testB & testB!="NA"){
		Arrows(tot$medianT[i], pmig_HH[i], tot$medianT[i], pmig_Hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(0,1,0, 0.75), lwd = 2.75)
	}
	if(testA != testB & testB=="NA"){
		Arrows(tot$medianT[i], pmig_HH[i], tot$medianT[i], pmig_Hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(1,0,0, 0.75), lwd = 2.75)
	}
}


# plot hH 
plot((tot$medianT), pmig_hH, ylab ="" , xlab = "estimated Tsplit", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(0, 30), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
mtext(side = 3, text = "hetero Ne", cex = 1.5)

div_NA = (sort(tot[is.na(bilanhH),]$medianT))
div_Mig = (sort(tot[which(bilanhH=="Migration"),]$medianT))
div_NoMig = (sort(tot[which(bilanhH=="NoMigration"),]$medianT))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
xA2 = div_NA[1]
xB2 = div_NA[length(div_NA)]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(seuil2, seuil1), lty=2, lwd=1.25)

par(las=1)
# formes des points en fonction du statu d'espces
formes = rep(0, nrow(tot))
formes[which(tot$species_status_NG == 0)] = 22 
formes[which(tot$species_status_NG == 0.5)] = 22
formes[which(tot$species_status_NG == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(tot))
pattern=c("Mhomo_Nhetero", "Hetero")
selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(tot))))
heteroM = apply(tot[, selectedCol], FUN="sum", MARGIN=1)


couleurs[which(pmig_hH>= seuil1 & heteroM >= seuil1)] = "purple"
couleurs[which(pmig_hH>= seuil1 & heteroM <= seuil1)] = "turquoise"
couleurs[which(pmig_hH<= seuil2)] = "red"
couleurs[which(pmig_hH> seuil2 & pmig_hH<seuil1)] = grey(0.25)

points(tot$medianT, pmig_hH, pch=formes, col="black", bg=couleurs, cex=1.8)

#
for(i in 1:length(bilanHH)){
	testA = as.character(bilanHH[i])
	testB = as.character(bilanhH[i])
	if(is.na(testA)==T){testA="NA"}
	if(is.na(testB)==T){testB="NA"}
	if(testA != testB & testB!="NA"){
		Arrows(tot$medianT[i], pmig_HH[i], tot$medianT[i], pmig_hH[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(0,1,0, 0.75), lwd = 2.75)
	}
	if(testA != testB & testB=="NA"){
		Arrows(tot$medianT[i], pmig_HH[i], tot$medianT[i], pmig_hH[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(1,0,0, 0.75), lwd = 2.75)
	}
}


# plot hh 
plot((tot$medianT), pmig_hh, ylab ="" , xlab = "estimated Tsplit", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(0, 30), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)
mtext(side = 3, text = "homo M + Ne", cex = 1.5)

div_NA = (sort(tot[is.na(bilanhh),]$medianT))
div_Mig = (sort(tot[which(bilanhh=="Migration"),]$medianT))
div_NoMig = (sort(tot[which(bilanhh=="NoMigration"),]$medianT))
xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
xA2 = div_NA[1]
xB2 = div_NA[length(div_NA)]

rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(seuil2, seuil1), lty=2, lwd=1.25)

par(las=1)
# formes des points en fonction du statu d'espces
formes = rep(0, nrow(tot))
formes[which(tot$species_status_NG == 0)] = 22 
formes[which(tot$species_status_NG == 0.5)] = 22
formes[which(tot$species_status_NG == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(tot))
pattern=c("Mhomo_Nhomo", "Hetero")
selectedCol = which(Reduce('&', lapply(pattern, grepl, colnames(tot))))
heteroM = apply(tot[, selectedCol], FUN="sum", MARGIN=1)


couleurs[which(pmig_hh>= seuil1 & heteroM >= seuil1)] = "purple"
couleurs[which(pmig_hh>= seuil1 & heteroM <= seuil1)] = "turquoise"
couleurs[which(pmig_hh<= seuil2)] = "red"
couleurs[which(pmig_hh> seuil2 & pmig_hh<seuil1)] = grey(0.25)

points(tot$medianT, pmig_hh, pch=formes, col="black", bg=couleurs, cex=1.8)

#
for(i in 1:length(bilanHH)){
	testA = as.character(bilanHH[i])
	testB = as.character(bilanhh[i])
	if(is.na(testA)==T){testA="NA"}
	if(is.na(testB)==T){testB="NA"}
	if(testA != testB & testB!="NA"){
		Arrows(tot$medianT[i], pmig_HH[i], tot$medianT[i], pmig_hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(0,1,0, 0.75), lwd = 2.75)
	}
	if(testA != testB & testB=="NA"){
		Arrows(tot$medianT[i], pmig_HH[i], tot$medianT[i], pmig_hh[i], arr.length=0.1, arr.type="triangle", arr.width=0.1, lty=4, col=rgb(1,0,0, 0.75), lwd = 2.75)
	}
}

dev.print(pdf, "/home/croux/Documents/popPhylABC/soumission_v2/figure_Tsplit_migration.pdf", bg="white")
dev.off()

