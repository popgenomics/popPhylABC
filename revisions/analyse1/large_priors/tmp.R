library(shape)
setwd("/home/croux/Documents/popPhylABC/figures/popPhylABC/figures")
x = read.csv("../tables/table_S2.csv", h=T)

tmp = rep(NA, nrow(x))

tmp[which(x$Pongoing_migration_Mhetero_Nhetero>=0.8)] = "Migration"
tmp[which(x$Pongoing_migration_Mhetero_Nhetero<0.2)] = "NoMigration"

bilanHH = tmp

pmig_HH = x$Pongoing_migration_Mhetero_Nhetero


# plot MHetero_NHetero
dev.new(width=9, height=15)
par(las=2, mar=c(5,4,1,1), mfrow=c(3, 1))
plot(log10(x$netdivAB_avg), pmig_HH, ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression(P["ongoing migration"]), line=2.25, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)

div_NA = log10(sort(x[is.na(bilanHH),]$netdivAB_avg))
div_Mig = log10(sort(x[which(bilanHH=="Migration"),]$netdivAB_avg))
div_NoMig = log10(sort(x[which(bilanHH=="NoMigration"),]$netdivAB_avg))

xA1 = div_NoMig[1]
xB1 = div_Mig[length(div_Mig)]
xA2 = div_NA[1]
xB2 = div_NA[length(div_NA)]


rect(xA1, -1, xB1, 2, col = rgb(0, 0, 0, 0.2), border=NA)
rect(xA2, -1, xB2, 2, col = rgb(0, 0, 0, 0.2), border=NA)
abline(h=c(0.2, 0.8), lty=2, lwd=1.25)

par(las=1)
# formes des points en fonction du statu d'espces
formes = rep(0, nrow(x))
formes[which(x$species_status_NG == 0)] = 22
formes[which(x$species_status_NG == 0.5)] = 22
formes[which(x$species_status_NG == 1)] = 21

# couleurs en fonction de NHetero, MHetero, NHomo, MHomo.
couleurs = rep(0, nrow(x))

heteroM = apply(x[,c(39, 40, 43, 44, 47, 48)], FUN="sum", MARGIN=1)

#homoM_homoN: species with homogeneity for both introgression and Ne
couleurs[which(pmig_HH>= 0.8 & heteroM >= 0.8)] = "purple"
couleurs[which(pmig_HH>= 0.8 & heteroM <= 0.8)] = "turquoise"
couleurs[which(pmig_HH<= 0.2)] = "red"
couleurs[which(pmig_HH> 0.2 & pmig_HH<0.8)] = grey(0.25)

couleurs_contour = couleurs
couleurs_fond = couleurs
pNe = apply(x[, grep("v2", colnames(x))], MARGIN=1, FUN="sum")
couleurs_fond[which(pNe < 0.8)] = rgb(0, 0, 0, 0)

points(log10(x$netdivAB_avg), pmig_HH, pch=formes, col=couleurs_contour, bg=couleurs_fond, cex=1.8, lwd=2.5)

rm(list=ls())














library(vioplot)
setwd("/home/croux/Documents/popPhylABC/figures/popPhylABC/revisions/large_priors")

SI = read.table("SI", h=T)

AM = read.table("AM", h=T)

IM = read.table("IM", h=T)

SC = read.table("SC", h=T)

PAN = read.table("PAN", h=T)

#par(mfrow=c(2,2))
#plot(density(c(0, 1, SI$IM + SI$SC + SI$PAN)), main="SI", xlab=expression(P["ongoing migration"]))
#plot(density(c(0, 1, AM$IM + AM$SC + AM$PAN)), main="AM", xlab=expression(P["ongoing migration"]))
#plot(density(c(0, 1, IM$IM + IM$SC + IM$PAN)), main="IM", xlab=expression(P["ongoing migration"]))
#plot(density(c(0, 1, SC$IM + SC$SC + SC$PAN)), main="SC", xlab=expression(P["ongoing migration"]))

mig_given_mig = c(IM$IM + IM$SC + IM$PAN, SC$IM + SC$SC + SC$PAN)
mig_given_nomig = c(SI$IM + SI$SC + SI$PAN, AM$IM + AM$SC + AM$PAN)
nomig_given_mig = c(IM$SI + IM$AM , SC$SI + SC$AM, PAN$SI + PAN$AM)
nomig_given_nomig = c(SI$SI + SI$AM , AM$SI + AM$AM)


# robustness when the target is SI for 2,000 PODs
res = NULL
for(i in 1:nrow(SI)){
	pmig = SI$IM[i] + SI$SC[i] + SI$PAN[i]
	pnomig = SI$SI[i] + SI$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = F
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = T
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

SI = cbind(SI, res)


# robustness when the target is AM for 2,000 PODs
res = NULL
for(i in 1:nrow(AM)){
	pmig = AM$IM[i] + AM$SC[i] + AM$PAN[i]
	pnomig = AM$SI[i] + AM$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = F
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = T
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

AM = cbind(AM, res)


# robustness when the target is IM for 2,000 PODs
res = NULL
for(i in 1:nrow(IM)){
	pmig = IM$IM[i] + IM$SC[i] + IM$PAN[i]
	pnomig = IM$SI[i] + IM$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = T
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = F
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

IM = cbind(IM, res)


# robustness when the target is SC for 2,000 PODs
res = NULL
for(i in 1:nrow(SC)){
	pmig = SC$IM[i] + SC$SC[i] + SC$PAN[i]
	pnomig = SC$SI[i] + SC$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = T
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = F
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

SC = cbind(SC, res)

# robustness when the target is PAN for 2,000 PODs
res = NULL
for(i in 1:nrow(PAN)){
	pmig = PAN$IM[i] + PAN$SC[i] + PAN$PAN[i]
	pnomig = PAN$SI[i] + PAN$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = T
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = F
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

PAN = cbind(PAN, res)


# robustness as a function of divergence
tmpSI = cbind(SI$netdivAB_avg, SI$IM, SI$SC, SI$PAN, SI$robustness, SI$correctly_infered)
tmpAM = cbind(AM$netdivAB_avg, AM$IM, AM$SC, AM$PAN, AM$robustness, AM$correctly_infered)
tmpIM = cbind(IM$netdivAB_avg, IM$IM, IM$SC, IM$PAN, IM$robustness, IM$correctly_infered)
tmpSC = cbind(SC$netdivAB_avg, SC$IM, SC$SC, SC$PAN, SC$robustness, SC$correctly_infered)
tmpPAN = cbind(PAN$netdivAB_avg, PAN$IM, PAN$SC, PAN$PAN, PAN$robustness, PAN$correctly_infered)

res = rbind(tmpSI, tmpAM, tmpIM, tmpSC, tmpPAN)
colnames(res) = c("netdivAB_avg", "IM", "SC", "PAN", "robustness", "correctly_infered")
res = as.data.frame(res)

good = which(res$correctly_infered==1)
pasGood = which(res$correctly_infered==0)

couleurs = vector(length=nrow(res))
couleurs[good] = "green"
couleurs[pasGood] = "red"

shapes = vector(length=nrow(res))
shapes[good] = 1
shapes[pasGood] = 4

#dev.new(width=9, height=12)
#par(las=2, mar=c(5,4,1,1), mfrow=c(2, 1))
plot(log10(res$netdivAB_avg), res$robustness, ylab ="" , xlab = "net synonymous divergence", col = "white", cex.lab = 1.25, xaxt = "n", yaxt = "n", xlim=c(log10(1e-5), log10(1)), cex.axis=1.25)
par(las=3)
mtext(side=2, text=expression("robustness"), line=2.5, cex=1.5)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), cex.axis=1.25)

points(log10(res$netdivAB_avg)[good], res$robustness[good], col = couleurs[good], pch = shapes[good])
points(log10(res$netdivAB_avg)[pasGood], res$robustness[pasGood], col = couleurs[pasGood], pch = shapes[pasGood])

nomigration = na.omit(log10(c(tmpAM[,1], tmpSI[,1])))
migration = na.omit(log10(c(tmpSC[,1], tmpIM[,1])))
par(las=3)
#boxplot(nomigration, migration, names = c("no\nmigration", "migration"), horizontal=T, ylim=c(log10(1e-5), log10(1)), xaxt = "n", xlab = "net synonymous divergence", cex.lab = 1.25, cex.axis = 1.25)
vioplot(nomigration, migration, names = c("no\nmigration", "migration"), horizontal=T, ylim=c(log10(1e-5), log10(1)), col="white", lwd=2)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)

