#! /usr/bin/env Rscript
a = read.table("migHomo_NHomo.txt", h=T)
b = read.table("migHetero_NHomo.txt", h=T)
c = read.table("migHomo_NHetero.txt", h=T)
d = read.table("migHetero_NHetero.txt", h=T)

#A = paste(a$spA, a$spB, sep="_")
#B = paste(b$spA, b$spB, sep="_")
#C = paste(c$spA, c$spB, sep="_")
#D = paste(d$spA, d$spB, sep="_")

SIa = grep("SI",colnames(a))
SIb = grep("SI",colnames(b))
SIc = grep("SI",colnames(c))
SId = grep("SI",colnames(d))

AMa = grep("AM",colnames(a))
AMb = grep("AM",colnames(b))
AMc = grep("AM",colnames(c))
AMd = grep("AM",colnames(d))

IMa = grep("IM",colnames(a))
IMb = grep("IM",colnames(b))
IMc = grep("IM",colnames(c))
IMd = grep("IM",colnames(d))

SCa = grep("SC",colnames(a))
SCb = grep("SC",colnames(b))
SCc = grep("SC",colnames(c))
SCd = grep("SC",colnames(d))

PANa = grep("Panmixia",colnames(a))
PANb = grep("Panmixia",colnames(b))
PANc = grep("Panmixia",colnames(c))
PANd = grep("Panmixia",colnames(d))

isolationA = a[,SIa]
isolationB = b[,SIb]
isolationC = apply(c[,SIc], MARGIN=1, FUN="sum")
isolationD = apply(d[,SId], MARGIN=1, FUN="sum")

ancientA = a[,AMa]
ancientB = apply(b[,AMb], MARGIN=1, FUN="sum")
ancientC = apply(c[,AMc], MARGIN=1, FUN="sum")
ancientD = apply(d[,AMd], MARGIN=1, FUN="sum")

islandA = a[,IMa]
islandB = apply(b[,IMb], MARGIN=1, FUN="sum")
islandC = apply(c[,IMc], MARGIN=1, FUN="sum")
islandD = apply(d[,IMd], MARGIN=1, FUN="sum")

secondaryA = a[,SCa]
secondaryB = apply(b[,SCb], MARGIN=1, FUN="sum")
secondaryC = apply(c[,SCc], MARGIN=1, FUN="sum")
secondaryD = apply(d[,SCd], MARGIN=1, FUN="sum")

panmixieA = a[, PANa]
panmixieB = b[, PANb]
panmixieC = apply(c[, PANc], MARGIN=1, FUN="sum")
panmixieD = apply(d[, PANd], MARGIN=1, FUN="sum")

migA = islandA + secondaryA + panmixieA
migB = islandB + secondaryB + panmixieB
migC = islandC + secondaryC + panmixieC
migD = islandD + secondaryD + panmixieD

noMigA = isolationA + ancientA
noMigB = isolationB + ancientB
noMigC = isolationC + ancientC
noMigD = isolationD + ancientD

nMigA = length(which(migA > 0.8))
nMigB = length(which(migB > 0.8))
nMigC = length(which(migC > 0.8))
nMigD = length(which(migD > 0.8))

nNoMigA = length(which(noMigA > 0.8))
nNoMigB = length(which(noMigB > 0.8))
nNoMigC = length(which(noMigC > 0.8))
nNoMigD = length(which(noMigD > 0.8))

nPanA = length(which(panmixieA > 0.8))
nPanB = length(which(panmixieB > 0.8))
nPanC = length(which(panmixieC > 0.8))
nPanD = length(which(panmixieD > 0.8))

nAmbigA = length(which(migA<0.8 & noMigA<0.8 & panmixieA<0.8))
nAmbigB = length(which(migB<0.8 & noMigB<0.8 & panmixieB<0.8))
nAmbigC = length(which(migC<0.8 & noMigC<0.8 & panmixieC<0.8))
nAmbigD = length(which(migD<0.8 & noMigD<0.8 & panmixieD<0.8))

# plot resultat 1
tmp = matrix(c(nNoMigA, nNoMigB, nNoMigC, nNoMigD, nMigA, nMigB, nMigC, nMigD, nAmbigA, nAmbigB, nAmbigC, nAmbigD), nrow=4)
dev.new(width=6, height=5.6)
barplot(tmp, beside=T, names = c("Current\nisolation", "Current\nintrogression", "NA"), col=grey(c(1, 0.75, 0.5, 0)), legend.text=c("homo M + N", "hetero M", "hetero N", "hetero M + N"), args.legend=list(bty="n"), ylab = "# of pair of species")
#tmp=dev.off()
cat(paste("#no migration (M_homo, N_homo): ", nNoMigA,
"\n#no migration (M_hetero, N_homo): ", nNoMigB,
"\n#no migration (M_homo, N_hetero): ", nNoMigC,
"\n#no migration (M_hetero, N_hetero): ", nNoMigD,
"\n#migration (M_homo, N_homo): ", nMigA,
"\n#migration (M_hetero, N_homo): ", nMigB,
"\n#migration (M_homo, N_hetero): ", nMigC,
"\n#migration (M_hetero, N_hetero): ", nMigD,
"\n#panmixia (M_homo, N_homo): ", nPanA,
"\n#panmixia (M_hetero, N_homo): ", nPanB,
"\n#panmixia (M_homo, N_hetero): ", nPanC,
"\n#panmixia (M_hetero, N_hetero): ", nPanD,
"\n#ambiguous (M_homo, N_homo): ", nAmbigA,
"\n#ambiguous (M_hetero, N_homo): ", nAmbigB,
"\n#ambiguous (M_homo, N_hetero): ", nAmbigC,
"\n#ambiguous (M_hetero, N_hetero): ", nAmbigD,
"\n\n", collapse=""))

dev.print(pdf, "figureS3.pdf", bg="white")
tmp = dev.off()

nSIa = length(which(isolationA > 0.8))
nSIb = length(which(isolationB > 0.8))
nSIc = length(which(isolationC > 0.8))
nSId = length(which(isolationD > 0.8))

nAMa = length(which(ancientA > 0.8))
nAMb = length(which(ancientB > 0.8))
nAMc = length(which(ancientC > 0.8))
nAMd = length(which(ancientD > 0.8))

nIMa = length(which(islandA > 0.8))
nIMb = length(which(islandB > 0.8))
nIMc = length(which(islandC > 0.8))
nIMd = length(which(islandD > 0.8))

nSCa = length(which(secondaryA > 0.8))
nSCb = length(which(secondaryB > 0.8))
nSCc = length(which(secondaryC > 0.8))
nSCd = length(which(secondaryD > 0.8))

nPANa = length(which(panmixieA > 0.8))
nPANb = length(which(panmixieB > 0.8))
nPANc = length(which(panmixieC > 0.8))
nPANd = length(which(panmixieD > 0.8))

nAmbigA = length(which(isolationA < 0.8 & ancientA < 0.8 & islandA < 0.8 & secondaryA < 0.8 & panmixieA < 0.8))
nAmbigB = length(which(isolationB < 0.8 & ancientB < 0.8 & islandB < 0.8 & secondaryB < 0.8 & panmixieB < 0.8))
nAmbigC = length(which(isolationC < 0.8 & ancientC < 0.8 & islandC < 0.8 & secondaryC < 0.8 & panmixieC < 0.8))
nAmbigD = length(which(isolationD < 0.8 & ancientD < 0.8 & islandD < 0.8 & secondaryD < 0.8 & panmixieD < 0.8))

# plot resultat 2
tmp = matrix(c(nSIa, nSIb, nSIc, nSId, nAMa, nAMb, nAMc, nAMd, nIMa, nIMb, nIMc, nIMd, nSCa, nSCb, nSCc, nSCd, nPANa, nPANb, nPANc, nPANd, nAmbigA, nAmbigB, nAmbigC, nAmbigD), nrow=4)
dev.new(width=8, height=5.6)
barplot(tmp, beside=T, names = c("SI", "AM", "IM", "SC", "PAN", "NA"), col=grey(c(1, 0.75, 0.5, 0)), legend.text=c("homo M + N", "hetero M", "hetero N", "hetero M + N"), args.legend=list(x="topleft", bty="n"), ylab = "# of pair of species")
dev.print(pdf, "figureS4.pdf", bg="white")
dev.off()


# unploted resultat 3
# migA = islandA + secondaryA + panmixieA
# migB = islandB + secondaryB + panmixieB
# migC = islandC + secondaryC + panmixieC
# migD = islandD + secondaryD + panmixieD
# 
# noMigA = isolationA + ancientA
# noMigB = isolationB + ancientB
# noMigC = isolationC + ancientC
# noMigD = isolationD + ancientD
# 
# boxplot(log10(a$netdivAB_avg[which(noMigA>0.8)]), log10(b$netdivAB_avg[which(noMigB>0.8)]), log10(c$netdivAB_avg[which(noMigC>0.8)]), log10(d$netdivAB_avg[which(noMigD>0.8)]),
# log10(a$netdivAB_avg[which(migA>0.8)]), log10(b$netdivAB_avg[which(migB>0.8)]), log10(c$netdivAB_avg[which(migC>0.8)]), log10(d$netdivAB_avg[which(migD>0.8)]),
# log10(a$netdivAB_avg[which(migA < 0.8 & noMigA < 0.8)]), log10(b$netdivAB_avg[which(migB < 0.8 & noMigB < 0.8)]), log10(c$netdivAB_avg[which(migC < 0.8 & noMigC < 0.8)]), log10(d$netdivAB_avg[which(migD < 0.8 & noMigC < 0.8)]),
# at = c(1,2,3,4, 6,7,8,9, 11,12,13,14), outline = F, axis=F, col=grey(c(1, 0.75, 0.5, 0, 1, 0.75, 0.5, 0, 1, 0.75, 0.5, 0)), xaxt="n", yaxt="n")
# legend("topright", c("refABC", "MigHetero", "NHetero", "MigHetero_NHetero"), bty="n", fill=grey(c(1, 0.75, 0.5, 0)))
# mtext("Current\nisolation", side=1, at=2.5, line=1.7, cex=1.2)
# mtext("Current\nintrogression", side=1, at=7.5, line=1.7, cex=1.2)
# mtext("NA", side=1, at=12.5, line=1, cex=1.2)
# par(las=0)
# mtext("synonymous divergence", side=2, line=3, cex=1.2)
# par(las=1)
# axis(side = 2, at = c(log10(0.35), -1, -2, -3, -4, log10(5e-5)), labels = c(0.35, 0.1, 0.01, 0.001, 0.0001, 0.00005))
# dev.print(pdf, "figureS2.pdf", bg="white")
# dev.off()
#tmp = dev.off()


