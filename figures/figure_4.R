library(shape)

x = read.csv("../tables/table_S2.csv", h=T)

tmp = rep(NA, nrow(x))

tmp[which(x$Pongoing_migration_Mhetero_Nhetero>=0.8)] = "Migration"
tmp[which(x$Pongoing_migration_Mhetero_Nhetero<0.2)] = "NoMigration"

bilanHH = tmp

pmig_HH = x$Pongoing_migration_Mhetero_Nhetero 

# plot MHetero_NHetero
dev.new(width=9, height=5)
par(las=2, mar=c(5,4,4,12))
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

points(log10(x$netdivAB_avg), pmig_HH, pch=formes, col="black", bg=couleurs, cex=1.8)

dev.print(pdf, "figure4.pdf", bg="white")
dev.off()

