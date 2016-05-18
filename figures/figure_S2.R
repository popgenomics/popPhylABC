x = read.table("migHetero_NHomo.txt", h=T)

# get data from lepus and eunicella
lepus = read.table("../data_popPhylABC/Lepus/SumStats.txt",h=T)
pointBLepus = c(log10(mean(lepus$netdivAB_avg)), sd(lepus$FST_avg))
pointALepus = c(log10(mean(lepus$netdivAB_avg)), mean(lepus$FST_avg))
lepus$FST_avg[which(lepus$FST_avg<0)] = 0

eunicella1 = read.table("../data_popPhylABC/Eunicella1/SumStats.txt", h=T)
pointBEunicella1 = c(log10(mean(eunicella1$netdivAB_avg)), sd(eunicella1$FST_avg))
pointAEunicella1 = c(log10(mean(eunicella1$netdivAB_avg)), mean(eunicella1$FST_avg))
eunicella1$FST_avg[which(eunicella1$FST_avg<0)] = 0

crepidula2 = read.table("../data_popPhylABC/Crepidula2/SumStats.txt", h=T)
pointBCrepidula2 = c(log10(mean(crepidula2$netdivAB_avg)), sd(crepidula2$FST_avg))
pointACrepidula2 = c(log10(mean(crepidula2$netdivAB_avg)), mean(crepidula2$FST_avg))
crepidula2$FST_avg[which(crepidula2$FST_avg<0)] = 0

# plot
#pdf("figure2.pdf", bg="white")
dev.new(width=9, height=3.5)
par(mfrow=c(1,3), las=2)
plot(log10(x$netdivAB_avg), x$FST_avg, pch=16, col=grey(0.5), cex=1.2, ylim=c(0,1), xlim=c(log10(1e-5), log10(1)), xlab = "net synonymous divergence", ylab = expression(paste("mean", " ", F["ST"], sep=" ")), cex.lab = 1.2, xaxt="n", cex.axis=1.2)
points(pointALepus[1], pointALepus[2], pch=16, col="blue", cex=1.3)
points(pointAEunicella1[1], pointAEunicella1[2], pch=16, col="green", cex=1.3)
points(pointACrepidula2[1], pointACrepidula2[2], pch=16, col="red", cex=1.3)
legend("topleft", legend=c("Lepus", "Eunicella", "Crepidula"), lty=1, col=c("blue", "green", "red"), bty="n", cex=1.1, lwd=2)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.2)

plot(log10(x$netdivAB_avg), x$FST_std, pch=16, col=grey(0.5), cex=1.2, xlab = "net synonymous divergence", ylab = expression(paste("std", " ", F["ST"], sep=" ")), cex.lab = 1.2, xaxt="n", ylim=c(0,0.5), xlim=c(log10(1e-5), log10(1)), cex.axis=1.2)
points(pointBLepus[1], pointBLepus[2], pch=16, col="blue", cex=1.3)
points(pointBEunicella1[1], pointBEunicella1[2], pch=16, col="green", cex=1.3)
points(pointBCrepidula2[1], pointBCrepidula2[2], pch=16, col="red", cex=1.3)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.2)


#dev.off()


lepus$FST_avg[which(lepus$FST_avg<0)] = 0

eunicella1$FST_avg[which(eunicella1$FST_avg<0)] = 0

crepidula2$FST_avg[which(crepidula2$FST_avg<0)] = 0

par(las=1)
plot(density(lepus$FST_avg, bw=0.05), col="blue", xlim=c(0, 1), main="", xlab=expression(F["ST"]), cex.lab=1.3, lwd=1.3)
lines(density(eunicella1$FST_avg, bw=0.1), col="green", lwd=1.3)
lines(density(crepidula2$FST_avg, bw=0.05), col="red", lwd=1.3)

dev.print(pdf, "figureS2.pdf", bg="white")
dev.off()

