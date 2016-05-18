x = read.table("migHetero_NHomo.txt", h=T)
x[which(x$FST_avg<0),]$FST_avg = 0
# plot
#pdf("figure2.pdf", bg="white")
dev.new(width=10.5, height=5.5)
par(mfrow=c(1,2), las=2)
plot(log10(x$netdivAB_avg), x$FST_avg, pch=16, col=grey(0.5), cex=1.2, ylim=c(0,1), xlim=c(log10(1e-5), log10(1)), xlab = "net synonymous divergence", ylab = "", cex.lab = 1.5, xaxt="n", yaxt="n", cex.axis=1.35)
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.25)
axis(side = 2, cex.axis=1.35)
par(las=3)
mtext(side = 2, text = expression(paste("mean", " ", F["ST"], sep=" ")), cex=1.5, line=2.5)

par(las=1)
plot(log10(x$netdivAB_avg), x$FST_std, pch=16, col=grey(0.5), cex=1.2, xlab = "net synonymous divergence", ylab = "", cex.lab = 1.5, xaxt="n", yaxt="n", ylim=c(0,0.5), xlim=c(log10(1e-5), log10(1)), cex.axis=1.25) 
par(las=1)
axis(side = 1, at = c(0, -1, -2, -3, -4, -5), labels = c(1, 0.1, 0.01, expression(paste("10"^"-3")), expression(paste("10"^"-4")) , expression(paste("10"^"-5"))), cex.axis=1.35)
axis(side = 2, cex.axis=1.35)
par(las=3)
mtext(side = 2, text = expression(paste("std", " ", F["ST"], sep=" ")), cex=1.5, line=2.5)


#dev.off()

dev.print(pdf, "figure3.pdf", bg="white")
dev.off()

