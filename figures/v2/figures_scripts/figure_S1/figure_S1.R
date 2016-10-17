# 4, 6, 50, 100 copies
# MHomo_NHomo
dev.new(width=5.5, height=7.5)
par(mfrow=c(4,2), mar=c(4, 3.5, 3, 1.5))
# 4copies
SI = read.table("4copies/OBS_SIv1_IMv1_AMv1_SCv1_SI")
AM = read.table("4copies/OBS_SIv1_IMv1_AMv1_SCv1_AM")
IM = read.table("4copies/OBS_SIv1_IMv1_AMv1_SCv1_IM")
SC = read.table("4copies/OBS_SIv1_IMv1_AMv1_SCv1_SC")

MAXY = max(c(density(SI[,1])$y, density(AM[,3])$y, density(IM[,2])$y, density(SC[,4])$y))*1.15

nomig_nomig = c(SI[,1]+SI[,3], AM[,1]+AM[,3])
mig_mig = c(IM[,2]+IM[,4], SC[,2]+SC[,4])
nomig_nomig = nomig_nomig[which(nomig_nomig>0.75)]
mig_mig = mig_mig[which(mig_mig>0.75)]
babar(nomig_nomig, mig_mig, legende=F, nameA="Current isolation", nameB="Current introgression", legx="topleft", xl="Posterior probability")

plot(density(SI[,1]), xlim=c(0,1), ylim=c(0,5), main="4copies", xlab="Posterior probability")
lines(density(AM[,3]), col="red")
lines(density(IM[,2]), col="blue")
lines(density(SC[,4]), col="green")


# 6copies
SI = read.table("6copies/OBS_SIv1_IMv1_AMv1_SCv1_SI")
AM = read.table("6copies/OBS_SIv1_IMv1_AMv1_SCv1_AM")
IM = read.table("6copies/OBS_SIv1_IMv1_AMv1_SCv1_IM")
SC = read.table("6copies/OBS_SIv1_IMv1_AMv1_SCv1_SC")

MAXY = max(c(density(SI[,1])$y, density(AM[,3])$y, density(IM[,2])$y, density(SC[,4])$y))*1.15

nomig_nomig = c(SI[,1]+SI[,3], AM[,1]+AM[,3])
mig_mig = c(IM[,2]+IM[,4], SC[,2]+SC[,4])
nomig_nomig = nomig_nomig[which(nomig_nomig>0.75)]
mig_mig = mig_mig[which(mig_mig>0.75)]
babar(nomig_nomig, mig_mig, legende=F, nameA="Current isolation", nameB="Current introgression", legx="topleft", xl="Posterior probability")

plot(density(SI[,1]), xlim=c(0,1), ylim=c(0,5), main="6copies", xlab="Posterior probability")
lines(density(AM[,3]), col="red")
lines(density(IM[,2]), col="blue")
lines(density(SC[,4]), col="green")


# 50copies
SI = read.table("50copies/OBS_SIv1_IMv1_AMv1_SCv1_SI")
AM = read.table("50copies/OBS_SIv1_IMv1_AMv1_SCv1_AM")
IM = read.table("50copies/OBS_SIv1_IMv1_AMv1_SCv1_IM")
SC = read.table("50copies/OBS_SIv1_IMv1_AMv1_SCv1_SC")

MAXY = max(c(density(SI[,1])$y, density(AM[,3])$y, density(IM[,2])$y, density(SC[,4])$y))*1.15

nomig_nomig = c(SI[,1]+SI[,3], AM[,1]+AM[,3])
mig_mig = c(IM[,2]+IM[,4], SC[,2]+SC[,4])
nomig_nomig = nomig_nomig[which(nomig_nomig>0.75)]
mig_mig = mig_mig[which(mig_mig>0.75)]
babar(nomig_nomig, mig_mig, legende=F, nameA="Current isolation", nameB="Current introgression", legx="topleft", xl="Posterior probability")

plot(density(SI[,1]), xlim=c(0,1), ylim=c(0,5), main="50copies", xlab="Posterior probability")
lines(density(AM[,3]), col="red")
lines(density(IM[,2]), col="blue")
lines(density(SC[,4]), col="green")


# 100copies
SI = read.table("100copies/OBS_SIv1_IMv1_AMv1_SCv1_SI")
AM = read.table("100copies/OBS_SIv1_IMv1_AMv1_SCv1_AM")
IM = read.table("100copies/OBS_SIv1_IMv1_AMv1_SCv1_IM")
SC = read.table("100copies/OBS_SIv1_IMv1_AMv1_SCv1_SC")

MAXY = max(c(density(SI[,1])$y, density(AM[,3])$y, density(IM[,2])$y, density(SC[,4])$y))*1.15

nomig_nomig = c(SI[,1]+SI[,3], AM[,1]+AM[,3])
mig_mig = c(IM[,2]+IM[,4], SC[,2]+SC[,4])
nomig_nomig = nomig_nomig[which(nomig_nomig>0.75)]
mig_mig = mig_mig[which(mig_mig>0.75)]
babar(nomig_nomig, mig_mig, legende=F, nameA="Current isolation", nameB="Current introgression", legx="topleft", xl="Posterior probability")

plot(density(SI[,1]), xlim=c(0,1), ylim=c(0,5), main="100copies", xlab="Posterior probability")
lines(density(AM[,3]), col="red")
lines(density(IM[,2]), col="blue")
lines(density(SC[,4]), col="green")

dev.print(pdf, "figureS1.pdf", bg="white")
dev.off()

