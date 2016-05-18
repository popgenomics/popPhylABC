# 4, 6, 50, 100 copies
# MHomo_NHomo
par(mfrow=c(4,2))
# 4copies
SI = read.table("/home/croux/Documents/test_robustesse_popPhyl/4copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_SI")
AM = read.table("/home/croux/Documents/test_robustesse_popPhyl/4copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_AM")
IM = read.table("/home/croux/Documents/test_robustesse_popPhyl/4copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_IM")
SC = read.table("/home/croux/Documents/test_robustesse_popPhyl/4copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_SC")

MAXY = max(c(density(SI[,1])$y, density(AM[,3])$y, density(IM[,2])$y, density(SC[,4])$y))*1.15

plot(density(SI[,1]), xlim=c(0,1), ylim=c(0,5), main="4copies", xlab="Posterior probability")
lines(density(AM[,3]), col="red")
lines(density(IM[,2]), col="blue")
lines(density(SC[,4]), col="green")

nomig_nomig = c(SI[,1]+SI[,3], AM[,1]+AM[,3])
mig_mig = c(IM[,2]+IM[,4], SC[,2]+SC[,4])
nomig_nomig = nomig_nomig[which(nomig_nomig>0.75)]
mig_mig = mig_mig[which(mig_mig>0.75)]
babar(nomig_nomig, mig_mig, legende=F, nameA="Current isolation", nameB="Current introgression", legx="topleft", xl="Posterior probability")

# 6copies
SI = read.table("/home/croux/Documents/test_robustesse_popPhyl/6copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_SI")
AM = read.table("/home/croux/Documents/test_robustesse_popPhyl/6copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_AM")
IM = read.table("/home/croux/Documents/test_robustesse_popPhyl/6copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_IM")
SC = read.table("/home/croux/Documents/test_robustesse_popPhyl/6copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_SC")

MAXY = max(c(density(SI[,1])$y, density(AM[,3])$y, density(IM[,2])$y, density(SC[,4])$y))*1.15

plot(density(SI[,1]), xlim=c(0,1), ylim=c(0,5), main="6copies", xlab="Posterior probability")
lines(density(AM[,3]), col="red")
lines(density(IM[,2]), col="blue")
lines(density(SC[,4]), col="green")

nomig_nomig = c(SI[,1]+SI[,3], AM[,1]+AM[,3])
mig_mig = c(IM[,2]+IM[,4], SC[,2]+SC[,4])
nomig_nomig = nomig_nomig[which(nomig_nomig>0.75)]
mig_mig = mig_mig[which(mig_mig>0.75)]
babar(nomig_nomig, mig_mig, legende=F, nameA="Current isolation", nameB="Current introgression", legx="topleft", xl="Posterior probability")

# 50copies
SI = read.table("/home/croux/Documents/test_robustesse_popPhyl/50copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_SI")
AM = read.table("/home/croux/Documents/test_robustesse_popPhyl/50copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_AM")
IM = read.table("/home/croux/Documents/test_robustesse_popPhyl/50copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_IM")
SC = read.table("/home/croux/Documents/test_robustesse_popPhyl/50copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_SC")

MAXY = max(c(density(SI[,1])$y, density(AM[,3])$y, density(IM[,2])$y, density(SC[,4])$y))*1.15

plot(density(SI[,1]), xlim=c(0,1), ylim=c(0,5), main="50copies", xlab="Posterior probability")
lines(density(AM[,3]), col="red")
lines(density(IM[,2]), col="blue")
lines(density(SC[,4]), col="green")

nomig_nomig = c(SI[,1]+SI[,3], AM[,1]+AM[,3])
mig_mig = c(IM[,2]+IM[,4], SC[,2]+SC[,4])
nomig_nomig = nomig_nomig[which(nomig_nomig>0.75)]
mig_mig = mig_mig[which(mig_mig>0.75)]
babar(nomig_nomig, mig_mig, legende=F, nameA="Current isolation", nameB="Current introgression", legx="topleft", xl="Posterior probability")

# 100copies
SI = read.table("/home/croux/Documents/test_robustesse_popPhyl/100copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_SI")
AM = read.table("/home/croux/Documents/test_robustesse_popPhyl/100copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_AM")
IM = read.table("/home/croux/Documents/test_robustesse_popPhyl/100copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_IM")
SC = read.table("/home/croux/Documents/test_robustesse_popPhyl/100copies/MHomo_NHomo/OBS_SIv1_IMv1_AMv1_SCv1_SC")

MAXY = max(c(density(SI[,1])$y, density(AM[,3])$y, density(IM[,2])$y, density(SC[,4])$y))*1.15

plot(density(SI[,1]), xlim=c(0,1), ylim=c(0,5), main="100copies", xlab="Posterior probability")
lines(density(AM[,3]), col="red")
lines(density(IM[,2]), col="blue")
lines(density(SC[,4]), col="green")

nomig_nomig = c(SI[,1]+SI[,3], AM[,1]+AM[,3])
mig_mig = c(IM[,2]+IM[,4], SC[,2]+SC[,4])
nomig_nomig = nomig_nomig[which(nomig_nomig>0.75)]
mig_mig = mig_mig[which(mig_mig>0.75)]
babar(nomig_nomig, mig_mig, legende=F, nameA="Current isolation", nameB="Current introgression", legx="topleft", xl="Posterior probability")

dev.print(pdf, "/home/croux/Documents/tmp/scripts_figures/fig_MHomo_NHomo.pdf", bg="white")
dev.off()

