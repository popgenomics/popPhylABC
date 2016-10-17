SI = read.table("probas_SIv1_IMv1_AMv1_SCv1_SI")
AM = read.table("probas_SIv1_IMv1_AMv1_SCv1_AM")
IM = read.table("probas_SIv1_IMv1_AMv1_SCv1_IM")
SC = read.table("probas_SIv1_IMv1_AMv1_SCv1_SC")

MAXY = max(c(density(SI[ ,1])$y, density(AM[ ,3])$y, density(IM[ ,2])$y, density(SC[ ,4])$y))*1.15

nomig_nomig = c(SI[ ,1] + SI[ ,3], AM[ ,1] + AM[ ,3]) # P(nomig | nomig)
mig_mig = c(IM[ ,2] + IM[ ,4], SC[ ,2] + SC[ ,4]) # P(mig | mig)
mig_nomig = c(SI[ ,2] + SI[ ,4], AM[ ,2] + AM[ ,4]) # P(mig | nomig)
nomig_mig = c(IM[ ,1] + IM[ ,3], SC[ ,1] + SC[ ,3]) # P(nomig | mig)

dev.new(width = 6.5, height = 6)
babar(mig_nomig, mig_mig, legende=F, nameA="", nameB="", legx="topleft", xl="P(ongoing migration)", space=1)
dev.print(pdf, "figure2B.pdf", bg="white")
dev.off()

