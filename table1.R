# Recognzed populations
# low divergence
length(which(div<0.02 & statu==0 & mig>seuil2  & mig<seuil1))
length(which(div<0.02 & statu==0 & mig<seuil2))
length(which(div<0.02 & statu==0 & mig>seuil1))

length(which(div<0.02 & statu==0 & mig>seuil1 & migHomo>seuil1))
length(which(div<0.02 & statu==0 & mig>seuil1 & migHetero>seuil1))
length(which(div<0.02 & statu==0 & mig>seuil1 & migHomo<seuil1 & migHetero<seuil1))

# high divergence
length(which(div>0.02 & statu==0 & mig>seuil2  & mig<seuil1))
length(which(div>0.02 & statu==0 & mig<seuil2))
length(which(div>0.02 & statu==0 & mig>seuil1))

length(which(div>0.02 & statu==0 & mig>seuil1 & migHomo>seuil1))
length(which(div>0.02 & statu==0 & mig>seuil1 & migHetero>seuil1))
length(which(div>0.02 & statu==0 & mig>seuil1 & migHomo<seuil1 & migHetero<seuil1))

# Recognzed species
# low divergence
length(which(div<0.02 & statu==1 & mig>seuil2  & mig<seuil1))
length(which(div<0.02 & statu==1 & mig<seuil2))
length(which(div<0.02 & statu==1 & mig>seuil1))

length(which(div<0.02 & statu==1 & mig>seuil1 & migHomo>seuil1))
length(which(div<0.02 & statu==1 & mig>seuil1 & migHetero>seuil1))
length(which(div<0.02 & statu==1 & mig>seuil1 & migHomo<seuil1 & migHetero<seuil1))

# high divergence
length(which(div>0.02 & statu==1 & mig>seuil2  & mig<seuil1))
length(which(div>0.02 & statu==1 & mig<seuil2))
length(which(div>0.02 & statu==1 & mig>seuil1))

length(which(div>0.02 & statu==1 & mig>seuil1 & migHomo>seuil1))
length(which(div>0.02 & statu==1 & mig>seuil1 & migHetero>seuil1))
length(which(div>0.02 & statu==1 & mig>seuil1 & migHomo<seuil1 & migHetero<seuil1))


