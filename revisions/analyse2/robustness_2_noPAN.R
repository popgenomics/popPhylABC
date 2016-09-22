noMig = c("SI_NHetero", "SI_NHomo", "AM_MHetero_NHetero", "AM_MHetero_NHomo", "AM_MHomo_NHetero", "AM_MHomo_NHomo")
mig = c("IM_MHetero_NHetero", "IM_MHetero_NHomo", "IM_MHomo_NHetero", "IM_MHomo_NHomo", "SC_MHetero_NHetero", "SC_MHetero_NHomo", "SC_MHomo_NHetero", "SC_MHomo_NHomo")


robust = function(x, M1, M2){
	# x = proba to support M1 for a given dataset
	# M1 = best supported model
	# M2 = less supported model
	X = density(M1, from=x, to=x, width = 0.15)$y[1]
	Y = density(M2, from=x, to=x, width = 0.15)$y[1]
	return(X/(X+Y))
} 
# enter the empirical distributions
empirical_noMig_tmp = NULL
empirical_mig_tmp = NULL
for(i in noMig){
	x = read.table(paste("summary_", i, ".txt", sep=""), h=T)
	SI = grep("SI", colnames(x))
	AM = grep("AM", colnames(x))
	IM = grep("IM", colnames(x))
	SC = grep("SC", colnames(x))
#	PAN = grep("PAN", colnames(x))
	empirical_noMig_tmp = rbind(empirical_noMig_tmp, x[, c(SI, AM, IM, SC)])
#	empirical_noMig_tmp = rbind(empirical_noMig_tmp, x[, c(SI, AM, IM, SC, PAN)])

}

for(i in mig){
	x = read.table(paste("summary_", i, ".txt", sep=""), h=T)
	SI = grep("SI", colnames(x))
	AM = grep("AM", colnames(x))
	IM = grep("IM", colnames(x))
	SC = grep("SC", colnames(x))
#	PAN = grep("PAN", colnames(x))
	empirical_mig_tmp = rbind(empirical_mig_tmp, x[, c(SI, AM, IM, SC)])
#	empirical_mig_tmp = rbind(empirical_mig_tmp, x[, c(SI, AM, IM, SC, PAN)])
}

SI = grep("SI", colnames(empirical_noMig_tmp))
AM = grep("AM", colnames(empirical_noMig_tmp))
IM = grep("IM", colnames(empirical_noMig_tmp))
SC = grep("SC", colnames(empirical_noMig_tmp))
#PAN = grep("PAN", colnames(empirical_noMig_tmp))
empirical_noMig = cbind(apply(empirical_noMig_tmp[, c(SI, AM)], MARGIN=1, FUN="sum"), apply(empirical_noMig_tmp[, c(IM, SC)], MARGIN=1, FUN="sum")) # noMig mig
#empirical_noMig = cbind(apply(empirical_noMig_tmp[, c(SI, AM)], MARGIN=1, FUN="sum"), apply(empirical_noMig_tmp[, c(IM, SC, PAN)], MARGIN=1, FUN="sum")) # noMig mig

SI = grep("SI", colnames(empirical_mig_tmp))
AM = grep("AM", colnames(empirical_mig_tmp))
IM = grep("IM", colnames(empirical_mig_tmp))
SC = grep("SC", colnames(empirical_mig_tmp))
#PAN = grep("PAN", colnames(empirical_mig_tmp))
empirical_mig = cbind(apply(empirical_mig_tmp[, c(SI, AM)], MARGIN=1, FUN="sum"), apply(empirical_mig_tmp[, c(IM, SC)], MARGIN=1, FUN="sum")) # noMig mig
#empirical_mig = cbind(apply(empirical_mig_tmp[, c(SI, AM)], MARGIN=1, FUN="sum"), apply(empirical_mig_tmp[, c(IM, SC, PAN)], MARGIN=1, FUN="sum")) # noMig mig
empirical_mig = empirical_mig/apply(empirical_mig, MARGIN=1, FUN="sum")
empirical_noMig = empirical_noMig/apply(empirical_noMig, MARGIN=1, FUN="sum")

# treat the datasets
# loop over datasets simulated with ongoing migration
for(i in mig){
	print(i)
	x = read.table(paste("summary_", i, ".txt", sep=""), h=T)
	SI = grep("pSI", colnames(x))
	AM = grep("pAM", colnames(x))
	IM = grep("pIM", colnames(x))
	SC = grep("pSC", colnames(x))
#	PAN = grep("pPAN", colnames(x))

#	pMig = apply(x[, c(IM, SC, PAN)], MARGIN=1, FUN="sum")
	pMig = apply(x[, c(IM, SC)], MARGIN=1, FUN="sum")
	pNoMig = apply(x[, c(AM, SI)], MARGIN=1, FUN="sum")
	pMig2 = pMig/(pMig+pNoMig)
	pNoMig2 = pMig/(pMig+pNoMig)
	pMig = pMig2
	pNoMig = pNoMig2

	rob = NULL
	for(j in 1:length(pMig)){
		if(pMig[j] > pNoMig[j]){
			rob = c(rob, robust(x = pMig[j], M1 = empirical_mig[, 2], M2 = empirical_noMig[, 2]))
		}else{
			rob = c(rob, robust(x = pNoMig[j], M1 = empirical_noMig[, 1], M2 = empirical_mig[, 1]))
		}
	}
	names_tmp = colnames(x)
	x = cbind(x, pNoMig, pMig, rob)
	colnames(x) = c(names_tmp, "pNoMigration", "pMigration", "robustness")
	write.table(x, file = paste("summary_3_noPAN_", i, "_robustness.txt", sep=""), row.names = F, col.names = T, quote = F)
}

# loop over datasets simulated with current isolation 
for(i in noMig){
	print(i)
	x = read.table(paste("summary_", i, ".txt", sep=""), h=T)
	SI = grep("pSI", colnames(x))
	AM = grep("pAM", colnames(x))
	IM = grep("pIM", colnames(x))
	SC = grep("pSC", colnames(x))
#	PAN = grep("pPAN", colnames(x))

	pMig = apply(x[, c(IM, SC)], MARGIN=1, FUN="sum")
#	pMig = apply(x[, c(IM, SC, PAN)], MARGIN=1, FUN="sum")
	pNoMig = apply(x[, c(AM, SI)], MARGIN=1, FUN="sum")
	pMig2 = pMig/(pMig+pNoMig)
	pNoMig2 = pMig/(pMig+pNoMig)
	pMig = pMig2
	pNoMig = pNoMig2

	rob = NULL
	for(j in 1:length(pNoMig)){
		if(pNoMig[j] > pMig[j]){
			rob = c(rob, robust(x = pNoMig[j], M1 = empirical_noMig[, 1], M2 = empirical_mig[, 1]))
		}else{
			rob = c(rob, robust(x = pMig[j], M1 = empirical_mig[, 2], M2 = empirical_noMig[, 2]))
		}
	}
	names_tmp = colnames(x)
	x = cbind(x, pNoMig, pMig, rob)
	colnames(x) = c(names_tmp, "pNoMigration", "pMigration", "robustness")
	write.table(x, file = paste("summary_3_noPAN_", i, "_robustness.txt", sep=""), row.names = F, col.names = T, quote = F)
}


### test
plot(density(empirical_mig[,1], width = 0.15), xlim=c(0, 1), ylim=c(0,0.1))
lines(density(empirical_noMig[,1], width = 0.15), col="red")

plot(density(empirical_noMig[,2], width = 0.15), xlim=c(0, 1), ylim=c(0,0.1))
lines(density(empirical_mig[,2], width = 0.15), col="red")


