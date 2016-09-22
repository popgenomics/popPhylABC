noMig = c("empirical_all_SI_NHetero", "empirical_all_SI_NHomo", "empirical_all_AM_MHetero_NHetero", "empirical_all_AM_MHetero_NHomo", "empirical_all_AM_MHomo_NHetero", "empirical_all_AM_MHomo_NHomo")
mig = c("empirical_all_IM_MHetero_NHetero", "empirical_all_IM_MHetero_NHomo", "empirical_all_IM_MHomo_NHetero", "empirical_all_IM_MHomo_NHomo", "empirical_all_SC_MHetero_NHetero", "empirical_all_SC_MHetero_NHomo", "empirical_all_SC_MHomo_NHetero", "empirical_all_SC_MHomo_NHomo", "empirical_all_PAN_NHetero", "empirical_all_PAN_NHomo")


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
	x = read.table(paste("empirical/", i, sep=""), h=T)
	empirical_noMig_tmp = rbind(empirical_noMig_tmp, x)
}

for(i in mig){
	x = read.table(paste("empirical/", i, sep=""), h=T)
	empirical_mig_tmp = rbind(empirical_mig_tmp, x)
}

SI = grep("SI", colnames(empirical_noMig_tmp))
AM = grep("AM", colnames(empirical_noMig_tmp))
IM = grep("IM", colnames(empirical_noMig_tmp))
SC = grep("SC", colnames(empirical_noMig_tmp))
PAN = grep("PAN", colnames(empirical_noMig_tmp))

empirical_noMig = cbind(apply(empirical_noMig_tmp[, c(SI, AM)], MARGIN=1, FUN="sum"), apply(empirical_noMig_tmp[, c(IM, SC, PAN)], MARGIN=1, FUN="sum")) # noMig mig
empirical_mig = cbind(apply(empirical_mig_tmp[, c(SI, AM)], MARGIN=1, FUN="sum"), apply(empirical_mig_tmp[, c(IM, SC, PAN)], MARGIN=1, FUN="sum")) # noMig mig

# treat the datasets
# loop over datasets simulated with ongoing migration
for(i in mig){
	i = paste(strsplit(i, "_")[[1]][-c(1,2)], collapse="_")
	print(i)
	x = read.table(paste("summary_", i, ".txt", sep=""), h=T)
	SI = grep("pSI", colnames(x))
	AM = grep("pAM", colnames(x))
	IM = grep("pIM", colnames(x))
	SC = grep("pSC", colnames(x))
	PAN = grep("pPAN", colnames(x))

	pMig = apply(x[, c(IM, SC, PAN)], MARGIN=1, FUN="sum")
	pNoMig = apply(x[, c(AM, SI)], MARGIN=1, FUN="sum")

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
	write.table(x, file = paste("summary_2_", i, "_robustness.txt", sep=""), row.names = F, col.names = T, quote = F)
}

# loop over datasets simulated with current isolation 
for(i in noMig){
	i = paste(strsplit(i, "_")[[1]][-c(1,2)], collapse="_")
	print(i)
	x = read.table(paste("summary_", i, ".txt", sep=""), h=T)
	SI = grep("pSI", colnames(x))
	AM = grep("pAM", colnames(x))
	IM = grep("pIM", colnames(x))
	SC = grep("pSC", colnames(x))
	PAN = grep("pPAN", colnames(x))

	pMig = apply(x[, c(IM, SC, PAN)], MARGIN=1, FUN="sum")
	pNoMig = apply(x[, c(AM, SI)], MARGIN=1, FUN="sum")

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
	write.table(x, file = paste("summary_2_", i, "_robustness.txt", sep=""), row.names = F, col.names = T, quote = F)
}

