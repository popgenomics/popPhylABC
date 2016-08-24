SI = read.table("SI", h=T)

AM = read.table("AM", h=T)

IM = read.table("IM", h=T)

SC = read.table("SC", h=T)

PAN = read.table("PAN", h=T)

par(mfrow=c(2,2))
plot(density(c(0, 1, SI$IM + SI$SC + SI$PAN)), main="SI")
plot(density(c(0, 1, AM$IM + AM$SC + AM$PAN)), main="AM")
plot(density(c(0, 1, IM$IM + IM$SC + IM$PAN)), main="IM")
plot(density(c(0, 1, SC$IM + SC$SC + SC$PAN)), main="SC")

mig_given_mig = c(IM$IM + IM$SC + IM$PAN, SC$IM + SC$SC + SC$PAN)
mig_given_nomig = c(SI$IM + SI$SC + SI$PAN, AM$IM + AM$SC + AM$PAN)
nomig_given_mig = c(IM$SI + IM$AM , SC$SI + SC$AM, PAN$SI + PAN$AM)
nomig_given_nomig = c(SI$SI + SI$AM , AM$SI + AM$AM)


# robustness when the target is SI for 2,000 PODs
res = NULL
for(i in 1:nrow(SI)){
	pmig = SI$IM[i] + SI$SC[i] + SI$PAN[i]
	pnomig = SI$SI[i] + SI$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = F
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = T
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

SI = cbind(SI, res)


# robustness when the target is AM for 2,000 PODs
res = NULL
for(i in 1:nrow(AM)){
	pmig = AM$IM[i] + AM$SC[i] + AM$PAN[i]
	pnomig = AM$SI[i] + AM$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = F
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = T
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

AM = cbind(AM, res)


# robustness when the target is IM for 2,000 PODs
res = NULL
for(i in 1:nrow(IM)){
	pmig = IM$IM[i] + IM$SC[i] + IM$PAN[i]
	pnomig = IM$SI[i] + IM$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = T
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = F
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

IM = cbind(IM, res)


# robustness when the target is SC for 2,000 PODs
res = NULL
for(i in 1:nrow(SC)){
	pmig = SC$IM[i] + SC$SC[i] + SC$PAN[i]
	pnomig = SC$SI[i] + SC$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = T
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = F
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

SC = cbind(SC, res)

# robustness when the target is PAN for 2,000 PODs
res = NULL
for(i in 1:nrow(PAN)){
	pmig = PAN$IM[i] + PAN$SC[i] + PAN$PAN[i]
	pnomig = PAN$SI[i] + PAN$AM[i]
	if(pmig>pnomig){
		A = density(mig_given_mig, from=pmig, to=pmig)$y[1]
		B = density(mig_given_nomig, from=pmig, to=pmig)$y[1]
		good_inference = T
	}else{
		A = density(nomig_given_nomig, from=pnomig, to=pnomig)$y[1]
		B = density(nomig_given_mig, from=pnomig, to=pnomig)$y[1]
		good_inference = F
	}
	res = rbind(res, c(good_inference, A/(A+B)))
}
colnames(res) = c("correctly_infered", "robustness")

PAN = cbind(PAN, res)

