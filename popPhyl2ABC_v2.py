#!/usr/bin/python
from Bio.SeqIO import parse
from numpy.random import choice
import sys
#import os
from os import path
from os import popen
from os import remove
from os import system

#fastaFile = "Abatus_3sp.macse_DNA"
#nameA = "Abatus_cordatus"
#nameB = "Tripylus_abatoides" 
#mu = 0.00000002763
#rhovertheta = 0.5 
#N = 100000
#ind1 = "GA30F"
#ind2 = "GA30H"
#ind3 = "GA30J"
#ind4 = "GA30L"
#./popPhyl2ABC.py Chelonoidis_2vs2.fas Chelonoidis_porteri_vandenburghi Chelonoidis_becki 0.00000002763 0.5 100000

if(len(sys.argv) != 11):
	print('Command line:\npopPhyl2ABC_v2.py\tfasta_from_popPhyl\tspeciesA\tspeciesB\tmu\trho_over_theta\tN\tname_ind1\tname_ind2\tname_ind3\tname_ind4\n')
	print('Exemple:\npopPhyl2ABC_v2.py Abatus_3sp.macse_DNA Abatus_cordatus Tripylus_abatoides 0.00000002763 0.5 100000 GA30F GA30H GA30J GA30L')
	exit(0)


lociSizeLimit = 30
numberIndLimit = 4

fastaFile = sys.argv[1]
nameA = sys.argv[2]
nameB = sys.argv[3]
mu = float(sys.argv[4])
rhovertheta = float(sys.argv[5])
N = int(sys.argv[6])
ind1 = sys.argv[7]
ind2 = sys.argv[8]
ind3 = sys.argv[9]
ind4 = sys.argv[10]

input = parse(fastaFile, "fasta")

seq = {}

for i in input:
	if ind1 in i.id or ind2 in i.id or ind3 in i.id or ind4 in i.id:
		name = i.id.split("|")[0]
		if name not in seq:
			seq[name] = {}
			seq[name]["spA"] = []
			seq[name]["spB"] = []
			seq[name]["Ltot"] = len(i.seq)
			seq[name]["msLike"] = ""
		if name in seq and nameA in i.id:
			seq[name]["spA"].append(i.seq)
		if name in seq and nameB in i.id:
			seq[name]["spB"].append(i.seq)
input.close()

# clean the alignment by removing loci with few individuals
toRemove = []
for i in seq:
	if len(seq[i]['spA']) < numberIndLimit or len(seq[i]['spB']) < numberIndLimit:
		toRemove.append(i)

if len(toRemove) != 0:
	for i in toRemove:
		seq.pop(i)


for i in seq.keys():
	# prepare input for polydNdS
	tmp = ""
	cnt = -1
	for j in seq[i]['spA']:
		cnt += 1
		tmp += ">{0}\n{1}\n".format("ind_" + str(cnt), j)
	for j in seq[i]['spB']:
		cnt += 1
		tmp += ">{0}\n{1}\n".format("ind_" + str(cnt), j)
	output = open("tmp.fas", "w")
	output.write(tmp)
	output.close()
	# launch polydNdS
	polydNdS = popen("polydNdS -i tmp.fas -P -N | grep 'Mean # of synonymous' | awk '{print $NF}'")
	# get the synonymous length
	Lsyno_tmp = polydNdS.read()
	if(Lsyno_tmp==""):
		seq[i]["nInd"] = -9
		seq[i]["nSegSites"] = -9
		seq[i]["toRemove"] = 1
		continue
	seq[i]["Lsyno"] = int(float(Lsyno_tmp.strip()))
	# remove loci with few synonymous sites
	if(seq[i]["Lsyno"] <= lociSizeLimit):
		seq[i]["nInd"] = -9
		seq[i]["nSegSites"] = -9
		seq[i]["toRemove"] = 1
		continue
	# only treat loci with a minimum synonymous length lociSizeLimit
	if(seq[i]["Lsyno"] > lociSizeLimit):
		seq[i]["toRemove"] = 0
		# read polydNdS' output file
		input = open("tmp.fas.synonymous", "r")
		# first line
		j = input.readline().strip().split("\t") # nind \t nSegsites\n
		seq[i]["nInd"] = int(j[0])
		seq[i]["nSegSites"] = int(j[1])
		# only treat loci with at least 8 sequenced haplotype
		if seq[i]["nInd"] == 8:
			# second line: positions of segregating sites
			j = input.readline().strip().split("\t")
			j = [ int(k) for k in j ]
			j = [ round(k/(max(j)*1.0), 4) for k in j ]
			# if segregating sites
			if seq[i]["nSegSites"] > 0:
				seq[i]["msLike"] = "// 78 1.58854 1.58854 111 56 22 0.795745 1.02097 0.574558 1.49299 4.90703 4.90703 6.22109\nsegsites: {0}\npositions: {1}\n".format(seq[i]["nSegSites"], "\t".join([ str(k) for k in j ]))
				# remove the first line full of '?'
				j = input.readline()
				# first haplotype is treated as the ancestral genotype (only 0s)
				ancestral = input.readline().strip().split("\t")[1::]
				seq[i]["msLike"] += seq[i]["nSegSites"] * "0" + "\n"
				# treat the other individuals
				for j in input:
					j = j.strip().split("\t")[1::]
					for k in range(seq[i]["nSegSites"]):
						if j[k] == ancestral[k]:
							seq[i]["msLike"] += "0"
						if j[k] != ancestral[k]:
							seq[i]["msLike"] += "1"
					seq[i]["msLike"] += "\n"
		filesToDelete = ["tmp.fas", "tmp.fas.3rdpositions", "tmp.fas.4fold", "tmp.fas.all_silent", "tmp.fas.exons", "tmp.fas.introns_flanking", "tmp.fas.replacement", "tmp.fas.synonymous"]
		for j in filesToDelete:
			if path.isfile(j):
				remove(j)


for i in seq:
	if seq[i]['nSegSites'] == 0:
		seq[i]['msLike'] = "// 78 1.58854 1.58854 111 56 22 0.795745 1.02097 0.574558 1.49299 4.90703 4.90703 6.22109\nsegsites: 0\n\n"
		seq[i]['toRemove'] = 1 # if we want to remove with NO POLYMORPHISM AND NO DIVERGENCE

# remove the small / badly covered loci (too many Ns)
toRemove = []
for i in seq:
	if seq[i]['toRemove'] == 1:
		toRemove.append(i)

for i in toRemove:
	seq.pop(i)


# if there are more than 900 loci:
if len(seq) > 900:
	selectedLoci = choice(seq.keys(), 900, replace=False) # randomly select 900 loci
else:
	selectedLoci = seq.keys() # keep all loci


# prepare the output:
bpfile_length = []
bpfile_nA = []
bpfile_nB = []
bpfile_mu = []
bpfile_rho = []

locus_ms = "./msnsam tbs 20 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 1 2 -eN tbs tbs\n3579 27011 59243\n\n"

spinput = "\n{0}\n".format(len(selectedLoci))

for i in selectedLoci:
	loc = seq[i]
	# bpfile
	bpfile_length.append(loc['Lsyno'])
	bpfile_nA.append(len(loc['spA']))
	bpfile_nB.append(len(loc['spB']))
	bpfile_mu.append(4 * N * mu * loc['Lsyno'])
	bpfile_rho.append(4 * N * mu * loc['Lsyno'] * rhovertheta)
	# spinput
	spinput += "{0}\n{1}\n{2}\n".format(len(loc['spA']), len(loc['spB']), loc['Lsyno'])
	# msfile
	locus_ms += loc['msLike']
	locus_ms += "\n"


spinput_observation = spinput + "1\nlocus.ms\n"
spinput_simulation = spinput + "100000\nmyfifo\n"


bpfile = "#\t{0}\t{1}\t{2}\t{3}\n".format(nameA, nameB, N, mu)
bpfile += "\t".join([ str(k) for k in bpfile_length] ) + "\n"
bpfile += "\t".join([ str(k) for k in bpfile_nA] ) + "\n"
bpfile += "\t".join([ str(k) for k in bpfile_nB] ) + "\n"
bpfile += "\t".join([ str(k) for k in bpfile_mu] ) + "\n"
bpfile += "\t".join([ str(k) for k in bpfile_rho] ) + "\n"

output = open("bpfile", "w")
output.write(bpfile)
output.close()

output = open("spinput.txt", "w")
output.write(spinput_observation)
output.close()

output = open("locus.ms", "w")
output.write(locus_ms)
output.close()

system("mscalc")

output = open("spinput.txt", "w")
output.write(spinput_simulation)
output.close()


