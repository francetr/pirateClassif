#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

import sys
import readInput
import save

############
#	Help command : scriptClassif -h
#   Commandes to launch this script :
#	   python3 <path/toward/this/script> <path/toward/the/classif/file> <path/toward/the/baseline/file> <path/toward/the/fasta/file>
# example	   python3 code/scriptClassif.py sortie_pastec/TisoRepet1.classif base_reference.txt hat.fasta
#	   or  ./<path/toward/the>/scriptClassif.py <path/toward/the/classif/file> <path/toward/the/baseline/file> <path/toward/the/fasta/file>
############

def main():
	print("Start of the classification\n")
	####	Instanciation of a dictionnary that will contain 5 dictionnaries (5 keys, one for each categories: nonTE, TE, potentialChimeric and noCat)
	####	that will contain the results and a log
	seqClassified = {"nonTE" : {}, "potentialChimeric" : {},  "noCat" : {}, "TE" : {}, "log":{}}

	####	 retrieve the arguments used in command line
	print("####	Recovery of the arguments\n")
	args = readInput.retrieveArguments();
	####	store the arguments into different variables
	pastecFile=args.classif
	fastaFile=args.fasta
	identityThreshold=args.e
	baselineFile=args.baseline
	####	Reading of the baseline file ####
	try:
		print("####	Import of the Baseline file")
		baseline=readInput.readBaseline(baselineFile)
		print("####	Baseline succesfully imported\n")
	except IndexError:
		print("/!\	Error: The baseline file {} provided is badly written\n####	Classification aborted".format(baselineFile))
		sys.exit(1)

	####	Reading of the classif file ####
	try:
		print("####	Read of the classif file")
		pastec=readInput.readPastec(pastecFile, seqClassified, baseline, identityThreshold)
		print("####	End of reading of the classif file\n")
	except IndexError:
		print("/!\	Error: Error with the PASTEC file provided {}\n####	Classification aborted".format(pastecFile))
		sys.exit(1)

	# for seq in seqClassified["TE"].keys():
	# 	print(seq, seqClassified["TE"][seq])
	#
	# print("TE : %s, noCat : %s, nonTE : %s, chimeric: %s" %(seqClassified["TE"], seqClassified["noCat"], seqClassified["nonTE"], seqClassified["potentialChimeric"]))
	print("TE : %d, noCat : %d, nonTE : %d, chimeric: %d, total: %d" %(len(seqClassified["TE"]), len(seqClassified["noCat"]), len(seqClassified["nonTE"]), len(seqClassified["potentialChimeric"]), (len(seqClassified["TE"])+len(seqClassified["noCat"])+len(seqClassified["nonTE"])+len(seqClassified["potentialChimeric"]))))
	#
	# print(seqClassified["TE"], len(seqClassified["TE"]))
	# print(seqClassified["noCat"])
	#
	####	Reading of the fasta file ####
	try:
		print("####	Read of the FASTA file")
		fasta=readInput.readFasta(fastaFile)
		print("####	End of reading of the FASTA file\n")
	except IndexError:
		print("/!\ Error: Error with the Fasta file provided {}\n####	Classification aborted".format(fastaFile))
		sys.exit(1)
	save.save(fasta, seqClassified)


####################
#	   MAIN
####################

if __name__ == "__main__":
	main()
