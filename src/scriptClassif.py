#!/usr/bin/python2
# coding: utf8

""" @author: Tristan Frances """

import sys
import readInput
import save
from timeit import default_timer as timer

############
#	Help command : scriptClassif -h
#   Commandes to launch this script :
#	   python3 <path/toward/this/script> <path/toward/the/classif/file> <path/toward/the/fasta/file>
# example	   python3 code/scriptClassif.py ArabiTEdenovo.txt ArabiTEdenovo.fasta
#	   or  ./<path/toward/the>/scriptClassif.py <path/toward/the/classif/file> <path/toward/the/fasta/file>
############

def main():

	print("Start of the classification\n")
	####	Instanciation of a dictionnary with the id of the sequence as key and the values will contain 6 dictionnaries
	####	(6 keys, 1 for each categories: class, order, superFamily; 1 for the type of file in which the sequence will be saved (e.g. nonTE, TE,
	####	potentialChimeric and noCat); 1 for the log and 1 for the unknown_keywords founded during the search for superFamily name)
	seqClassified = {}

	####	 retrieve the arguments used in command line
	print("####	Recovery of the arguments\n")
	args = readInput.retrieveArguments();
	####	store the arguments into different variables
	pastecFile=args.classif
	print("The pastec file used is %s"%(pastecFile))
	fastaFile=args.fasta
	print("The FASTA file used is %s"%(fastaFile))
	identityThreshold=args.e
	print("The identity threshold used is {identity} %".format(identity=identityThreshold))
	baselineFile=args.baseline
	print("The baseline file used is %s\n"%(baselineFile))

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

	####	print the number sequences by saveType found and the number of unknown_keywords if founded
	cptTE, cptnonTE, cptChimeric, cptnoCat, cptError = 0, 0, 0, 0, 0
	for seq in seqClassified:
		# print(seqClassified[seq])
		if seqClassified[seq]["saveType"] == "TE":
			cptTE+=1
		if seqClassified[seq]["saveType"] == "nonTE":
			cptnonTE+=1
		if seqClassified[seq]["saveType"] == "potentialChimeric":
			cptChimeric+=1
		if seqClassified[seq]["saveType"] == "noCat":
			cptnoCat+=1
		if seqClassified[seq]["unknown_keyword"] != "\n":
			cptError+=1

	print("Number total of sequences found : %d" %((cptTE + cptnonTE + cptnoCat + cptChimeric)))
	print("Number total of TE : %d; with recognized TE : %d, uncategorized : %d, potentialChimeric : %d"%((cptTE + cptnoCat + cptChimeric), cptTE, cptnoCat, cptChimeric))
	print("Number of sequences found for nonTE : %d" %(cptnonTE))
	print("Number of sequence with unknown keywords found: %s\n"%(cptError))

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
	start = timer()
	main()
	end = timer()
	print("\nThe execution of the script has taken %s sec"%(end - start))
