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
	####	 Instanciation of dictionnaries that will contain the results
	nonTE, potentialChimeric,  noCat, TE = {}, {}, {}, {}

	####	 retrieve and check the arguments used in command line
	print("####	Recovery of the arguments\n")
	args = readInput.retrieveArguments();
	# print("####	Valid extension of the files\n")
	pastecFile=args.classif
	fastaFile=args.fasta
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
		pastec=readInput.readPastec(pastecFile, nonTE, potentialChimeric, noCat, TE, baseline)
		print("####	End of reading of the classif file\n")
	except IndexError:
		print("/!\	Error: Error with the PASTEC file provided {}\n####	Classification aborted".format(pastecFile))
		sys.exit(1)

	# # for seq in TE.keys():
	# # 	print(seq, TE[seq])
	#
	# print("TE : %s, noCat : %s, nonTE : %s, chimeric: %s" %(TE, noCat, nonTE, potentialChimeric))
	# print("TE : %d, noCat : %d, nonTE : %d, chimeric: %d, total: %d" %(len(TE), len(noCat), len(nonTE), len(potentialChimeric), (len(TE)+len(noCat)+len(nonTE)+len(potentialChimeric))))
	#
	# print(TE, len(TE))
	# print(noCat)
	#
	# ##	Reading of the fasta file ####
	# try:
	# 	print("####	Read of the FASTA file")
	# 	fasta=readInput.readFasta(fastaFile)
	# 	print("####	End of reading of the FASTA file\n")
	# 	# for key in fasta.keys():
	# 	# 	print(">%s :\n%s\n" %(key, fasta[key]["seq"]))
	# except IndexError:
	# 	print("/!\ Error: Error with the Fasta file provided {}\n####	Classification aborted".format(fastaFile))
	# 	sys.exit(1)
	# save.save(fasta, nonTE, potentialChimeric, noCat, TE)





####################
#	   MAIN
####################

if __name__ == "__main__":
	main()
