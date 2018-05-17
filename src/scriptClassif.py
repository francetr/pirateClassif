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
	####	 Instanciation of a dictionnary that will contain 4 dictionnaries (4 keys, one for each categories: nonTE, TE, potentialChimeric and noCat) that will contain the results
	seqClassified = {"nonTE" : {}, "potentialChimeric" : {},  "noCat" : {}, "TE" : {}}

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
		pastec=readInput.readPastec(pastecFile, seqClassified, baseline)
		print("####	End of reading of the classif file\n")
	except IndexError:
		print("/!\	Error: Error with the PASTEC file provided {}\n####	Classification aborted".format(pastecFile))
		sys.exit(1)

	# # for seq in sequencesClassified["TE"].keys():
	# # 	print(seq, sequencesClassified["TE"][seq])
	#
	# print("TE : %s, noCat : %s, nonTE : %s, chimeric: %s" %(sequencesClassified["TE"], sequencesClassified["noCat"], sequencesClassified["nonTE"], sequencesClassified["potentialChimeric"]))
	# print("TE : %d, noCat : %d, nonTE : %d, chimeric: %d, total: %d" %(len(sequencesClassified["TE"]), len(sequencesCategorizednoCat), len(sequencesClassified["nonTE"]), len(sequencesClassified["potentialChimeric"]), (len(sequencesClassified["TE"])+len(sequencesClassified["noCat"])+len(sequencesClassified["nonTE"])+len(sequencesClassified["potentialChimeric"]))))
	#
	# print(sequencesClassified["TE"], len(sequencesClassified["TE"]))
	# print(sequencesClassified["noCat"])
	#
	####	Reading of the fasta file ####
	# try:
	# 	print("####	Read of the FASTA file")
	# 	fasta=readInput.readFasta(fastaFile)
	# 	print("####	End of reading of the FASTA file\n")
	# 	# for key in fasta.keys():
	# 	# 	print(">%s :\n%s\n" %(key, fasta[key]["seq"]))
	# except IndexError:
	# 	print("/!\ Error: Error with the Fasta file provided {}\n####	Classification aborted".format(fastaFile))
	# 	sys.exit(1)
	# save.save(fasta, sequencesClassified)





####################
#	   MAIN
####################

if __name__ == "__main__":
	main()
