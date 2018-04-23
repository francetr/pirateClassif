#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

import sys
from Bio import SeqIO # For the fasta reading
import readInput

############
#	Help command : scriptClassif -h
#   Commandes to launch this script :
#	   python3 <path/toward/this/script> <path/toward/the/classif/file> <path/toward/the/baseline/file> <path/toward/the/fasta/file>
# example	   python3 code/scriptClassif.py sortie_pastec/TisoRepet1.classif base_reference.txt hat.fasta
#	   or  ./<path/toward/the>/scriptClassif.py <path/toward/the/classif/file> <path/toward/the/baseline/file> <path/toward/the/fasta/file>
############

# @author: Tristan Frances

def main():
	print("Start of the classification\n")
	####	 Instanciation of dictionnaries that will contain the results
	nonTE, potentialChimeric,  noCat, TE = {}, {}, {}, {}

	####	 retrieve and check the arguments used in command line
	print("####	Recovery of the arguments\n")
	args = readInput.retrieveArguments();
	# print("####	Valid extension of the files\n")
	pastecFile=args.classif
	baselineFile=args.baseline
	fastaFile=args.fasta
	####	Reading of the baseline file ####
	try:
		# take the second argument in the terminal command as file name
		print("####	Import of the Baseline file")
		baseline=readInput.readBaseline(baselineFile)
		# print(baseline)
		print("####	Baseline succesfully imported\n")
	except IndexError:
		print("/!\	Error: The baseline file provided is badly written\n####	Classification aborted")
		sys.exit(1)

	####	Reading of the classif file ####
	try:
		print("####	Read of the classif file")
		pastec=readInput.readPastec(pastecFile, nonTE, potentialChimeric, noCat, TE, baseline)
		# print(TE)
		print("####	End of reading of the classif file\n")
	except IndexError:
		print("/!\	Error: No PASTEC provided\n####	Classification aborted")
		sys.exit(1)

	# for seq in TE.keys():
	# 	print(seq, TE[seq])

	# print("TE : %d, noCat : %d, nonTE : %d, chimeric: %d, total: %d" %(len(TE), len(noCat), len(nonTE), len(potentialChimeric), (len(TE)+len(noCat)+len(nonTE)+len(potentialChimeric))))
	# print(TE, len(TE))
	# print(noCat)

	###	Reading of the fasta file ####
	# try:
	# 	print("####	Read of the FASTA file")
	# 	fasta=readFasta(fastaFile)
	# 	print("####	End of reading of the FASTA file\n")
	# 	for key in fasta.keys():
	# 		print(">%s :\n%s\n" %(key, fasta[key]["seq"]))
	# except IndexError:
	# 	print("####	No Fasta provided\n####	Classification aborted")
	# 	sys.exit(1)
	# save(fasta, nonTE, potentialChimeric, noCat, TE)

def save(FASTA, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""

	Save the sequence in a file.

	Return a dictionnary which contains the different catagories which caracterize the sequence.

	Keyword arguments:
	@type FASTA: string
	@param FASTA: name of the FASTA file containing the sequence that will be opened.
	@type NONTE: dictionnary
	@param NONTE: dictionnary for non transposable element (nonTE).
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@type NOCAT: dictionnary
	@param NOCAT: dictionnary for non categorized element (noCat).
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).

	@rtype: TODO
	"""
	#TODO
	# for sequenceName in TE:
	#	 if FASTA.id!=sequenceName:
	#	 print(FASTA.id, sequenceName)
	# return




####################
#	   MAIN
####################

if __name__ == "__main__":
	main()
