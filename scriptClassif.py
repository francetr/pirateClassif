#!/usr/bin/python3.5
# coding: utf8
import sys
from Bio import SeqIO # For the fasta reading
import argparse
import re
import itertools


############
#	Help command : scriptClassif -h
#   Commandes to launch this script :
#	   python3 <path/toward/this/script> <path/toward/the/classif/file> <path/toward/the/baseline/file> <path/toward/the/fasta/file>
# example	   python3 code/scriptClassif.py sortie_pastec/TisoRepet1.classif base_reference.txt hat.fasta
#	   or  ./<path/toward/the>/scriptClassif.py <path/toward/the/classif/file> <path/toward/the/baseline/file> <path/toward/the/fasta/file>
############

# @author: Tristan Frances

def retrieveArguments():
	"""
	Retrieve the arguments from the command line

	:return args: attributes NameSpace of object argParse that contain the two arguments name

	"""
	####	 Mananage the 2 arguments (PASTEC and FASTA file name) when the command is launched
	parser = argparse.ArgumentParser(prog="scriptClassif.py", description="This program is a part of the PiRATE project. It aims to automatized the step of TE classification")
	parser.add_argument("classif", type=str, help="classif file that comes from PASTEC")
	parser.add_argument("baseline", type=str, help="baseline file giving the different names possible for a superfamily")
	parser.add_argument("fasta", type=str, help="fasta file providing the sequence")
	args = parser.parse_args()
	# checkArguments(args.classif, args.fasta)
	return(args)

# def checkArguments(CLASSIFNAME, FASTANAME):
# 	"""
# 	Check the extension of the two arguments of command line. If it is not correct, stop the programm.
#
# 	Keyword argument
# 	:param CLASSIFNAME: name of the argument for the CLASSIF file
# 	:param FASTANAME: name of the argument for the FASTA file
#
# 	:return: None
#
# 	"""
# 	try:
# 		classifName=re.match(r'[\S]*[/]?[\w\-.]+(classif)$',CLASSIFNAME).groups()[0] ####	regex checking is classif file as good extension
# 	except AttributeError as e:
# 		print("####	Wrong extension for the classif file\n####	 Classification aborted")
# 		sys.exit(1)
#
# 	try:
# 		fastaName=re.match(r'[\S]*[/]?[\w\-.]+(fasta)$', FASTANAME).groups()[0] ####	regex checking is fasta file as good extension
# 	except AttributeError as e:
# 		print("####	Wrong extension for the fasta file\n####	Classification aborted")
# 		sys.exit(1)



def readPastec(PASTEC, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""
	Read a PASTEC file as input line by line, and then proceed to the categorization of each sequence

	Keyword argument
	:param PASTEC: name of the classif file that will be opened
	:param NONTE: dictionnary for non transposable element (nonTE)
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	:return: None

	"""
	try:
		####	Open the classif file
		with open(PASTEC, "r") as f:
			print("Begin of the categorization\n")
			####	Parse every line/sequence of the file
			for line in f:
				sequence=line.replace("\t\t", "\t").strip() # remove the first column double tab and carriage return of the line
				categorization(sequence, NONTE, POTENTIALCHIMERIC, NOCAT, TE)
				# print("Summary, sequences find by type :\nTE : {}\nnonTE : {}\nnoCat : {}\npotentialChimeric : {}\nTotal : {}".format(len(TE),\
				# len(nonTE), len(noCat), len(potentialChimeric), (len(TE)+len(nonTE)+len(noCat)+len(potentialChimeric))))
			print("End of the categorization\n")

	except FileNotFoundError:
		####	Prevent the opening if the file name is incorrect
		print("####	No such file as {}\n####	Classification aborted".format(PASTEC))
		sys.exit(1)

def readBaseline(BASELINE):
	"""
	Read a PASTEC file as input line by line, and then proceed to the categorization of each sequence

	Keyword argument
	:param BASELINE: name of the baseline file that will be opened.

	:return baselineDictionnary: Dictionnary containing the different names of possible superFamily. Key = name of superfamily used later; Values = list of possible names for this superfamily

	"""
	baselineDictionnary={}
	try:
		with open(BASELINE, "r") as f:
			for line in f:
				####	List that will contains the different possibles names of a superfamily
				listPossibleNames=[]
				####	Retrieve the reference name for the superfamily
				superFamilyNames = line.split("\t")
				####	Retrieve the possible names for the superfamily
				possibleNames = superFamilyNames[1].strip().split(":")
				for possibleName in possibleNames:
					####	Add the possibles names for a superFamily in the list
					listPossibleNames.append(possibleName)
				####	Complete the dictionnary with : Key = name of superfamily used later; Value = list of possible names for this superfamily
				baselineDictionnary[superFamilyNames[0]]=listPossibleNames
	except FileNotFoundError:
		print("####	No such file as {}\n####	Classification aborted".format(BASELINE))
	return(baselineDictionnary)

def readFasta(FASTA):
	"""
	Read a FASTA file as input and return it as a string. Use the SeqIO from the package Bio of BioPython

	Return a dictionnary with keys:
	id : id of the sequence
	name : name of the sequence
	description : description of the sequence
	number of features : number of features of the sequence
	seq : sequence concerned

	Keyword argument
	:param FASTA: name of the FASTA file containing the sequence that will be opened

	:return seqReturned: Dictionnary with key: id of the sequence and value : {"seq" : sequence of the FASTA sequence}

	"""
	seqReturned={}
	try:
		####	 Open the fasta file
		with open(fastaFile, "rU") as handle:
		####	Parse every line/sequence of the file
			for record in SeqIO.parse(handle, "fasta"):
				####	Save the sequence in a dictionnary
				seqReturned[record.id]={"seq":record.seq}
		return seqReturned
	except (FileNotFoundError, NameError):
		####	prevent the opening if the file name is incorrect
		print("####	No such file as {}\n####	Classification aborted".format(FASTA))
		sys.exit(1)

def categorization(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""
	Categorize the Transposable Element from the input argument

	Return a dictionnary which contains the different catagories which caracterize the sequence.

	Keyword argument
	:param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed. Usefull values of this list :
	- [0] : id of the sequence
	- [3] : potentiel chimeric
	- [4] : class of the sequence
	- [5] : order of the sequence
	- [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
	:param NONTE: dictionnary for non transposable element (nonTE)
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	:return: None

	"""
	features=SEQUENCE.split("\t")
	#TODO Treatment of the superfamily
	classDetermination(features, nonTE, potentialChimeric,  noCat, TE)


def classDetermination(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""
	Determine the class of a sequence.
	For this, complete four dictionnaries, passed onto arguments, that will contain the different catagories that caracterize the sequence.

	Keyword argument
	:param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed
	:param NONTE: dictionnary for non transposable element (nonTE)
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	:return: None

	"""
	####	Check first if the sequence is chimeric
	if SEQUENCE[3] != "PotentialChimeric":
		####	Check the class of the sequence if it is not chimeric : class I
		if SEQUENCE[4] == "I" :
			TE[SEQUENCE[0]]={"class":"I"}
			# print(SEQUENCE[0])
			orderDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE)
			pass
		####	classII
		elif SEQUENCE[4] == "II":
			TE[SEQUENCE[0]]={"class":"II"}
			orderDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE)
			# print(SEQUENCE[0])
			pass
		####	NoCat
		elif SEQUENCE[4] == "noCat":
			NOCAT[SEQUENCE[0]]={"class":"noCat"}
			# print("NoCat : %s %s" % (SEQUENCE[4], SEQUENCE[5]))
			pass
		####	 NonTE
		else:
			NONTE[SEQUENCE[0]]={"class":"nonTE"}
			# print("NA : %s %s" % (SEQUENCE[4], SEQUENCE[5]))
			pass
	else:
		POTENTIALCHIMERIC[SEQUENCE[0]]={"Class":"potentialChimeric"}
		# print("chimere : %s %s %s" % (SEQUENCE[0], SEQUENCE[4], SEQUENCE[5]))

def orderDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE):
	"""
	Determine the order of a sequence. Doing so it complete a dictionnary containing the transposable element.
	If the order wasn't found, the order's sequence is considered as unknown.

	Keyword argument
	:param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	:return: None

	"""
	if SEQUENCE[5] == "noCat" or "NA":
		TE[SEQUENCE[0]]["order"]="unknown"
	else:
		TE[SEQUENCE[0]]["order"]=SEQUENCE[5]
	####	 The order of the TE is determined : the super family of the TE will be determined
	superFamilyDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE)

def superFamilyDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE):
	"""
	Determine the super family of one sequence. Doing so it complete a dictionnary
	containing the transposable element.
	If the superFamily can't be defined, it will be unknown
	If there is 2 or more superfamily possible, the superFamily sequence will be considered as POTENTIALCHIMERIC

	:param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	:return: None

	"""
	#TODO deal with sequence that don't have DB comparisons
	try:
		####	 search if there is a 'coding' part in the 7th value of SEQUENCE => needed to defined the SEQUENCE superfamily
		codingRecord = re.search(r'coding=\(([^\)]+)\)', SEQUENCE[7]).groups()[0]
		# print(codingRecord)
		####	Split the coding part to obtain the different results according to the comparison with different databases (3 possibilties: TEBLRtx, TEBLRx and profiles)
		databaseRecords = codingRecord.split(';')
		####	Dictionnary that will contains (or not) the different superFamilies find for the concerned sequence
		matches={SEQUENCE[0]:{}}
		####	Scan the different results obtain for each comparisons
		for dr in databaseRecords:
			####	Split to obtain the name of the database used (TE_BLRx or TE_BLRtx)
			dbName = dr.split(':')[0].strip()
			matches[SEQUENCE[0]] = {dbName:[]}
			####	List to store the super family found for a sequence. Usefull to compare if there are differences between them
			superFamilyFound=[]

			####	Do nothing for the profil hmm database. TODO implement the code to search regex in the case it exists
			if dbName=='profiles': continue
			####	Split to obtain each comparison, and parse them in order to search interesting regex
			# TODO compare the superFamily find for each substr, !!!!!! troubles if the sequence don't have DB but just profiles => put in noCat
			for substr in dr.split(','):
				try:
					searchObj = re.search(r' ([^:]+):Class(I+):([^:]+):([^:]+): ([0-9\.]+)%', substr)
					superFamilyFound.append(searchObj.groups()[3])
					####	TODO If a sequence haven't been detemined, put it in NOCAT dictionnary and remove it from the TE. !!!! Troubles : don't take into account ? can be found in multiple time for 1 sequence
					# if searchObj.groups()[3] == "?":
						# TE[SEQUENCE[0]]["superFamily"]="unknown"
				except AttributeError:
					print('Issue on : '+ substr)
			superFamilyComparison(SEQUENCE, superFamilyFound)
			matches[SEQUENCE[0]][dbName].append(searchObj.groups())
			TE[SEQUENCE[0]]["superFamily"]=searchObj.groups()[3]

			# print(SEQUENCE[0], dbName, (matches[SEQUENCE[0]][dbName]))
	except AttributeError:
		####	If there is no coding part : the sequence is put in the NOCAT dictionnary and removed from the TE dictionnary
		# print(SEQUENCE[7])
		NOCAT[SEQUENCE[0]]=TE[SEQUENCE[0]]
		del TE[SEQUENCE[0]]
		return

def superFamilyComparison(SEQUENCE, SUPERFAMILYFOUND):
	"""

	Compare the different superFamilies found for one sequence. For this uses the method combinations of the package itertools allowing to combine 2 by 2 the different elements of a list.

	Return a dictionnary which contains the different catagories which caracterize the sequence.

	Keyword argument
	:param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed
	:param SUPERFAMILYFOUND: list of the the superfamilies found for one sequence during the superFamilyDetermination

	"""
	####	TODO Comparison with a reference base
	####	Assignation of 2 counter, 1 for the identical superfamily and the other for the different superfamily
	cptIdentique = 0; cptDifferent = 0
	####	Draft
	####	Construction of the different pair of family possible
	for familyA, familyB in itertools.combinations(SUPERFAMILYFOUND, 2):
		if SEQUENCE[0] == "TEDENOVO_Arabidopsis_thaliana1723" and familyA != familyB:
			# print("%s : A = %s different de B= %s" %(SEQUENCE[0], familyA, familyB))
			cptDifferent+=1
		elif SEQUENCE[0] == "TEDENOVO_Arabidopsis_thaliana1723" and familyA == familyB :
			# print("%s : A = %s identique a B = %s" %(SEQUENCE[0], familyA, familyB))
			cptIdentique+=1

	if SEQUENCE[0] == "TEDENOVO_Arabidopsis_thaliana1723":
		# print("different: {}, identique: {}".format(cptDifferent, cptIdentique))
		pass

def save(FASTA, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""

	Save the sequence in a file

	Return a dictionnary which contains the different catagories which caracterize the sequence.

	Keyword argument
	:param FASTA: name of the list of strings, which contain the nucleotid of the sequences
	:param NONTE: dictionnary for non transposable element (nonTE)
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

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
	# nom="/home/tfrances/Bureau/donnees/sortie_PASTEC/TisoMgescan.classif"

	print("Start of the classification\n")
	####	 Instanciation of dictionnaries that will contain the results
	nonTE, potentialChimeric,  noCat, TE = {}, {}, {}, {}

	####	 retrieve and check the arguments used in command line
	print("####	Check the extensions of the files")
	args = retrieveArguments();
	print("####	Valid extension of the files\n")
	pastecFile=args.classif
	baselineFile=args.baseline
	####	Reading of the classif file ####
	try:
		# take the second argument in the terminal command as file name
		print("####	Read of the Baseline file")
		baseline=readBaseline(baselineFile)
		print(baseline)
		print("####	End of reading of the Baseline file\n")
		print("####	Read of the classif file")
		pastec=readPastec(pastecFile, nonTE, potentialChimeric, noCat, TE)
		# print(TE)
		print("####	End of reading of the classif file\n")
	except IndexError:
		print("####	No PASTEC provided\n####	Classification aborted")
		sys.exit(1)


	print("TE : %d, noCat : %d, nonTE : %d, chimeric: %d" %(len(TE), len(noCat), len(nonTE), len(potentialChimeric)))
	# print(TE)

	# for key in TE.keys():
	#	 if (TE[key]["class"]!="PotentialChimeric"):
	#	 print(TE[key])
	#	 pass
	#
	####	Reading of the fasta file ####
	# try:
	#	 print("####	Read of the FASTA file")
	#	 fastaFile=sys.argv[2] # take the second argument in the terminal command as file name
	#	 fasta=readFasta(fastaFile)
	#	 print("####	End of reading of the FASTA file\n")
	# except IndexError:
	#	 print("####	No Fasta provided\n####	Classification aborted")
	#	 sys.exit(1)
	# for key in fasta.keys():
	#	 print(">%s :\n%s\n" %(key, fasta[key]["seq"]))
	# save(fasta, nonTE, potentialChimeric, noCat, TE)
