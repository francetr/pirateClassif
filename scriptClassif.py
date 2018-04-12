#!/usr/bin/python3
# coding: utf8
import sys
from Bio import SeqIO # For the fasta reading
import argparse
import re


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
	Retrieve the arguments from the command line.

	@return: attributes NameSpace of object argParse that contain the two arguments name.

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
# 	@param CLASSIFNAME: name of the argument for the CLASSIF file
# 	@param FASTANAME: name of the argument for the FASTA file
#
# 	@return: None
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



def readPastec(PASTEC, NONTE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
	"""
	Read a PASTEC file as input line by line, and then proceed to the categorization of each sequence.

	Keyword arguments:
	@param PASTEC: name of the classif file that will be opened.
	@param NONTE: dictionnary for non transposable element (nonTE).
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@param NOCAT: dictionnary for non categorized element (noCat).
	@param TE: dictionnary for transposable element (I or II).
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None.

	"""
	try:
		####	Open the classif file
		with open(PASTEC, "r") as f:
			print("Begin of the categorization\n")
			####	Parse every line/sequence of the file
			for line in f:
				sequence=line.replace("\t\t", "\t").strip() # remove the first column double tab and carriage return of the line
				categorization(sequence, NONTE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE)
				# print("Summary, sequences find by type :\nTE : {}\nnonTE : {}\nnoCat : {}\npotentialChimeric : {}\nTotal : {}".format(len(TE),\
				# len(nonTE), len(noCat), len(potentialChimeric), (len(TE)+len(nonTE)+len(noCat)+len(potentialChimeric))))
			print("End of the categorization\n")

	except FileNotFoundError:
		####	Prevent the opening if the file name is incorrect
		print("/!\	Error: No such file as {}\n####	Classification aborted".format(PASTEC))
		sys.exit(1)

def readBaseline(BASELINE):
	"""
	Read a PASTEC file as input line by line, and then proceed to the categorization of each sequence.

	Return a dictionnary which contains the sequence and the id of each sequences which are in the FASTA file.

	Keyword argument:
	@param BASELINE: name of the baseline file that will be opened.

	@return: Dictionnary containing the different names of possible superFamily. {Key = name of superfamily used later : Values = [list of possible names for this superfamily]}.
	"""
	baselineDictionnary={}
	try:
		with open(BASELINE, "r") as f:
			for line in f:
				####	List that will contains the different possibles names of a superfamily
				listPossibleNames=[]
				####	Recover the reference name for the superfamily
				superFamilyNames = line.split("\t")
				####	Recover the possible names for the superfamily
				allPossibleNames = superFamilyNames[1].strip().split(":")
				for possibleName in allPossibleNames:
					####	Add the possibles names for a superFamily in the list
					listPossibleNames.append(possibleName)
				####	Complete the dictionnary with : Key = name of superfamily used later; Value = list of possible names for this superfamily
				baselineDictionnary[superFamilyNames[0]]=listPossibleNames
	except FileNotFoundError:
		print("/!\	Error: No such file as {}\n####	Classification aborted".format(BASELINE))
	return(baselineDictionnary)

def readFasta(FASTA):
	"""
	Read a FASTA file as input and return it as a string. Use the SeqIO from the package Bio of BioPython.

	Return a dictionnary which contains the sequence and the id of each sequences which are in the FASTA file.

	Keyword argument:
	@param FASTA: name of the FASTA file containing the sequence that will be opened.

	@return: Dictionnary with {key = id of the sequence and value : value = {"seq" : sequence of the FASTA sequence}}
		- id : id of the sequence
		- name : name of the sequence
		- description : description of the sequence
		- number of features : number of features of the sequence
		- seq : sequence concerned

	"""
	seqReturned={}
	try:
		####	 Open the fasta file
		with open(fastaFile, "rU") as handle:
			print(handle)
		####	Parse every line/sequence of the file
			for record in SeqIO.parse(handle, "fasta"):
				####	Save the sequence in a dictionnary
				seqReturned[record.id]={"seq":record.seq}
		return seqReturned
	except (FileNotFoundError, NameError):
		####	prevent the opening if the file name is incorrect
		print("/!\	Error: No such file as {}\n####	Classification aborted".format(FASTA))
		sys.exit(1)

def categorization(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
	"""
	Categorize the Transposable Element from the input argument.

	Keyword arguments:
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed. Usefull values of this list :
		- [0] : id of the sequence
		- [3] : potentiel chimeric
		- [4] : class of the sequence
		- [5] : order of the sequence
		- [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
	@param NONTE: dictionnary for non transposable element (nonTE).
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@param NOCAT: dictionnary for non categorized element (noCat).
	@param TE: dictionnary for transposable element (I or II).
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None.

	"""
	features=SEQUENCE.split("\t")
	#TODO Treatment of the superfamily
	classDetermination(features, nonTE, potentialChimeric,  noCat, TE, BASELINE)


def classDetermination(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
	"""
	Determine the class of a sequence.
	For this, complete four dictionnaries, passed onto arguments, that will contain the different catagories that caracterize the sequence.

	Keyword arguments:
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@param NONTE: dictionnary for non transposable element (nonTE).
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@param NOCAT: dictionnary for non categorized element (noCat).
	@param TE: dictionnary for transposable element (I or II).
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None

	"""
	####	Check first if the sequence is chimeric
	if SEQUENCE[3] != "PotentialChimeric":
		####	Check the class of the sequence if it is not chimeric : class I
		if SEQUENCE[4] == "I" :
			TE[SEQUENCE[0]]={"class":"I"}
			# print(SEQUENCE[0])
			orderDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE)
			pass
		####	classII
		elif SEQUENCE[4] == "II":
			TE[SEQUENCE[0]]={"class":"II"}
			orderDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE)
			# print(SEQUENCE[0])
			pass
		####	NoCat
		elif SEQUENCE[4] == "noCat":
			NOCAT[SEQUENCE[0]]={"class":"unknown"}
			# print("NoCat : %s %s" % (SEQUENCE[4], SEQUENCE[5]))
			# pass
		####	 NonTE
		else:
			NONTE[SEQUENCE[0]]={"class":"nonTE"}
			# print("NA : %s %s" % (SEQUENCE[4], SEQUENCE[5]))
			pass
	else:
		POTENTIALCHIMERIC[SEQUENCE[0]]={"Class":"potentialChimeric"}
		# print("chimere : %s %s %s" % (SEQUENCE[0], SEQUENCE[4], SEQUENCE[5]))

def orderDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
	"""
	Determine the order of a sequence. Doing so it complete a dictionnary containing the transposable element.
	If the order wasn't found, the order's sequence is considered as unknown.

	Keyword arguments:
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@param NOCAT: dictionnary for non categorized element (noCat).
	@param TE: dictionnary for transposable element (I or II).
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None.

	"""
	if SEQUENCE[5] == "noCat" :
		####	 The order of the TE is not determined, the superfamily of the TE will not be determined
		TE[SEQUENCE[0]]["order"]="unknown"
	elif SEQUENCE[5] == "MITE":
		####	 The order of the TE is a MITE : the superfamily of the TE will not be determined
		TE[SEQUENCE[0]]["order"]=SEQUENCE[5]
	else:
		####	 The order of the TE is determined : the superfamily of the TE will be determined
		TE[SEQUENCE[0]]["order"]=SEQUENCE[5]
		superFamilyDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE)

def superFamilyDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
	"""
	Determine the superfamily of one sequence. Doing so it complete a dictionnary containing the transposable element.
	If the superFamily can't be defined, it will be unknown.
	@todo: Consider if POTENTIALCHIMERIC and NOCAT arguments are mandatory

	Keyword arguments:
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@param NOCAT: dictionnary for non categorized element (noCat).
	@param TE: dictionnary for transposable element (I or II).
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None.

	"""
	#TODO deal with sequence that don't have DB comparisons
	try:
		####	 search if there is a 'coding' part in the 7th value of SEQUENCE => needed to defined the SEQUENCE superfamily
		codingRecord = re.search(r'coding=\(([^\)]+)\)', SEQUENCE[7]).groups()[0]
		# print(codingRecord)
		####	Split the coding part to obtain the different results according to the comparison with different databases (3 possibilties: TEBLRtx, TEBLRx and profiles)
		databaseRecords = codingRecord.split(';')
		####	Search if there are different names contained in the codingRecord
		searchDifferentName(SEQUENCE, TE, databaseRecords, BASELINE)

		# print(SEQUENCE[0], dbName, (matches[SEQUENCE[0]][dbName]))
	except AttributeError:
		####	If there is no coding part : do nothing (the order has been previously defined as unknown)
		# print(SEQUENCE[7])
		# NOCAT[SEQUENCE[0]]=TE[SEQUENCE[0]]
		# del TE[SEQUENCE[0]]
		return

def searchDifferentName(SEQUENCE, TE, DATABASERECORDS, BASELINE):
	"""
	If the superFamily can't be defined, it will be unknown.
	The name is recovered from two : the hmm profiles (profiles part) and REPET (TE_BLRx and TE_BLRtx part).

	Keyword arguments:
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@param TE: dictionnary for transposable element (I or II) which will be completed with superfamily name.
	@param DATABASERECORDS: list containing the different results for the different database (TE_BLRx, TE_BLRtx, profiles).
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None.

	"""
	####	Dictionnary that will contains (or not) the different superFamilies find for the concerned sequence
	matches={SEQUENCE[0]:{}}
	superFamilyFound=[]
	####	Scan the different results obtain for each comparisons
	for dr in DATABASERECORDS:
		####	Split to obtain the name of the database used (TE_BLRx or TE_BLRtx)
		dbName = dr.split(':')[0].strip()
		matches[SEQUENCE[0]] = {dbName:[]}
		####	List to store the superfamily found for a sequence. Usefull to compare if there are differences between them


		####	Do nothing for the profil hmm database. TODO implement the code to search regex in the case it exists
		if dbName=='profiles': continue
		####	Split to obtain each comparison, and parse them in order to search interesting regex
		# TODO compare the superFamily find for each substr, !!!!!! troubles if the sequence don't have DB but just profiles => put in noCat
		for substr in dr.split(','):
			try:
				####	regex split into groups. Group's number correspondance : 0 : name of the sequence; 1 : type of class; 2 : type of order; 3 : superFamily name
				searchObj = re.search(r' ([^:]+):(?:Class)?(I+|\?):([^:]+):([^:]+): ([0-9\.]+)%', substr)
				superFamilyFound.append(searchObj.groups()[3])
				# print(SEQUENCE[0], searchObj.groups()[1], searchObj.groups()[2], searchObj.groups()[3])
				####	TODO If a sequence superfamily haven't been determined.
				if searchObj.groups()[3] == "?":
					TE[SEQUENCE[0]]["superFamily"]="unknown"
			except AttributeError:
				print('Issue on : '+ substr)
	# if multiple names for superfamily are found tor the sequence, proceed to the comparison between all the names found
	if len(superFamilyFound) > 1:
		superFamilyComparison(SEQUENCE, superFamilyFound, BASELINE)
	# matches[SEQUENCE[0]][dbName].append(searchObj.groups())
	# TE[SEQUENCE[0]]["superFamily"]=searchObj.groups()[3]


def superFamilyComparison(SEQUENCE, SUPERFAMILYFOUND, BASELINE):
	"""

	Compare the different superFamilies found for one sequence. For this uses the method combinations of the package itertools allowing to combine 2 by 2 the different elements of a list.

	Return a dictionnary which contains the different catagories which caracterize the sequence.

	Keyword arguments:
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@param SUPERFAMILYFOUND: list of the name of the superfamilies found for one sequence during the superFamilyDetermination.
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily.

	"""
	####	TODO Comparison with a reference base
	####	Dictionnary with key = name of superFamily : value = number of presence of this name in the SUPERFAMILYFOUND list
	superFamilyCount={}
	####	Draft
	####	Comparison between names found for a sequence and names in the BASELINE
	####	Run through the list of superFamily found
	for superFamily in SUPERFAMILYFOUND:
		####	Run through the reference name of the BASELINE
		for possibleName in BASELINE.keys():
			####	If the superfamily name is not in the BASELINE AND is not a key of the dictionnary count, its string is declared as value of this dictionnary and has its counter equal to 1
			if superFamily in BASELINE[possibleName] and (not(superFamily in superFamilyCount.keys())):
				superFamilyCount[superFamily]=1
			####	If the superfamily name is already in the dictionnary and exists in the BASELINE, its counter is incremented by 1
			elif (superFamily in BASELINE[possibleName]) and (superFamily in superFamilyCount.keys()):
				superFamilyCount[superFamily]+=1
	# print(SEQUENCE[0], superFamilyCount)

	name=""
	max = 0
	####	Parse the keys of the superfamily name dictionnary count
	for superFamily in superFamilyCount.keys():
		####	If the count of superfamily name is superior to max, max is equal to this count and name is the superFamily name
		if superFamilyCount[superFamily] > max:
			name=superFamily
			max=superFamilyCount[superFamily]
		####	Else if the count is equal to the max
		elif superFamilyCount[superFamily] == max:
			max=superFamilyCount[superFamily]
			name="potentialChimeric"
	# print(SEQUENCE[0], name, max)

	print("{}, {}, count: {}, len: {}".format(SEQUENCE[0], SUPERFAMILYFOUND, superFamilyCount, len(SUPERFAMILYFOUND)))

def save(FASTA, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""

	Save the sequence in a file.

	Return a dictionnary which contains the different catagories which caracterize the sequence.

	Keyword arguments:
	@param FASTA: name of the list of strings, which contain the nucleotid of the sequences.
	@param NONTE: dictionnary for non transposable element (nonTE).
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@param NOCAT: dictionnary for non categorized element (noCat).
	@param TE: dictionnary for transposable element (I or II).

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
	print("####	Recovery of the arguments\n")
	args = retrieveArguments();
	# print("####	Valid extension of the files\n")
	pastecFile=args.classif
	baselineFile=args.baseline
	fastaFile=args.fasta
	####	Reading of the baseline file ####
	try:
		# take the second argument in the terminal command as file name
		print("####	Import of the Baseline file")
		baseline=readBaseline(baselineFile)
		# print(baseline)
		print("####	Baseline succesfully imported\n")
	except IndexError:
		print("/!\	Error: The baseline file provided is badly written\n####	Classification aborted")
		sys.exit(1)

	####	Reading of the classif file ####
	try:
		print("####	Read of the classif file")
		pastec=readPastec(pastecFile, nonTE, potentialChimeric, noCat, TE, baseline)
		# print(TE)
		print("####	End of reading of the classif file\n")
	except IndexError:
		print("/!\	Error: No PASTEC provided\n####	Classification aborted")
		sys.exit(1)
	# for referenceName in baseline.keys():
	# 	for i in range(len(baseline[referenceName])):
	# 		print(referenceName, baseline[referenceName][i])

	print("TE : %d, noCat : %d, nonTE : %d, chimeric: %d, total: %d" %(len(TE), len(noCat), len(nonTE), len(potentialChimeric), (len(TE)+len(noCat)+len(nonTE)+len(potentialChimeric))))
	# print(TE, len(TE))
	# print(noCat)

	# for key in TE.keys():
	# 	if (TE[key]["class"]!="PotentialChimeric"):
	# 		print(TE[key])
	# 	pass

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
