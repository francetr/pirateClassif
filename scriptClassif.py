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
	@type PASTEC: string
	@param PASTEC: name of the classif file that will be opened.
	@type NONTE: dictionnary
	@param NONTE: dictionnary for non transposable element (nonTE).
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@type NOCAT: dictionnary
	@param NOCAT: dictionnary for non categorized element (noCat).
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).
	@type BASELINE: dictionnary
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
	@type BASELINE: string
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

	Return a dictionnary which contains the sequence and the id of each sequences that are in the FASTA file.

	Keyword argument:
	@type FASTA: string
	@param FASTA: name of the FASTA file containing the sequence that will be opened.

	@return: Dictionnary with {key = id of the sequence and value : value = {"seq" : sequence of the FASTA sequence} }
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
	@type SEQUENCE: list
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed. Usefull values of this list :
		- [0] : id of the sequence
		- [3] : potentiel chimeric
		- [4] : class of the sequence
		- [5] : order of the sequence
		- [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
	@type NONTE: dictionnary
	@param NONTE: dictionnary for non transposable element (nonTE).
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@type NOCAT: dictionnary
	@param NOCAT: dictionnary for non categorized element (noCat).
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).
	@type BASELINE: dictionnary
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
	@type NONTE: dictionnary
	@param NONTE: dictionnary for non transposable element (nonTE).
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@type NOCAT: dictionnary
	@param NOCAT: dictionnary for non categorized element (noCat).
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None.
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
	@type SEQUENCE: list
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@type NOCAT: dictionnary
	@param NOCAT: dictionnary for non categorized element (noCat).
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None.
	"""
	####	 The order of the TE is not determined, the superfamily of the TE will not be determined
	if SEQUENCE[5] == "noCat" :
		TE[SEQUENCE[0]]["order"]="unknown"
	####	 The order of the TE is a MITE : the superfamily of the TE will not be determined
	elif SEQUENCE[5] == "MITE":
		TE[SEQUENCE[0]]["order"]=SEQUENCE[5]
	####	 The order of the TE is determined : the superfamily of the TE will be determined
	else:
		TE[SEQUENCE[0]]["order"]=SEQUENCE[5]
		superFamilyDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE)

def superFamilyDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
	"""
	Determine the superfamily of one sequence. Doing so it complete a dictionnary containing the transposable element.
	If the superFamily can't be defined, it will be unknown.
	@todo: Consider if POTENTIALCHIMERIC and NOCAT arguments are mandatory

	Keyword arguments:
	@type SEQUENCE: list
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@type NOCAT: dictionnary
	@param NOCAT: dictionnary for non categorized element (noCat).
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).
	@type BASELINE: dictionnary
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
		####	If there is no coding part : declaration of superfamily as unknown
		TE[SEQUENCE[0]]["superfamily"] = "unknown"
		# NOCAT[SEQUENCE[0]]=TE[SEQUENCE[0]]
		# del TE[SEQUENCE[0]]
		return

def searchDifferentName(SEQUENCE, TE, DATABASERECORDS, BASELINE):
	"""
	If the superFamily can't be defined, it will be unknown.
	The name is recovered from two : the hmm profiles (profiles part) and REPET (TE_BLRx and TE_BLRtx part).

	Keyword arguments:
	@type SEQUENCE: list
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@type DATABASERECORDS: list
	@param DATABASERECORDS: list of string in which there are the superFamily name to search.
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None.
	"""
	####	Dictionnary that will contains (or not) the different superFamilies find for the concerned sequence
	matches={SEQUENCE[0]:{}}
	####	List to store the superfamily found for a sequence. Usefull to compare if there are differences between them
	superFamilyFound=[]
	####	Scan the different results obtain for each comparisons
	for dr in DATABASERECORDS:
		####	Split to obtain the name of the database used (TE_BLRx or TE_BLRtx)
		dbName = dr.split(':')[0].strip()
		# matches[SEQUENCE[0]] = {dbName:[]}

		####	Do nothing for the profil hmm database. TODO implement the code to search regex in the case it exists
		if dbName=='profiles':
			searchProfilesName(SEQUENCE, dr, superFamilyFound)
		####	Split to obtain each comparison, and parse them in order to search interesting regex
		elif dbName=="TE_BLRx" or "TE_BLRtx":
			searchRepBaseName(SEQUENCE, dr, superFamilyFound)
	####	String that will contain the final supefamily name of the sequence
	finalSuperFamily=""
	####	if multiple names for superfamily are found tor the sequence, proceed to the comparison between all the names found
	if len(superFamilyFound) > 1:
		finalSuperFamily=superFamilyComparison(SEQUENCE, superFamilyFound, BASELINE)
	#### One superFamily name have been found
	elif len(superFamilyFound) == 1:
		####	Replace the superfamily name by "unknown"
		if superFamilyFound[0] == "?":
			finalSuperFamily = "unknown"
		else:
			finalSuperFamily=superFamilyFound[0]
	####	No superFamily name have been found
	else:
		finalSuperFamily="unknown"
	# matches[SEQUENCE[0]][dbName].append(searchObj.groups())
	TE[SEQUENCE[0]]["superFamily"]=finalSuperFamily

def searchProfilesName(SEQUENCE, DATABASERECORD, SUPERFAMILYFOUND):
	"""

	Search the keywords in the profiles part of coding.

	Keyword arguments:
	@type SEQUENCE: list
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type SUPERFAMILYFOUND: list
	@param SUPERFAMILYFOUND: names of the superfamily found for one sequence during the superFamilyDetermination.

	@return: None
	"""
	for substr in DATABASERECORD.split(','):
		try:
			searchObj = re.search(r' ([^:]+):(?:Class)?(I+|\?):([^:]+):([^:]+): ([0-9\.]+)%', substr)
			# print(searchObj)
		except AttributeError:
			print('Issue on : '+ substr)


def searchRepBaseName(SEQUENCE, DATABASERECORD, SUPERFAMILYFOUND):
	"""

	Search the keywords in the TE_BLRx and TE_BLRtx part of coding.

	Keyword arguments:
	@type SEQUENCE: list
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type SUPERFAMILYFOUND: list
	@param SUPERFAMILYFOUND: names of the superfamily found for one sequence during the superFamilyDetermination.

	@return: None
	"""
	for substr in DATABASERECORD.split(','):
		try:
			####	regex split into groups. Group's number correspondance : 0 : name of the sequence; 1 : type of class; 2 : type of order; 3 : superFamily name
			searchObj = re.search(r' ([^:]+):(?:Class)?(I+|\?):([^:]+):([^:]+): ([0-9\.]+)%', substr)
			SUPERFAMILYFOUND.append(searchObj.groups()[3])
			# print(SEQUENCE[0], searchObj.groups()[1], searchObj.groups()[2], searchObj.groups()[3])
			####	TODO If a sequence superfamily haven't been determined.
			if searchObj.groups()[3] == "?":
				TE[SEQUENCE[0]]["superFamily"]="unknown"
		except AttributeError:
			print('Issue on : '+ substr)


def superFamilyComparison(SEQUENCE, SUPERFAMILYFOUND, BASELINE):
	"""

	Compare the different superFamilies found for one sequence. For this uses the method combinations of the package itertools allowing to combine 2 by 2 the different elements of a list.

	Keyword arguments:
	@type SEQUENCE: list
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@type SUPERFAMILYFOUND: list
	@param SUPERFAMILYFOUND: names of the superfamily found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily.

	@return: name of the superfamily corresponding to the sequence
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
			max=superFamilyCount[superFamily]
			name=superFamily
		####	Else if the count is equal to the max
		elif superFamilyCount[superFamily] == max:
			max=superFamilyCount[superFamily]
			name="potentialChimeric"
	#### if the name is ? convert it into unknown
	if name=="?":
		name="unknown"

	# print(SEQUENCE[0], name, max)
	return name
	# print("{}, {}, count: {}, len: {}".format(SEQUENCE[0], SUPERFAMILYFOUND, superFamilyCount, len(SUPERFAMILYFOUND)))

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

	@return: TODO
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

	# for seq in TE.keys():
	# 	print(TE[seq])

	print("TE : %d, noCat : %d, nonTE : %d, chimeric: %d, total: %d" %(len(TE), len(noCat), len(nonTE), len(potentialChimeric), (len(TE)+len(noCat)+len(nonTE)+len(potentialChimeric))))
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
