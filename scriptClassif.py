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

	@rtype: argParse
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
	It first search if the beginning of line contains the character > typical of non specific keyword.
	The keywords are then added to the dictionnary

	Return a dictionnary which contains the sequence and the id of each sequences which are in the FASTA file.

	Keyword argument:
	@type BASELINE: string
	@param BASELINE: name of the baseline file that will be opened.

	@rtype: dictionnary
	@return: Dictionnary containing the different names of possible superFamily.
	{B{Key} = specific : {B{Key} = name of specific keywords for superfamily used later : B{I{Values}} = [list of possible names for this superfamily]},
	B{Key} = non specific :{B{Key} = name of non specific keyword for superfamily used later : B{I{Values}} = [list of possible names for this superfamily]}}.
	"""
	baselineDictionnary={"specific":{},"non specific":{}}
	try:
		with open(BASELINE, "r") as f:
			for line in f:
				####	List that will contains the different possibles names of a superfamily
				listPossibleNames=[]
				####	check the content of the file and look for > character (It indicates if the keyword is specific or non specific)
				keywords=re.search(r'>([^\n]+)', line)
				####	For specific keywords
				if not keywords:
					####	Recover the reference name for the superfamily
					superFamilyNames = line.split("\t")
					####	Recover the possible names for the superfamily
					allPossibleNames = superFamilyNames[1].strip().split(":")
					for possibleName in allPossibleNames:
						####	Add the possibles names for a superFamily in the list
						listPossibleNames.append(possibleName)
					####	Complete the dictionnary with : Key = name of superfamily used later; Value = list of possible names for this superfamily
					baselineDictionnary["specific"][superFamilyNames[0]]=listPossibleNames
				####	For non specific keyword
				else:
					####	Recover the reference name for the superfamily
					superFamilyNames = keywords.groups()[0].split("\t")
					####	Recover the possible names for the superfamily
					allPossibleNames = superFamilyNames[1].strip().split(":")
					for possibleName in allPossibleNames:
						####	Add the possibles names for a superFamily in the list
						listPossibleNames.append(possibleName)
					####	Complete the dictionnary with : Key = name of superfamily used later; Value = list of possible names for this superfamily
					baselineDictionnary["non specific"][superFamilyNames[0]]=listPossibleNames				# ####	List that will contains the different possibles names of a superfamily

	except FileNotFoundError:
		print("/!\	Error: No such file as {}\n####	Classification aborted".format(BASELINE))
	return baselineDictionnary

def readFasta(FASTA):
	"""
	Read a FASTA file as input and return it as a string. Use the SeqIO from the package Bio of BioPython.

	Return a dictionnary which contains the sequence and the id of each sequences that are in the FASTA file.

	Keyword argument:
	@type FASTA: string
	@param FASTA: name of the FASTA file containing the sequence that will be opened.

	@rtype: dictionnary
	@return: Dictionnary with {B{key} = id of the sequence and value : B{I{value}} = {"seq" : sequence of the FASTA sequence} }
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


def classDetermination(FEATURES, NONTE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
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
	if FEATURES[3] != "PotentialChimeric":
		####	Check the class of the sequence if it is not chimeric : class I
		if FEATURES[4] == "I" :
			TE[FEATURES[0]]={"class":"I"}
			# print(FEATURES[0])
			orderDetermination(FEATURES, POTENTIALCHIMERIC, NOCAT, TE, BASELINE)
		####	classII
		elif FEATURES[4] == "II":
			TE[FEATURES[0]]={"class":"II"}
			orderDetermination(FEATURES, POTENTIALCHIMERIC, NOCAT, TE, BASELINE)
			# print(FEATURES[0])
		####	NoCat
		elif FEATURES[4] == "noCat":
			NOCAT[FEATURES[0]]={"class":"unknown"}
			# print("NoCat : %s %s" % (FEATURES[4], FEATURES[5]))
			# pass
		####	 NonTE
		else:
			NONTE[FEATURES[0]]={"class":"nonTE"}
			# print("NA : %s %s" % (FEATURES[4], FEATURES[5]))
	else:
		POTENTIALCHIMERIC[FEATURES[0]]={"Class":"potentialChimeric"}
		# print("chimere : %s %s %s" % (FEATURES[0], FEATURES[4], FEATURES[5]))

def orderDetermination(FEATURES, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
	"""
	Determine the order of a sequence. Doing so it complete a dictionnary containing the transposable element.
	If the order wasn't found, the order's sequence is considered as unknown.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
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
	if FEATURES[5] == "noCat" :
		TE[FEATURES[0]]["order"]="unknown"
	####	 The order of the TE is a MITE, a LARD or a TRIM : the superfamily of the TE will not be determined
	elif FEATURES[5] == "MITE" or FEATURES[5] == "LARD" or FEATURES[5] == "TRIM" :
		TE[FEATURES[0]]["order"]=FEATURES[5]
	####	 The order of the TE is determined : the superfamily of the TE will be determined
	else:
		TE[FEATURES[0]]["order"]=FEATURES[5]
		superFamilyDetermination(FEATURES, POTENTIALCHIMERIC, NOCAT, TE, BASELINE)

def superFamilyDetermination(FEATURES, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
	"""
	Determine the superfamily of one sequence. Doing so it complete a dictionnary containing the transposable element.
	If the superFamily can't be defined, it will be unknown.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@type NOCAT: dictionnary
	@param NOCAT: dictionnary for non categorized element (noCat).
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@todo: Consider if POTENTIALCHIMERIC and NOCAT arguments are mandatory

	@return: None.
	"""
	#TODO deal with sequence that don't have DB comparisons
	try:
		####	 search if there is a 'coding' part in the 7th value of FEATURES => needed to defined the FEATURES superfamily. Search for coding(<anything that is not ));>); NB : regex (?!) anything that is not in ()
		codingRecord = re.search(r'coding=\(((?!\)\);).+?)\);', FEATURES[7]).groups()[0]
		# codingRecord = re.search(r'coding=\(([^\)]+?)\);?', FEATURES[7]).groups()[0] Another possibility of regex
		####	Split the coding part to obtain the different results according to the comparison with different databases (3 possibilties: TEBLRtx, TEBLRx and profiles)
		databaseRecords = codingRecord.split(';')
		####	Search if there are different names contained in the codingRecord
		searchDifferentName(FEATURES, TE, databaseRecords, BASELINE)

		# print(FEATURES[0], dbName, (matches[FEATURES[0]][dbName]))
	except AttributeError as e:
		####	If there is no coding part : declaration of superfamily as unknown
		TE[FEATURES[0]]["superfamily"] = "unknown"
		# NOCAT[FEATURES[0]]=TE[FEATURES[0]]
		# del TE[FEATURES[0]]
		return

def searchDifferentName(FEATURES, TE, DATABASERECORDS, BASELINE):
	"""
	If the superFamily can't be defined, it will be unknown.
	The name is recovered from two : the hmm profiles (profiles part) and REPET (TE_BLRx and TE_BLRtx part).

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type DATABASERECORDS: list
	@param DATABASERECORDS: list of string in which there are the superFamily name to search.
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@return: None.
	"""
	####	Dictionnary that will contains (or not) the different superFamilies find for the concerned sequence
	# matches={FEATURES[0]:{}}
	####	List to store the superfamily found for a sequence. Usefull to compare if there are differences between them
	superFamilyFound=[]
	####	Scan the different results obtain for each comparisons
	for dr in DATABASERECORDS:
		####	Split to obtain the name of the database used (TE_BLRx or TE_BLRtx)
		dbName = dr.split(':')[0].strip()
		# matches[FEATURES[0]] = {dbName:[]}
		####	if profiles is found in coding, search for superFamily name using function searchProfilesName
		if dbName=='profiles':
			searchProfilesName(FEATURES, dr, superFamilyFound)
		####	if TE_BLRx or TE_BLRtx is found in coding, search for superFamily name using function searchRepBaseName
		elif dbName=="TE_BLRx" or "TE_BLRtx":
			searchRepBaseName(FEATURES, dr, superFamilyFound)
	####	String that will contain the final supefamily name of the sequence
	finalSuperFamily=""
	####	if multiple names for superfamily are found tor the sequence, proceed to the comparison between all the names found
	if len(superFamilyFound) > 1:
		finalSuperFamily=superFamilyComparison(FEATURES, superFamilyFound, BASELINE)
	####	One superFamily name have been found
	elif len(superFamilyFound) == 1:
		####	Replace the superfamily name by "unknown"
		if superFamilyFound[0] == "?":
			finalSuperFamily = "unknown"
		else:
			finalSuperFamily=superFamilyFound[0]
	####	No superFamily name have been found
	else:
		finalSuperFamily="unknown"
	# matches[FEATURES[0]][dbName].append(searchObj.groups())
	TE[FEATURES[0]]["superFamily"]=finalSuperFamily

def searchProfilesName(FEATURES, DATABASERECORD, SUPERFAMILYFOUND):
	"""

	Search the keywords in the profiles part of coding.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type SUPERFAMILYFOUND: list
	@param SUPERFAMILYFOUND: names of the superfamily found for one sequence during the superFamilyDetermination.

	@return: TODO
	"""
	#TODO try to differenciate different possible syntax according the first regex found. Ex: PiRATEdb_CACTA_Tase_NA or _RT_reina_NA_RT_NA or PF05699.9_Dimer_Tnp_hAT_NA_Tase_21.4
	# print(DATABASERECORD)
	####	first split to get rid of profiles:
	DATABASERECORD=DATABASERECORD.split("profiles:")
	####	Parse all the results of profiles and search for regex
	for substr in DATABASERECORD[1].split(','):
		# print(substr+"\n")
		try:
			searchObj = re.search(r'(?:profiles:)? _?([^_]+)_([^_]+)_([^_]+)_([^:]+): ([0-9\.]+)%\(([0-9\.]+)%\)', substr)
			# print(FEATURES[0], searchObj.groups()[0])
			####	Proceed to different operations according the first string of the first regex
			####	If first regex is PiRATEdb
			if searchObj.groups()[0] == "PiRATEdb":
				print("Base de Données PiRATE : ", FEATURES[0], searchObj.groups()[1], searchObj.groups()[2])
			# print(searchObj.groups())
		except AttributeError:
			print('Issue during searching profiles on : '+ substr + "\n")


def searchRepBaseName(FEATURES, DATABASERECORD, SUPERFAMILYFOUND):
	"""

	Search the keywords in the TE_BLRx and TE_BLRtx part of coding.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
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
			# print(FEATURES[0], searchObj.groups()[1], searchObj.groups()[2], searchObj.groups()[3])
			####	TODO If a sequence superfamily haven't been determined.
			if searchObj.groups()[3] == "?":
				TE[FEATURES[0]]["superFamily"]="unknown"
		except AttributeError:
			print('Issue during searching RepBase name on : '+ substr)


def superFamilyComparison(FEATURES, SUPERFAMILYFOUND, BASELINE):
	"""

	Compare the different superFamilies found for one sequence. For this uses the method combinations of the package itertools allowing to combine 2 by 2 the different elements of a list.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type SUPERFAMILYFOUND: list
	@param SUPERFAMILYFOUND: names of the superfamily found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily.

	@rtype: string
	@return: name of the superfamily corresponding to the sequence.
	"""
	####	TODO Comparison with a reference base
	####	Dictionnary with key = name of superFamily : value = number of presence of this name in the SUPERFAMILYFOUND list
	superFamilyCount={}
	####	Draft
	####	Comparison between names found for a sequence and names in the BASELINE
	####	Run through the list of superFamily found
	for superFamily in SUPERFAMILYFOUND:
		####	Run through the reference name of the BASELINE
		for possibleName in BASELINE["specific"].keys():

			####	If the superfamily name is not in the BASELINE AND is not a key of the dictionnary count, its string is declared as value of this dictionnary and has its counter equal to 1
			if superFamily in BASELINE["specific"][possibleName] and (not(superFamily in superFamilyCount.keys())):
				superFamilyCount[superFamily]=1
			####	If the superfamily name is already in the dictionnary and exists in the BASELINE, its counter is incremented by 1
			elif (superFamily in BASELINE["specific"][possibleName]) and (superFamily in superFamilyCount.keys()):
				superFamilyCount[superFamily]+=1
	# print(FEATURES[0], superFamilyCount)

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

	# print(FEATURES[0], name, max)
	return name
	# print("{}, {}, count: {}, len: {}".format(FEATURES[0], SUPERFAMILYFOUND, superFamilyCount, len(SUPERFAMILYFOUND)))

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
