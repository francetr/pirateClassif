#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

import re
import search

def initCategorization(SEQUENCE, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD):
	"""
	Categorize the Transposable Element from the input argument.

	Keyword arguments:
	@type SEQUENCE: list
	@param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed. Usefull values of this list :
		- [0] : id of the sequence
		- [1] : length
		- [3] : potentiel chimeric
		- [4] : class of the sequence
		- [5] : order of the sequence
		- [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 6 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 3 : for the results (class, order and superFamily)
		- 1 : for the classification summary (classification_summary)
		- 1 : for the unknown keyword (unknown_keyword)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: float
	@param IDENTITYTHRESHOLD: threshold for which the superFamily of a sequence will be named, if its frequency in the sequence is greater or equal to this threshold.

	@rtype: None
	"""
	features=SEQUENCE.split("\t")
	####	Instanciation of the key for the sequence which will contain 6 dictionnaries: 4 for the results (TE, nonTE, potentialChimeric and noCat), 1 for the type of file save and one for the log
	SEQCLASSIFIED[features[0]]={}
	classDetermination(features, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD)
	####	Add a carriage return at the end of the sequence log
	finalDegreeClassification(features, SEQCLASSIFIED)
	SEQCLASSIFIED[features[0]]["unknown_keyword"]+="\n"

def finalDegreeClassification(FEATURES, SEQCLASSIFIED):
	"""
	Write a final degree of classification for a sequence in the Classification summary.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: name of the list of strings, which contain the sequence to categorize, that will be parsed. Usefull values of this list :
		- [0] : id of the sequence
		- [1] : length
		- [3] : potentiel chimeric
		- [4] : class of the sequence
		- [5] : order of the sequence
		- [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 6 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 3 : for the results (class, order and superFamily)
		- 1 : for the classification summary (classification_summary)
		- 1 : for the unknown keyword (unknown_keyword)

	@rtype: None
	"""
	finalDegree=""
	####	Case potentialChimeric
	if SEQCLASSIFIED[FEATURES[0]]["saveType"] == "potentialChimeric":
		####	Case the sequence has been found as chimeric in the beginning
		if not "order" in SEQCLASSIFIED[FEATURES[0]]:
			finalDegree = SEQCLASSIFIED[FEATURES[0]]["class"]
		####	Case the sequence has been found as chimeric in the beginning
		else:
			finalDegree = SEQCLASSIFIED[FEATURES[0]]["class"]

	####	Case uncategorized
	elif SEQCLASSIFIED[FEATURES[0]]["saveType"] == "noCat":
		finalDegree = SEQCLASSIFIED[FEATURES[0]]["class"]
	####	Case non TE
	elif SEQCLASSIFIED[FEATURES[0]]["saveType"] == "nonTE":
		finalDegree = SEQCLASSIFIED[FEATURES[0]]["order"]
	####	Case TE
	else:
		####	No superFamily have been found
		if not "superFamily" in SEQCLASSIFIED[FEATURES[0]]:
			finalDegree = SEQCLASSIFIED[FEATURES[0]]["order"]
		####	A superFamily has been found
		else:
			finalDegree = SEQCLASSIFIED[FEATURES[0]]["superFamily"]
	SEQCLASSIFIED[FEATURES[0]]["log"] += "\tFINALDEGREE:\t%s\n"%(finalDegree)

def classDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD):
	"""
	Determine the class of a sequence.
	For this, complete four dictionnaries, passed onto arguments, that will contain the different catagories that caracterize the sequence.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: name of the list of strings, which contain the sequence to categorize, that will be parsed. Usefull values of this list :
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 6 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 3 : for the results (class, order and superFamily)
		- 1 : for the classification summary (classification_summary)
		- 1 : for the unknown keyword (unknown_keyword)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: float
	@param IDENTITYTHRESHOLD: threshold for which the superFamily of a sequence will be named, if its frequency in the sequence is greater or equal to this threshold.

	@rtype: None
	"""
	####	Boolean if the order need to be search
	searchOrder=False
	####	Check first if the sequence is chimeric
	if FEATURES[3] != "PotentialChimeric":
		####	Check the class of the sequence if it is not chimeric : class I or II
		if FEATURES[4] == "I" or FEATURES[4] == "II" :
			SEQCLASSIFIED[FEATURES[0]]={"saveType":"TE","class":FEATURES[4]}
			searchOrder=True
		####	NoCat
		elif FEATURES[4] == "noCat":
			SEQCLASSIFIED[FEATURES[0]]={"saveType":"noCat","class":"undefined"}
		####	 NonTE
		else:
			SEQCLASSIFIED[FEATURES[0]]={"saveType":"nonTE","class":"nonTE"}
			searchOrder=True
	else:
		SEQCLASSIFIED[FEATURES[0]]={"saveType":"potentialChimeric","class":"potentialChimeric"}
	####	Initialize a log for concerned sequences that will be put into a log file and add the length, the type and the class of the sequence
	SEQCLASSIFIED[FEATURES[0]]["log"]=str("{name}\t{length}\t{saveType}\t{seqClass}\t".format(name=FEATURES[0], length=FEATURES[1], saveType=SEQCLASSIFIED[FEATURES[0]]["saveType"], seqClass=SEQCLASSIFIED[FEATURES[0]]["class"]))
	####	Initialize a unknown_keyword for concerned sequences that will be put into an unknown_keyword file. It concerned unknown_keywords found during search of BLAST and PROTPROFILES keywords
	SEQCLASSIFIED[FEATURES[0]]["unknown_keyword"]=str("")
	####	The sequence is a TE so order need to be determined
	if searchOrder:
		orderDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD)
		# print(SEQCLASSIFIED[FEATURES[0]]["unknown_keyword"])


def orderDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD):
	"""
	Determine the order of a sequence. Doing so it complete a dictionnary containing the transposable element.
	If the order wasn't found, the order's sequence is considered as unknown.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: name of the list of strings, which contain the sequence to categorize, that will be parsed. Usefull values of this list :
		- [0] : id of the sequence
		- [1] : length
		- [3] : potentiel chimeric
		- [4] : class of the sequence
		- [5] : order of the sequence
		- [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 6 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 3 : for the results (class, order and superFamily)
		- 1 : for the classification summary (classification_summary)
		- 1 : for the unknown keyword (unknown_keyword)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: float
	@param IDENTITYTHRESHOLD: threshold for which the superFamily of a sequence will be named, if its frequency in the sequence is greater or equal to this threshold.

	@rtype: None
	"""
	####	Boolean to know if the superFamily need to be search
	searchSuperFamily=False
	####	 The order of the TE is not determined, the superfamily of the TE will not be determined
	if FEATURES[5] == "noCat" :
		SEQCLASSIFIED[FEATURES[0]]["order"]="undefined"
	####	 The order of the TE is a MITE, a LARD or a TRIM : the superfamily of the TE will not be determined
	elif FEATURES[5] == "MITE" or FEATURES[5] == "LARD" or FEATURES[5] == "TRIM" or FEATURES[5] == "PotentialHostGene" or FEATURES[5] == "SSR" :
		SEQCLASSIFIED[FEATURES[0]]["order"]=FEATURES[5]
	####	 The order of the TE is determined : the superfamily of the TE will be determined
	else:
		SEQCLASSIFIED[FEATURES[0]]["order"]=FEATURES[5]
		searchSuperFamily=True
	####	Add the order of the sequence onto the log file
	SEQCLASSIFIED[FEATURES[0]]["log"]+=str("{order}\t".format(order=SEQCLASSIFIED[FEATURES[0]]["order"]))

	####	The superFamily need to be determined
	if searchSuperFamily:
		superFamilyDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD)

def superFamilyDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD):
	"""
	Determine the superfamily of one sequence. Doing so it complete a dictionnary containing the transposable element.
	If the superFamily can't be defined, it will be unknown.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: name of the list of strings, which contain the sequence to categorize, that will be parsed. Usefull values of this list :
		- [0] : id of the sequence
		- [1] : length
		- [3] : potentiel chimeric
		- [4] : class of the sequence
		- [5] : order of the sequence
		- [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 6 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 3 : for the results (class, order and superFamily)
		- 1 : for the classification summary (classification_summary)
		- 1 : for the unknown keyword (unknown_keyword)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: float
	@param IDENTITYTHRESHOLD: threshold for which the superFamily of a sequence will be named, if its frequency in the sequence is greater or equal to this threshold.

	@rtype: None
	"""
	try:
		####	 search if there is a 'coding' part in the 7th value of FEATURES => needed to defined the FEATURES superfamily. Search for coding(<anything that is not ));>); NB : regex (?!) anything that is not in ()
		codingRecord = re.search(r'coding=\(((?!\)\);).+?)\);', FEATURES[7]).groups()[0]
		####	Split the coding part to obtain the different results according to the comparison with different databases (3 possibilties: TEBLRtx, TEBLRx and profiles)
		databaseRecords = codingRecord.split(';')
		####	Search if there are different names contained in the codingRecord
		finalSuperFamilyName = search.searchDifferentName(FEATURES, SEQCLASSIFIED, databaseRecords, BASELINE, IDENTITYTHRESHOLD)
		####	Associate the superFamily name with the corresponding sequence into the right dictionnary
		associateSuperFamily(FEATURES, SEQCLASSIFIED, finalSuperFamilyName)

	except AttributeError:
		####	If there is no coding part : declaration of superfamily as undefined
		SEQCLASSIFIED[FEATURES[0]]["superfamily"] = "undefined"
		return

def associateSuperFamily(FEATURES, SEQCLASSIFIED, FINALSUPERFAMILYNAME):
	"""
	Associate the superFamily name for the sequence in the corresponding dictionnary (potentialChimeric if find potentialChimeric in FINALSUPERFAMILYNAME; else it will be TE)

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: name of the list of strings, which contain the sequence to categorize, that will be parsed.
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 6 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 3 : for the results (class, order and superFamily)
		- 1 : for the classification summary (classification_summary)
		- 1 : for the unknown keyword (unknown_keyword)
	@type FINALSUPERFAMILYNAME: string
	@param FINALSUPERFAMILYNAME: name of the superFamily find according the different keywords founded

	@rtype: None
	"""
	####	Regex for serching the name
	name = re.search(r'([^_]+)', FINALSUPERFAMILYNAME).groups()[0]

	# if the name searched is potentialChimeric, we convert the saveType (type of file onto which the sequence will be saved)
	if name == "potentialChimeric":
		####	Pass the saveType from TE in potnetialChimeric
		SEQCLASSIFIED[FEATURES[0]]["saveType"] = "potentialChimeric"
		####	Rewrite the correct saveType of the sequence in the log
		SEQCLASSIFIED[FEATURES[0]]["log"]=SEQCLASSIFIED[FEATURES[0]]["log"].replace("\tTE\t", "\tpotentialChimeric\t")

	####	Define the name of the superFamily sequence as FINALSUPERFAMILYNAME
	SEQCLASSIFIED[FEATURES[0]]["superFamily"] = FINALSUPERFAMILYNAME
