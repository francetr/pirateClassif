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
		- [3] : potentiel chimeric
		- [4] : class of the sequence
		- [5] : order of the sequence
		- [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@rtype: None
	"""
	features=SEQUENCE.split("\t")
	classDetermination(features, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD)
	SEQCLASSIFIED["log"][features[0]]+="\n"


def classDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD):
	"""
	Determine the class of a sequence.
	For this, complete four dictionnaries, passed onto arguments, that will contain the different catagories that caracterize the sequence.

	Keyword arguments:
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: integer
	@param IDENTITYTHRESHOLD: percentage for which the superFamily name of a sequence will be choosen.

	@rtype: None
	"""
	####	Initialize a log for concerned sequences that will be put into a log file
	SEQCLASSIFIED["log"][FEATURES[0]]=str(FEATURES[0] + "\t" + FEATURES[4] + "\t")
	####	Check first if the sequence is chimeric
	if FEATURES[3] != "PotentialChimeric":
		####	Check the class of the sequence if it is not chimeric : class I
		if FEATURES[4] == "I" :
			SEQCLASSIFIED["TE"][FEATURES[0]]={"class":"I"}
			orderDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD)
		####	classII
		elif FEATURES[4] == "II":
			SEQCLASSIFIED["TE"][FEATURES[0]]={"class":"II"}
			orderDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD)
		####	NoCat
		elif FEATURES[4] == "noCat":
			SEQCLASSIFIED["noCat"][FEATURES[0]]={"class":"undefined"}
		####	 NonTE
		else:
			SEQCLASSIFIED["nonTE"][FEATURES[0]]={"class":"nonTE"}
	else:
		SEQCLASSIFIED["potentialChimeric"][FEATURES[0]]={"class":"potentialChimeric"}

def orderDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD):
	"""
	Determine the order of a sequence. Doing so it complete a dictionnary containing the transposable element.
	If the order wasn't found, the order's sequence is considered as unknown.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: integer
	@param IDENTITYTHRESHOLD: percentage for which the superFamily name of a sequence will be choosen.

	@rtype: None
	"""
	SEQCLASSIFIED["log"][FEATURES[0]] += str("\t" + FEATURES[5] + "\t")
	####	 The order of the TE is not determined, the superfamily of the TE will not be determined
	if FEATURES[5] == "noCat" :
		SEQCLASSIFIED["TE"][FEATURES[0]]["order"]="undefined"
	####	 The order of the TE is a MITE, a LARD or a TRIM : the superfamily of the TE will not be determined
	elif FEATURES[5] == "MITE" or FEATURES[5] == "LARD" or FEATURES[5] == "TRIM" :
		SEQCLASSIFIED["TE"][FEATURES[0]]["order"]=FEATURES[5]
	####	 The order of the TE is determined : the superfamily of the TE will be determined
	else:
		SEQCLASSIFIED["TE"][FEATURES[0]]["order"]=FEATURES[5]
		superFamilyDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD)

def superFamilyDetermination(FEATURES, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD):
	"""
	Determine the superfamily of one sequence. Doing so it complete a dictionnary containing the transposable element.
	If the superFamily can't be defined, it will be unknown.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: integer
	@param IDENTITYTHRESHOLD: percentage for which the superFamily name of a sequence will be choosen.

	@todo: Consider if POTENTIALCHIMERIC and NOCAT arguments are mandatory

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
		SEQCLASSIFIED["TE"][FEATURES[0]]["superfamily"] = "undefined"
		return

def associateSuperFamily(FEATURES, SEQCLASSIFIED, FINALSUPERFAMILYNAME):
	"""
	Associate the superFamily name for the sequence in the corresponding dictionnary (potentialChimeric if find potentialChimeric in FINALSUPERFAMILYNAME; else it will be TE)

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)
	@type FINALSUPERFAMILYNAME: string
	@param FINALSUPERFAMILYNAME: name of the superFamily find according the different keywords founded

	@rtype: None
	"""
	####	Regex for serching the name
	name = re.search(r'([^_]+)', FINALSUPERFAMILYNAME).groups()[0]
	# if the name searched is potentialChimeric, we had the sequence in the potentialChimeric dictionnary and remove it from the TE ditctionnary
	if name == "potentialChimeric":
		####	Add the sequence in the potentialchimeric dictionnary
		SEQCLASSIFIED["potentialChimeric"][FEATURES[0]] = SEQCLASSIFIED["TE"][FEATURES[0]]
		####	Define the name of the superFamily sequence as FINALSUPERFAMILYNAME
		SEQCLASSIFIED["potentialChimeric"][FEATURES[0]]["superFamily"] = FINALSUPERFAMILYNAME
		####	Remove the sequence from TE dictionnary
		del SEQCLASSIFIED["TE"][FEATURES[0]]
	####	Name is not potentialChimeric : can be either undefined or superFamily name (Copia, Gypsy, etc..)
	else:
		####	Define the name of the superFamily sequence as FINALSUPERFAMILYNAME
		SEQCLASSIFIED["TE"][FEATURES[0]]["superFamily"] = FINALSUPERFAMILYNAME
