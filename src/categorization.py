#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

import re
import search

def initCategorization(SEQUENCE, SEQCLASSIFIED, BASELINE):
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
	classDetermination(features, SEQCLASSIFIED, BASELINE)


def classDetermination(FEATURES, SEQCLASSIFIED, BASELINE):
	"""
	Determine the class of a sequence.
	For this, complete four dictionnaries, passed onto arguments, that will contain the different catagories that caracterize the sequence.

	Keyword arguments:
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@rtype: None
	"""
	####	Check first if the sequence is chimeric
	if FEATURES[3] != "PotentialChimeric":
		####	Check the class of the sequence if it is not chimeric : class I
		if FEATURES[4] == "I" :
			SEQCLASSIFIED["TE"][FEATURES[0]]={"class":"I"}
			# print(FEATURES[0])
			orderDetermination(FEATURES, SEQCLASSIFIED, BASELINE)
		####	classII
		elif FEATURES[4] == "II":
			SEQCLASSIFIED["TE"][FEATURES[0]]={"class":"II"}
			orderDetermination(FEATURES, SEQCLASSIFIED, BASELINE)
		####	NoCat
		elif FEATURES[4] == "noCat":
			SEQCLASSIFIED["noCat"][FEATURES[0]]={"class":"undefined"}
			# print("NoCat : %s %s" % (FEATURES[4], FEATURES[5]))
			# pass
		####	 NonTE
		else:
			SEQCLASSIFIED["nonTE"][FEATURES[0]]={"class":"nonTE"}
			# print("NA : %s %s" % (FEATURES[4], FEATURES[5]))
	else:
		SEQCLASSIFIED["potentialChimeric"][FEATURES[0]]={"class":"potentialChimeric"}
		# print("chimere : %s %s %s" % (FEATURES[0], FEATURES[4], FEATURES[5]))

def orderDetermination(FEATURES, SEQCLASSIFIED, BASELINE):
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

	@rtype: None
	"""
	####	 The order of the TE is not determined, the superfamily of the TE will not be determined
	if FEATURES[5] == "noCat" :
		SEQCLASSIFIED["TE"][FEATURES[0]]["order"]="unknown"
	####	 The order of the TE is a MITE, a LARD or a TRIM : the superfamily of the TE will not be determined
	elif FEATURES[5] == "MITE" or FEATURES[5] == "LARD" or FEATURES[5] == "TRIM" :
		SEQCLASSIFIED["TE"][FEATURES[0]]["order"]=FEATURES[5]
	####	 The order of the TE is determined : the superfamily of the TE will be determined
	else:
		SEQCLASSIFIED["TE"][FEATURES[0]]["order"]=FEATURES[5]
		superFamilyDetermination(FEATURES, SEQCLASSIFIED, BASELINE)

def superFamilyDetermination(FEATURES, SEQCLASSIFIED, BASELINE):
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

	@todo: Consider if POTENTIALCHIMERIC and NOCAT arguments are mandatory

	@rtype: None
	"""
	try:
		####	 search if there is a 'coding' part in the 7th value of FEATURES => needed to defined the FEATURES superfamily. Search for coding(<anything that is not ));>); NB : regex (?!) anything that is not in ()
		codingRecord = re.search(r'coding=\(((?!\)\);).+?)\);', FEATURES[7]).groups()[0]
		####	Split the coding part to obtain the different results according to the comparison with different databases (3 possibilties: TEBLRtx, TEBLRx and profiles)
		databaseRecords = codingRecord.split(';')
		####	Search if there are different names contained in the codingRecord
		search.searchDifferentName(FEATURES, SEQCLASSIFIED, databaseRecords, BASELINE)

	except AttributeError:
		####	If there is no coding part : declaration of superfamily as undefined
		SEQCLASSIFIED["TE"][FEATURES[0]]["superfamily"] = "undefined"
		return
