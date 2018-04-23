#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

import re
import search

def initCategorization(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
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

	@rtype: None
	"""
	features=SEQUENCE.split("\t")
	#TODO Treatment of the superfamily
	classDetermination(features, NONTE, POTENTIALCHIMERIC,  NOCAT, TE, BASELINE)


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

	@rtype: None
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

	@rtype: None
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

	@rtype: None
	"""
	try:
		####	 search if there is a 'coding' part in the 7th value of FEATURES => needed to defined the FEATURES superfamily. Search for coding(<anything that is not ));>); NB : regex (?!) anything that is not in ()
		codingRecord = re.search(r'coding=\(((?!\)\);).+?)\);', FEATURES[7]).groups()[0]
		# codingRecord = re.search(r'coding=\(([^\)]+?)\);?', FEATURES[7]).groups()[0] Another possibility of regex
		####	Split the coding part to obtain the different results according to the comparison with different databases (3 possibilties: TEBLRtx, TEBLRx and profiles)
		databaseRecords = codingRecord.split(';')
		####	Search if there are different names contained in the codingRecord
		search.searchDifferentName(FEATURES, TE, databaseRecords, BASELINE)

		# print(FEATURES[0], dbName, (matches[FEATURES[0]][dbName]))
	except AttributeError:
		####	If there is no coding part : declaration of superfamily as unknown
		TE[FEATURES[0]]["superfamily"] = "unknown"
		# NOCAT[FEATURES[0]]=TE[FEATURES[0]]
		# del TE[FEATURES[0]]
		return
