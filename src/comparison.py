#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

def superFamilyComparison(FEATURES, SUPERFAMILYFOUND, BASELINE):
	"""
	Compare the different superFamilies found for one sequence. For this uses the method combinations of the package itertools allowing to combine 2 by 2 the different elements of a list.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type SUPERFAMILYFOUND: dictionnary
	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily.

	@rtype: string
	@return: name of the superfamily corresponding to the sequence.
	"""
	####	TODO Comparison with a reference base
	####	Draft

	####	First compare differents specific keywords founded for the sequence
	# specificResult = specificComparison(SUPERFAMILYFOUND)

	####	Second compare differents non specific keyword founded for the sequence
	# nonSpecificResult = nonSpecificComparison(SUPERFAMILYFOUND)

	####	Then compare the concordance between the non specific keywords and the specific keywords. Then the result is assigned to a string
	name = compareKeywordsFounded(SUPERFAMILYFOUND["blast"], SUPERFAMILYFOUND["protProfiles"], BASELINE["protProfiles"])

	# print(FEATURES[0], SUPERFAMILYFOUND)
	####	String for the final superFamily name
	return name
	# print("{}, {}, count: {}, len: {}".format(FEATURES[0], SUPERFAMILYFOUND, superFamilyCount, len(SUPERFAMILYFOUND)))

def percentageCalculation(KEYWORDFOUND):
	"""
	Calculus of the percentage of each keyword (can be either specific or non specific, choose as arg) found for a given sequence.

	Keyword arguments:
	@type KEYWORDFOUND: dictionnary
	@param KEYWORDFOUND: All the keywords found (can be either specific or non specific) and their occurrence for a given sequence.

	@rtype: dictionnary
	@return: Percentage foud for each superFamily name in a given sequence.
	"""
	# TODO Check if usefull
	####	Calculus of the percentage of each superFamilyName found for this sequence
	tot = 0
	####	Dictionnary that will contains the percentage of the superFamilyNames found
	percent = {}
	####	Calculus of the total superFamilyName
	for key in KEYWORDFOUND.keys():
		tot+=KEYWORDFOUND[key]

	####	Calculus of the percentage for each superFamilyName found
	for key in KEYWORDFOUND.keys():
		percentage=round(KEYWORDFOUND[key]/tot*100, 1)
		percent[key]=percentage
	return percent

def compareKeywordsFounded(BLAST, PROTPROFILES, BASELINE):
	"""
	Comparison between the specific and the non specific keywords to check if it's in accordance with Wicker classification.
	This classification is implemented in the BASELINE dictionnary. Key = specific keyword : Value = non specific keywords possible

	Keyword arguments:
	@type BLAST: dictionnary
	@param BLAST: All the specific keywords found and their proportion for a given sequence.
	@type PROTPROFILES: dictionnary
	@param PROTPROFILES: All the non specific keywords found and their proportion for a given sequence.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily.

	@rtype: string
	@return: Name of the superFamily with different matches found.
	"""
	name = ""
	####	For this part, 7 cases have been identified
	####	First THREE CASES, Consider NO blast keywords have been founded
	if len(BLAST) == 0:
		####	NO type of proteine profiles keyword have been founded
		if len(PROTPROFILES) == 0:
			name = "unknown"

		####	ONE or MULTIPLE type of proteine profiles keyword have been founded
	elif len(PROTPROFILES) >= 1:
			# name = compareWickerClassification(BLAST, PROTPROFILES, BASELINE)
			pass

	####	Second THREE CASES, Consider ONE blast keywords have been founded
	elif len(BLAST) == 1:
		####	NO type of proteine profiles keyword have been founded
		if len(PROTPROFILES) == 0:
			name = "TE_"
			####	Retrieve the blast superFamily keyword
			for superFamily in BLAST:
				name += str(superFamily)

		####	ONE or MULTIPLE type of proteine profiles keyword have been founded
		elif len(PROTPROFILES) >= 1:
			name = compareWickerClassification(BLAST, PROTPROFILES, BASELINE)
			pass

	####	Last case, MULTIPLE blast keyword have been founded
	elif len(BLAST) > 1:
		name = "potentialChimeric"
		for superFamily in BLAST:
			name += str("_"+superFamily+"_")

	return name

def compareWickerClassification(BLAST, PROTPROFILES, BASELINE):
	"""
	Comparison between the specific and the non specific keywords to check if it's in accordance with Wicker classification.
	Procceed in 2 steps : first check the non spe keyword length
	This classification is implemented in the BASELINE dictionnary. Key = specific keyword : Value = non specific keywords possible

	Keyword arguments:
	@type BLAST: dictionnary
	@param BLAST: All the specific keywords found and their proportion for a given sequence.
	@type PROTPROFILES: dictionnary
	@param PROTPROFILES: All the non specific keywords found and their proportion for a given sequence.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily.

	@rtype: string
	@return: Name of the superFamily with differents matches found.
	"""
	####	Check the coherence between specific and non specific keywords
	name = ""
	keywordsCompared = []
	####	First parse PROTPROFILES keyowrds founded
	for protProfilesFounded in PROTPROFILES:
		keywordsCompared.append(BASELINE[protProfilesFounded])

	print(keywordsCompared)
	####	Second compare the PROTPROFILES keywords founded


	return name


# def specificComparison(SUPERFAMILYFOUND):
# 	"""
# 	Search if there one or multiple specific keyword(s) founded for one sequence.
#
# 	Keyword arguments:
# 	@type SUPERFAMILYFOUND: dictionnary
# 	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
#
# 	@rtype: dictionnary
# 	@return: dictionnary containing the name of superFamily found and its count.
# 	{B{Key} = specific : B{I{Values}} = count of this keyword}\n
# 	3 return possible.
#
# 		1. If no keyword => None: 0
#
# 		2. If 1 keyword => name of keyword : count of this keyword
#
# 		3. If more than 1 keyword => keywords founded : proportion of the keyword
# 	"""
# 	resultSpecific={}
#
# 	####	If NO specific keyword has been found
# 	if len(SUPERFAMILYFOUND["specific"]) == 0:
# 		resultSpecific = {"undefined": 0}
#
# 	####	If ONE same specific keyword has been found (ex: Mariner, Mariner)
# 	elif len(SUPERFAMILYFOUND["specific"]) == 1:
# 		for superFamilyName in SUPERFAMILYFOUND["specific"].keys():
# 			resultSpecific = {superFamilyName: SUPERFAMILYFOUND["specific"][superFamilyName]}
#
# 	####	If MULTIPLE different specific keywords have been found, considered as potential chimeric (ex: Mariner, Copia)
# 	elif len(SUPERFAMILYFOUND["specific"]) >= 1:
# 		####	Convert the counter of keywords founded into a proportion
# 		percentName = percentageCalculation(SUPERFAMILYFOUND["specific"])
# 		resultSpecific = percentName
#
# 	return resultSpecific
#
# def nonSpecificComparison(SUPERFAMILYFOUND):
# 	"""
# 	Search if there one or multiple non specific keyword(s) founded for one sequence.
#
# 	Keyword arguments:
# 	@type SUPERFAMILYFOUND: dictionnary
# 	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
#
# 	@rtype: tuple
# 	@return: name of the superfamily corresponding to the sequence and number of non specific keyword(s) founded for this sequence.
# 	3 return possible.
#
# 		1. If no keyword => None: 0
#
# 		2. If 1 keyword => name of keyword : count of this keyword
#
# 		3. If more than 1 keyword => keywords founded : proportion of the keyword
#
# 	"""
# 	resultNonSpecific={}
# 	####	If NO non specific keyword has been found
# 	if len(SUPERFAMILYFOUND["nonSpecific"]) == 0:
# 		resultNonSpecific = {"undefined": 0}
#
# 	####	If ONE same non specific keyword has been found (ex: Mariner, Mariner)
# 	elif len(SUPERFAMILYFOUND["nonSpecific"]) == 1:
# 		for superFamilyName in SUPERFAMILYFOUND["nonSpecific"].keys():
# 			resultNonSpecific= {superFamilyName: SUPERFAMILYFOUND["nonSpecific"][superFamilyName]}
#
# 	####	If MULTIPLE different non specific keywords have been found, considered as potential chimeric (ex: Mariner, Copia)
# 	elif len(SUPERFAMILYFOUND["nonSpecific"]) >= 1:
# 			####	Convert the counter of non specific keywords founded into a proportion
# 			percentName = percentageCalculation(SUPERFAMILYFOUND["nonSpecific"])
# 			resultNonSpecific= percentName
#
# 	return resultNonSpecific
