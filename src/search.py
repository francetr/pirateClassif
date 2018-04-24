#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

import re
import comparison

def searchDifferentName(FEATURES, TE, DATABASERECORDS, BASELINE):
	"""
	If the superFamily can't be defined, it will be unknown.
	The name is recovered from two : the hmm profiles (profiles part) and REPET (TE_BLRx and TE_BLRtx part).

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence. Usefull to know the sequence concerned.
	@type DATABASERECORDS: list
	@param DATABASERECORDS: list of string in which there are the superFamily name to search.
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@rtype: None.
	"""
	####	Dictionnary that will contains (or not) the different superFamilies find for the concerned sequence
	# matches={FEATURES[0]:{}}
	####	Dictionnary with key = specific : value = { key = specifc keyword found : value = number of presence of this keyword},
	####					 key = nonSpecific : value = { key = non specifc keyword found : value = number of presence of this keyword}}
	superFamilyFound={"specific":{}, "nonSpecific":{}}

	####	Scan the different results obtain for each comparisons
	for dr in DATABASERECORDS:
		####	Split to obtain the name of the database used (TE_BLRx or TE_BLRtx)
		dbName = dr.split(':')[0].strip()
		# print(FEATURES[0], dbName)
		# matches[FEATURES[0]] = {dbName:[]}
		####	if profiles is found in coding, search for superFamily name using function searchProfilesName
		if dbName=='profiles':
			searchProfilesName(FEATURES, dr, superFamilyFound, BASELINE)
		####	if TE_BLRx or TE_BLRtx is found in coding, search for superFamily name using function searchRepBaseName
		elif dbName=="TE_BLRx" or "TE_BLRtx":
			searchRepBaseName(FEATURES, dr, superFamilyFound, BASELINE)
	print(FEATURES[0], superFamilyFound)
	####	String that will contain the final supefamily name of the sequence
	finalSuperFamily = ""
	####	if no keywords have been found, finalSuperFamily will be unknown
	if len(superFamilyFound["specific"]) == 0 and len(superFamilyFound["nonSpecific"]) == 0 :
		finalSuperFamily = "unknown"
	####	else, a comparison between the keywords founded is done
	else:
		comparison.superFamilyComparison(FEATURES, superFamilyFound, BASELINE)
		# finalSuperFamily=comparison.superFamilyComparison(FEATURES, superFamilyFound, BASELINE)
	# ####	One superFamily name have been found
	# elif len(superFamilyFound) == 1:
	# 	####	Replace the superfamily name by "unknown"
	# 	if superFamilyFound[0] == "?":
	# 		finalSuperFamily = "unknown"
	# 	else:
	# 		finalSuperFamily=superFamilyFound[0]
	# ####	No superFamily name have been found
	# else:
	# 	finalSuperFamily="unknown"
	# # matches[FEATURES[0]][dbName].append(keywordSearch.groups())
	# TE[FEATURES[0]]["superFamily"]=finalSuperFamily

def searchProfilesName(FEATURES, DATABASERECORD, SUPERFAMILYFOUND, BASELINE):
	"""

	Search the keywords in the profiles part of coding.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence. Usefull to know the sequence concerned. Usefull to know the sequence concerned.
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type SUPERFAMILYFOUND: dictionnary
	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@rtype: TODO
	"""
	# print(DATABASERECORD)
	####	first split to get rid of profiles:
	DATABASERECORD=DATABASERECORD.split("profiles:")
	####	Parse all the results of profiles and search for regex
	for substr in DATABASERECORD[1].split(','):
		# print(substr+"\n")
		try:
			####	APPROACH : Find keywords by match with BASELINE keyword
			####	Counter for the keyword found
			cptKeyword=0
			####	Boolean allowing to know if the keyword is in the baseline
			keywordFound=False
			####	parse the possible name which are in the BASELINE
			for nonSpecificName in BASELINE["nonSpecific"]:
				####	Search in the BASELINE file if there is a keyword matching with the profiles. NB: need to escape the keyword because of special character like *
				keywordSearch = re.search(r'_%s_'%(re.escape(nonSpecificName)), substr)
				####	Check if the keywords in profiles is in the BASELINE (convert string into lower case berfore comparing)
				if keywordSearch:
					keywordFound=True
					####	If the keyword hasn't been found in the SUPERFAMILYFOUND dictionnary, its counter is 1
					if not nonSpecificName in SUPERFAMILYFOUND["nonSpecific"]:
						cptKeyword+=1
						if nonSpecificName == "Tase" or nonSpecificName == "Tase*":
							SUPERFAMILYFOUND["nonSpecific"]["Tase"]=cptKeyword
						else:
							SUPERFAMILYFOUND["nonSpecific"][nonSpecificName]=cptKeyword
					####	If the keyword has already been found in the SUPERFAMILYFOUND dictionnary, its counter is incremented by 1
					else:
						if nonSpecificName == "Tase" or nonSpecificName == "Tase*":
							SUPERFAMILYFOUND["nonSpecific"]["Tase"]+=1
						else:
							SUPERFAMILYFOUND["nonSpecific"][nonSpecificName]+=1

			####	If there is no matches between the string and the baseline, print the string
			if not keywordFound:
				print("/!\	Sequence : %s No match in the string %s with BASELINE"%(FEATURES[0], substr))

			# print(FEATURES[0], SUPERFAMILYFOUND)
		except AttributeError:
			print('Issue during searching profiles on : '+ substr)


def searchRepBaseName(FEATURES, DATABASERECORD, SUPERFAMILYFOUND, BASELINE):
	"""

	Search the keywords in the TE_BLRx and TE_BLRtx part of coding. Then check if the keyword is founded into the BASELINE file.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence. Usefull to know the sequence concerned. Usefull to know the sequence concerned.
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type SUPERFAMILYFOUND: dictionnary
	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@rtype: None
	"""
	####	Counter for the keyword found
	cptKeyword=0
	for substr in DATABASERECORD.split(','):
		try:
			####	regex split into groups. Group's number correspondance : 0 : name of the sequence; 1 : type of class; 2 : type of order; 3 : superFamily name
			keywordSearch = re.search(r' ([^:]+):(?:Class)?(I+|\?):([^:]+):([^:]+): ([0-9\.]+)%', substr)
			keyword = keywordSearch.groups()[3]

			####	Boolean allowing to know if the keyword is in the baseline
			keywordFound=False

			####	Parse the key of the specific keywords in the BASELINE file
			for specificName in BASELINE["specific"]:
				####	Check if the keywords in profiles is in the BASELINE
				if keyword in BASELINE["specific"][specificName]:
					keywordFound=True
					####	If the keyword hasn't been found in the SUPERFAMILYFOUND dictionnary, its counter is 1
					if not keyword in SUPERFAMILYFOUND["specific"]:
						cptKeyword+=1
						SUPERFAMILYFOUND["specific"][specificName]=cptKeyword
					####	If the keyword has already been found in the SUPERFAMILYFOUND dictionnary, its counter is incremented by 1
					else:
						SUPERFAMILYFOUND["specific"][specificName]+=1

			####	If there is no matches between the string and the baseline, print the string
			if not keywordFound:
				print("/!\	Sequence : %s No match in the string %s with BASELINE"%(FEATURES[0], substr))

		except AttributeError:
			print('Issue during searching RepBase name on : '+ substr)
