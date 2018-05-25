#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

import re
import comparison

def searchDifferentName(FEATURES, SEQCLASSIFIED, DATABASERECORDS, BASELINE, IDENTITYTHRESHOLD):
	"""
	If the superFamily can't be defined, it will be unknown.
	The name is recovered from two : the hmm profiles (profiles part) and REPET (TE_BLRx and TE_BLRtx part).

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence. Usefull to know the sequence concerned.
	@type DATABASERECORDS: list
	@param DATABASERECORDS: list of string in which there are the superFamily name to search.
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: integer
	@param IDENTITYTHRESHOLD: percentage for which the superFamily name of a sequence will be choosen.

	@rtype: string
	@return: SuperFamily name of the sequence
	"""
	####	Dictionnary that will contains (or not) the different superFamilies find for the concerned sequence
	####		key = blast : value = { key = blast keyword found : value = number of presence of this keyword},
	####		key = protProfiles : value = { key = proteine Profiles keyword found : value = number of presence of this keyword}}
	superFamilyFound={"blast":{}, "protProfiles":{}}

	####	Scan the different results obtain for each comparisons
	for dr in DATABASERECORDS:
		####	Split to obtain the name of the database used (TE_BLRx or TE_BLRtx)
		dbName = dr.split(':')[0].strip()
		####	if profiles is found in coding, search for superFamily name using function searchProfilesName
		if dbName=='profiles':
			searchProfilesName(FEATURES, dr, superFamilyFound["protProfiles"], BASELINE)
		####	if TE_BLRx or TE_BLRtx is found in coding, search for superFamily name using function searchRepBaseName
		elif dbName=="TE_BLRx" or "TE_BLRtx":
			searchBlastName(FEATURES, dr, superFamilyFound["blast"], BASELINE)
	####	Convert the count of profiles keywords founded into percentage
	percentageCalculation(superFamilyFound["protProfiles"])
	####	Convert the count of blast name founded into percentage
	percentageCalculation(superFamilyFound["blast"])
	####	String that will contain the final supefamily name of the sequence
	finalSuperFamilyName = ""
	####	Do a comparison between the keywords founded
	finalSuperFamilyName = comparison.superFamilyComparison(superFamilyFound, BASELINE, IDENTITYTHRESHOLD)
	print("id: %s; BLAST : %s; PROTPROFILES : %s; FINALNAME : %s "%(FEATURES[0], superFamilyFound["blast"], superFamilyFound["protProfiles"], finalSuperFamilyName))
	return finalSuperFamilyName

def searchProfilesName(FEATURES, DATABASERECORD, PROFILESFOUND, BASELINE):
	"""

	Search the keywords in the profiles part of coding.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence. Usefull to know the sequence concerned. Usefull to know the sequence concerned.
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type PROFILESFOUND: dictionnary
	@param PROFILESFOUND: count of the keywords of the proteines profiles found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different profiles keyword possible for a given superfamily (usefull for the function superFamilyComparison).

	@rtype: None
	"""
	####	first split to get rid of profiles:
	DATABASERECORD=DATABASERECORD.split("profiles:")
	####	Parse all the results of profiles and search for regex
	for substr in DATABASERECORD[1].split(','):
		try:
			####	APPROACH : Find keywords by match with BASELINE keyword
			####	Boolean allowing to know if the keyword is in the baseline
			keywordFound=False
			####	parse the possible name which are in the BASELINE
			for protProfilesKeyword in BASELINE["protProfiles"]:
				####	Search in the BASELINE file if there is a keyword matching with the profiles. NB: need to escape the keyword because of special character like *
				keywordSearch = re.search(r'_%s_'%(re.escape(protProfilesKeyword)), substr, re.IGNORECASE)
				####	Check if the keywords in profiles is in the BASELINE (convert string into lower case berfore comparing)
				if keywordSearch:
					keywordFound=True
					####	If the keyword hasn't been found in the PROFILESFOUND dictionnary, its counter is 1
					if not protProfilesKeyword in PROFILESFOUND:
						####	Consider that keyword Tase and Tase* have the same counter
						if protProfilesKeyword == "Tase" or protProfilesKeyword == "Tase*":
							PROFILESFOUND["Tase"]=1
						####	Counter for all the rest of the possibles keywords
						else:
							PROFILESFOUND[protProfilesKeyword]=1

					####	If the keyword has already been found in the PROFILESFOUND dictionnary, its counter is incremented by 1
					else:
						if protProfilesKeyword == "Tase" or protProfilesKeyword == "Tase*":
							PROFILESFOUND["Tase"]+=1
						else:
							PROFILESFOUND[protProfilesKeyword]+=1

			####	If there is no matches between the string and the baseline, print the string
			if not keywordFound:
				print("/!\	Sequence : %s No match in the string %s with BASELINE"%(FEATURES[0], substr))

		except AttributeError:
			print('Issue during searching profiles on : '+ substr)


def searchBlastName(FEATURES, DATABASERECORD, BLASTFOUND, BASELINE):
	"""

	Search the keywords in the TE_BLRx and TE_BLRtx part of coding. Then check if the keyword is founded into the BASELINE file.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence. Usefull to know the sequence concerned. Usefull to know the sequence concerned.
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type BLASTFOUND: dictionnary
	@param BLASTFOUND: count of the blast name found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@rtype: None
	"""
	for substr in DATABASERECORD.split(','):
		try:
			####	Approach: use regex to find the keyword
			####	regex split into groups. Group's number correspondance : 0 : name of the sequence; 1 : type of class; 2 : type of order; 3 : superFamily name
			keywordSearch = re.search(r' ?([^:]+):(?:Class)?(I+|\?):([^:]+):([^:]+): ([0-9\.]+)%', substr)
			keyword = keywordSearch.groups()[3]
			####	Boolean allowing to know if the keyword is in the baseline
			keywordFound=False
			####	Parse the key of the specific keywords in the BASELINE file
			for blastKeyword in BASELINE["blast"]:
				####	Check if the keywords in profiles is in the BASELINE
				if keyword in BASELINE["blast"][blastKeyword]:
					keywordFound=True
					####	If the keyword hasn't been found in the BLASTFOUND dictionnary, we creates its value in BLASTFOUND and define its counter as 1
					if not blastKeyword in BLASTFOUND:
						BLASTFOUND[blastKeyword]=1
					####	If the keyword has already been found in the BLASTFOUND dictionnary, its counter is incremented by 1
					else:
						BLASTFOUND[blastKeyword]+=1

			####	If there is no matches between the string and the baseline, print the string
			if not keywordFound:
				print("/!\	Sequence : %s No match in the string %s with BASELINE"%(FEATURES[0], substr))

		except AttributeError:
			print('Issue during searching RepBase name on : '+ substr)

def percentageCalculation(KEYWORDFOUND):
	"""
	Convert the count of each keyword of KEYWORDFOUND (can be either blast or protProfiles) into a percentage.

	Keyword arguments:
	@type KEYWORDFOUND: dictionnary
	@param KEYWORDFOUND: All the keywords found (can be either blast or proteines) and their occurrence for a given sequence.

	@rtype: None
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
		KEYWORDFOUND[key]=percentage
