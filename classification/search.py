#!/usr/bin/python2
# coding: utf8

""" @author: Tristan Frances """

import re
import comparison
import math

def searchDifferentName(FEATURES, SEQCLASSIFIED, DATABASERECORDS, BASELINE, IDENTITYTHRESHOLD):
	"""
	If the superFamily can't be defined, it will be unknown.
	The name is recovered from two : the hmm profiles (profiles part) and REPET (TE_BLRx and TE_BLRtx part).

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type DATABASERECORDS: list
	@param DATABASERECORDS: list of string in which there are the superFamily name to search.
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 8 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 4 : for the results (length, completeness, class and finalDegree), that we'll find for each sequences
		- 3 : for the order, for the predictedSuperFamily and 1 for the proofs. (These 2 dic are just for TE or potentialChimeric)
		- 1 : for the unknown keyword (unknown_keyword)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: float
	@param IDENTITYTHRESHOLD: threshold for which the superFamily of a sequence will be named, if its frequency in the sequence is greater or equal to this threshold.

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
			searchProfilesName(FEATURES, SEQCLASSIFIED, dr, superFamilyFound["protProfiles"], BASELINE)
		####	if TE_BLRx or TE_BLRtx is found in coding, search for superFamily name using function searchRepBaseName
		elif dbName=="TE_BLRx" or "TE_BLRtx":
			searchBlastName(FEATURES, SEQCLASSIFIED, dr, superFamilyFound["blast"], BASELINE)
	####	Convert the count of profiles keywords founded into percentage
	percentageCalculation(superFamilyFound["protProfiles"])
	####	Convert the count of blast name founded into percentage
	percentageCalculation(superFamilyFound["blast"])
	####	String that will contain the final supefamily name of the sequence
	finalSuperFamilyName = ""
	####	Do a comparison between the keywords founded
	finalSuperFamilyName = comparison.superFamilyComparison(superFamilyFound, BASELINE, IDENTITYTHRESHOLD)
	####	 save the different proofs found for the name determination of the sequence
	SEQCLASSIFIED[FEATURES[0]]["superFamilyProofs"] = superFamilyFound
	return finalSuperFamilyName

def searchProfilesName(FEATURES, SEQCLASSIFIED, DATABASERECORD, PROFILESFOUND, BASELINE):
	"""
	Search the keywords in the profiles part of coding.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@param SEQCLASSIFIED: Used in case a keyword hasn't be found. dictionnary storing the result of the classification into 6 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 3 : for the results (class, order and superFamily)
		- 1 : for the classification summary (classification_summary)
		- 1 : for the unknown keyword (unknown_keyword)
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
							PROFILESFOUND["Tase"]=[1]
						####	Counter for all the rest of the possibles keywords
						else:
							PROFILESFOUND[protProfilesKeyword]=[1]

					####	If the keyword has already been found in the PROFILESFOUND dictionnary, its counter is incremented by 1
					else:
						if protProfilesKeyword == "Tase" or protProfilesKeyword == "Tase*":
							PROFILESFOUND["Tase"][0]+=1
						else:
							PROFILESFOUND[protProfilesKeyword][0]+=1

			####	If there is no matches between the string and the baseline, print the string
			if not keywordFound:
				print("/!\	Sequence : %s No match in the string %s with BASELINE"%(FEATURES[0], substr))
				SEQCLASSIFIED[FEATURES[0]]["unknown_keyword"] += str("/!\ Proteine profile search for sequence : %s, no match in the string %s with BASELINE. Counting it as undefined\n"%(FEATURES[0], substr))
				####	Add an undefined keywords for PROFILESFOUND
				if not "undefined" in PROFILESFOUND:
					PROFILESFOUND["undefined"]=[1]
				else:
					PROFILESFOUND["undefined"][0]+=1

		except AttributeError:
			print('Issue during searching profiles on : '+ substr)
			SEQCLASSIFIED[FEATURES[0]]["unknown_keyword"] += str("/!\ Proteine profiles search for sequence : %s, issue with the string %s with BASELINE\n"%(FEATURES[0], substr))
			if not "undefined" in PROFILESFOUND:
				PROFILESFOUND["undefined"]=[1]
			else:
				PROFILESFOUND["undefined"][0]+=1

def searchBlastName(FEATURES, SEQCLASSIFIED, DATABASERECORD, BLASTFOUND, BASELINE):
	"""
	Search the keywords in the TE_BLRx and TE_BLRtx part of coding. Then check if the keyword is founded into the BASELINE file.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 8 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 4 : for the results (length, completeness, class and finalDegree), that we'll find for each sequences
		- 3 : for the order, for the predictedSuperFamily and 1 for the proofs. (These 2 dic are just for TE or potentialChimeric)
		- 1 : for the unknown keyword (unknown_keyword)	@type DATABASERECORD: string
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
			keywordFound=False
			####	Boolean allowing to know if the keyword is in the baseline
			####	Parse the key of the specific keywords in the BASELINE file
			for blastKeyword in BASELINE["blast"]:
				####	Check if the keywords in profiles is in the BASELINE
				if keyword in BASELINE["blast"][blastKeyword]:
					keywordFound=True
					####	If the keyword hasn't been found in the BLASTFOUND dictionnary, we creates its value in BLASTFOUND and define its counter as 1
					if not blastKeyword in BLASTFOUND:
						BLASTFOUND[blastKeyword]=[1]
					####	If the keyword has already been found in the BLASTFOUND dictionnary, its counter is incremented by 1
					else:
						BLASTFOUND[blastKeyword][0]+=1

			####	If there is no matches between the string and the baseline, print the string
			if not keywordFound:
				print("/!\	Sequence : %s No match in the string %s with BASELINE"%(FEATURES[0], substr))
				SEQCLASSIFIED[FEATURES[0]]["unknown_keyword"] += str("/!\ BLAST search for sequence : %s, no match in the string %s with BASELINE\n"%(FEATURES[0], substr))
				####	Add an undefined keywords for BLASTFOUND
				if not "undefined" in BLASTFOUND:
					BLASTFOUND["undefined"]=[1]
				else:
					BLASTFOUND["undefined"][0]+=1

		except AttributeError:
			print('Issue during BLAST search name on : '+ substr)
			SEQCLASSIFIED[FEATURES[0]]["unknown_keyword"] += str("/!\ BLAST search for sequence : %s, issue for string %s with BASELINE\n"%(FEATURES[0], substr))
			if not "undefined" in BLASTFOUND:
				BLASTFOUND["undefined"]=[1]
			else:
				BLASTFOUND["undefined"][0]+=1

def percentageCalculation(KEYWORDFOUND):
	"""
	Convert the count of each keyword of KEYWORDFOUND (can be either blast or protProfiles) into a percentage.

	Keyword arguments:
	@type KEYWORDFOUND: dictionnary
	@param KEYWORDFOUND: All the keywords found (can be either blast or proteines) and their occurrence for a given sequence.

	@rtype: None
	"""
	####	Calculus of the percentage of each superFamilyName found for this sequence
	tot = float(0.0)

	####	Calculus of the total keywords found in KEYWORDFOUND (can be BLAST or PROTPROFILES)
	for key in KEYWORDFOUND.keys():
		tot+=KEYWORDFOUND[key][0]

	####	Calculus of the percentage (round to 2 digits after the comma) for each superFamilyName found
	for key in KEYWORDFOUND.keys():
		percentage=round(float(KEYWORDFOUND[key][0]/tot*100), 2)
		####	Add the value of percentage in a list where index [0] = counter of keyword and index [1] = percentage of this keyword
		####	/!\ It is possible that the percentage is != 100% because we round up to the last digit
		KEYWORDFOUND[key].append(percentage)
