#!/usr/bin/python3
# coding: utf8
import re
import comparison

""" @author: Tristan Frances """

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

	@rtype: None.
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
			searchProfilesName(FEATURES, dr, superFamilyFound, BASELINE)
		####	if TE_BLRx or TE_BLRtx is found in coding, search for superFamily name using function searchRepBaseName
		elif dbName=="TE_BLRx" or "TE_BLRtx":
			searchRepBaseName(FEATURES, dr, superFamilyFound)
	####	String that will contain the final supefamily name of the sequence
	finalSuperFamily=""
	####	if multiple names for superfamily are found tor the sequence, proceed to the comparison between all the names found
	if len(superFamilyFound) > 1:
		finalSuperFamily=comparison.superFamilyComparison(FEATURES, superFamilyFound, BASELINE)
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

def searchProfilesName(FEATURES, DATABASERECORD, SUPERFAMILYFOUND, BASELINE):
	"""

	Search the keywords in the profiles part of coding.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type SUPERFAMILYFOUND: list
	@param SUPERFAMILYFOUND: names of the superfamily found for one sequence during the superFamilyDetermination.

	@rtype: TODO
	"""
	#TODO try to differenciate different possible syntax according the first regex found. Ex:
	# case 1: PiRATEdb_CACTA_Tase_NA
	# case 2: _RT_reina_NA_RT_NA
	# case 3: PF05699.9_Dimer_Tnp_hAT_NA_Tase_21.4
	# print(DATABASERECORD)
	####	first split to get rid of profiles:
	DATABASERECORD=DATABASERECORD.split("profiles:")
	keywordFound=False
	####	Parse all the results of profiles and search for regex
	for substr in DATABASERECORD[1].split(','):
		# print(substr+"\n")
		try:
			####	SECOND APPROACH : find keywords by match with BASELINE keyword
			####	Split the substr with _ : the profiles keywords are separated by a
			substr=substr.split("_")
			####	parse the profiles keywords one by one
			for word in substr:
				####	parse the possible name which are in the BASELINE
				for possibleName in BASELINE["nonSpecific"]:
					####	check if the keywords in profiles is in the BASELINE (convert string into lower case berfore comparing)
					if word.lower() == possibleName.lower():
						SUPERFAMILYFOUND.append(word)
						keywordFound=True
				# if not keywordFound:
				# 	print("%s, %s not in BASELINE"%(FEATURES[0], word))

			# print(FEATURES[0], SUPERFAMILYFOUND)
			# print(FEATURES[0], searchObj.groups()[0])
			####	FIRST APPROACH : find keywords by position
			####	Proceed to different operations according the first string of the first regex
			####	If first regex is PiRATEdb or mydatabase (NB maybe should delete this one)
			# searchObj = re.search(r' (_?[^_]+)_([^_]+)_([^_]+)_([^:]+)(?:_[0-9\.]+)?: ([0-9\.]+)%\(([0-9\.]+)%\)', substr)
			# if searchObj.groups()[0] == "PiRATEdb" or searchObj.groups()[0] == "mydatabase":
			# 	# print("Base de Données PiRATE : ", FEATURES[0], searchObj.groups()[1], searchObj.groups()[2])
			# 	pass
			# ####	if second regex begin by PF... (case 2)
			# elif re.search(r'^PF[.\w]+', searchObj.groups()[0]):
			# 	# print(substr)
			# 	# print("Base de Données PF : ", FEATURES[0], searchObj.groups()[-3])
			# 	pass
			# ####	if third regex begin by _... (case 3)
			# elif re.search(r'^_[.\w]+', searchObj.groups()[0]) :
			# 	# print(substr)
			# 	# print("Base de Données _ : ", FEATURES[0], searchObj.groups()[-3])
			# 	pass
			# print(searchObj.groups())
		except AttributeError:
			print('Issue during searching profiles on : '+ substr)


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

	@rtype: None
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
