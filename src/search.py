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
			searchRepBaseName(dr, superFamilyFound, BASELINE)
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
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type SUPERFAMILYFOUND: dictionnary
	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@rtype: TODO
	"""
	#TODO try to differenciate different possible syntax according the first regex found. Ex:
	# case 1: PiRATEdb_CACTA_Tase_NA
	# case 2: _RT_reina_NA_RT_NA
	# case 3: PF05699.9_Dimer_Tnp_hAT_NA_Tase_21.4
	# print(DATABASERECORD)
	####	first split to get rid of profiles:
	DATABASERECORD=DATABASERECORD.split("profiles:")
	####	Parse all the results of profiles and search for regex
	for substr in DATABASERECORD[1].split(','):
		# print(substr+"\n")
		try:
			####	SECOND APPROACH : find keywords by match with BASELINE keyword
			####	Counter for the keyword found
			cptKeyword=0

			####	parse the possible name which are in the BASELINE
			for nonSpecificName in BASELINE["nonSpecific"]:
				####	Search in the BASELINE file if there is a keyword matching with the profiles
				keywordSearch = re.search(r'_%s_'%(re.escape(nonSpecificName)), substr)
				####	Check if the keywords in profiles is in the BASELINE (convert string into lower case berfore comparing)
				if keywordSearch:
					print(keywordSearch)
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
					# ####	If the keyword hasn't been found in the SUPERFAMILYFOUND dictionnary, its counter is 1
					# if not nonSpecificName in SUPERFAMILYFOUND["nonSpecific"]:
					# 	cptKeyword+=1
					# 	if nonSpecificName == "Tase" or nonSpecificName == "Tase*":
					# 		SUPERFAMILYFOUND["nonSpecific"]["Tase"]=cptKeyword
					# 	else:
					# 		SUPERFAMILYFOUND["nonSpecific"][nonSpecificName]=cptKeyword
					# ####	If the keyword has already been found in the SUPERFAMILYFOUND dictionnary, its counter is incremented by 1
					# else:
					# 	if nonSpecificName == "Tase" or nonSpecificName == "Tase*":
					# 		SUPERFAMILYFOUND["nonSpecific"]["Tase"]+=1
					# 	else:
					# 		SUPERFAMILYFOUND["nonSpecific"][nonSpecificName]+=1

						# keywordFound=True

				# if not keywordFound:
				# 	print("%s, %s not in BASELINE"%(FEATURES[0], word))

			# print(FEATURES[0], SUPERFAMILYFOUND)
			# print(FEATURES[0], keywordSearch.groups()[0])
			####	FIRST APPROACH : find keywords by position
			####	Proceed to different operations according the first string of the first regex
			####	If first regex is PiRATEdb or mydatabase (NB maybe should delete this one)
			# keywordSearch = re.search(r' (_?[^_]+)_([^_]+)_([^_]+)_([^:]+)(?:_[0-9\.]+)?: ([0-9\.]+)%\(([0-9\.]+)%\)', substr)
			# if keywordSearch.groups()[0] == "PiRATEdb" or keywordSearch.groups()[0] == "mydatabase":
			# 	# print("Base de Données PiRATE : ", FEATURES[0], keywordSearch.groups()[1], keywordSearch.groups()[2])
			# 	pass
			# ####	if second regex begin by PF... (case 2)
			# elif re.search(r'^PF[.\w]+', keywordSearch.groups()[0]):
			# 	# print(substr)
			# 	# print("Base de Données PF : ", FEATURES[0], keywordSearch.groups()[-3])
			# 	pass
			# ####	if third regex begin by _... (case 3)
			# elif re.search(r'^_[.\w]+', keywordSearch.groups()[0]) :
			# 	# print(substr)
			# 	# print("Base de Données _ : ", FEATURES[0], keywordSearch.groups()[-3])
			# 	pass
			# print(keywordSearch.groups())
		except AttributeError:
			print('Issue during searching profiles on : '+ substr)


def searchRepBaseName(DATABASERECORD, SUPERFAMILYFOUND, BASELINE):
	"""

	Search the keywords in the TE_BLRx and TE_BLRtx part of coding. Then check if the keyword is founded into the BASELINE file.

	Keyword arguments:
	@type DATABASERECORD: string
	@param DATABASERECORD: profiles that will be parsed to find the keywords.
	@type SUPERFAMILYFOUND: dictionnary
	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).

	@rtype: None
	"""
	cptKeyword=0
	for substr in DATABASERECORD.split(','):
		try:
			####	regex split into groups. Group's number correspondance : 0 : name of the sequence; 1 : type of class; 2 : type of order; 3 : superFamily name
			keywordSearch = re.search(r' ([^:]+):(?:Class)?(I+|\?):([^:]+):([^:]+): ([0-9\.]+)%', substr)
			keyword = keywordSearch.groups()[3]
			####	Counter for the keyword found


			####	Parse the key of the specific keywords in the BASELINE file
			for specificName in BASELINE["specific"]:
				####	Check if the keywords in profiles is in the BASELINE
				if keyword in BASELINE["specific"][specificName]:
					####	If the keyword hasn't been found in the SUPERFAMILYFOUND dictionnary, its counter is 1
					if not keyword in SUPERFAMILYFOUND["specific"]:
						cptKeyword+=1
						SUPERFAMILYFOUND["specific"][specificName]=cptKeyword
					####	If the keyword has already been found in the SUPERFAMILYFOUND dictionnary, its counter is incremented by 1
					else:
						SUPERFAMILYFOUND["specific"][specificName]+=1

		except AttributeError:
			print('Issue during searching RepBase name on : '+ substr)
