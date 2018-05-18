#!/usr/bin/python3
# coding: utf8
from collections import Counter

""" @author: Tristan Frances """

def superFamilyComparison(FEATURES, SUPERFAMILYFOUND, SUPERFAMILYASSOCIATED):
	"""
	Compare the different superFamilies found for one sequence. For this uses the method combinations of the package itertools allowing to combine 2 by 2 the different elements of a list.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type SUPERFAMILYFOUND: dictionnary
	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
	@type SUPERFAMILYASSOCIATED: dictionnary
	@param SUPERFAMILYASSOCIATED: dictionnary containing different superfamily names possible for a given superfamily.

	@rtype: string
	@return: name of the superfamily corresponding to the sequence.
	"""
	####	Compare all the keywords founded for one sequence (providing from blast or protProfiles) . Then the result is assigned to a string
	name = compareKeywordsFounded(SUPERFAMILYFOUND["blast"], SUPERFAMILYFOUND["protProfiles"], SUPERFAMILYASSOCIATED["protProfiles"])
	####	String for the final superFamily name
	return name

def compareKeywordsFounded(BLAST, PROTPROFILES, SUPERFAMILYASSOCIATED):
	"""
	Comparison between the blast and the proteines profiles keywords to check if it's in accordance with Wicker classification.
	This classification is implemented in the SUPERFAMILYASSOCIATED dictionnary. Key = blast keyword : Value = protProfile keywords possible

	Keyword arguments:
	@type BLAST: dictionnary
	@param BLAST: All the blast keywords found and their proportion for a given sequence.
	@type PROTPROFILES: dictionnary
	@param PROTPROFILES: All the proteines profiles keywords found and their proportion for a given sequence.
	@type SUPERFAMILYASSOCIATED: dictionnary
	@param SUPERFAMILYASSOCIATED: dictionnary containing different superfamily names possible for a given proteine profile.

	@rtype: string
	@return: Name of the superFamily with different matches found.
	"""
	name = ""
	####	For this part, 7 cases have been identified
	####	First THREE CASES, Consider NO BLAST keywords have been founded
	if len(BLAST) == 0:
		####	NO type of PROPROFILES keyword have been founded
		if len(PROTPROFILES) == 0:
			name = "unknown"

		####	ONE type of PROPROFILES keyword have been founded
		elif len(PROTPROFILES) == 1:
			name = compareSingleProtProfiles(PROTPROFILES, SUPERFAMILYASSOCIATED)

		####	MULTIPLE type of PROPROFILES keyword have been founded
		elif len(PROTPROFILES) > 1:
			name = compareMultipleProtProfiles(PROTPROFILES, SUPERFAMILYASSOCIATED)

	####	Second THREE CASES, Consider ONE blast keywords have been founded
	elif len(BLAST) == 1:
		####	NO type of PROPROFILES keyword have been founded
		if len(PROTPROFILES) == 0:
			####	Retrieve the BLAST superFamily keyword
			for superFamily in BLAST:
				name += str(superFamily)

		####	ONE or MULTIPLE type of PROPROFILES keyword have been founded
		elif len(PROTPROFILES) >= 1:
			name = compareWickerClassification(BLAST, PROTPROFILES, SUPERFAMILYASSOCIATED)

	####	Last case, MULTIPLE BLAST keywords have been founded
	elif len(BLAST) > 1:
		name = "potentialChimeric"
		for superFamily in BLAST:
			name += str("_"+superFamily)

	return name

def compareSingleProtProfiles(PROTPROFILE, SUPERFAMILYASSOCIATED):
	"""
	Comparison between the blast and the proteines profiles keywords to check if it's in accordance with Wicker classification.
	This classification is implemented in the SUPERFAMILYASSOCIATED dictionnary. Key = blast keyword : Value = protProfile keywords possible

	Keyword arguments:
	@type PROTPROFILE: dictionnary
	@param PROTPROFILE: Proteine profile keyword found and its proportion for a given sequence.
	@type SUPERFAMILYASSOCIATED: dictionnary
	@param SUPERFAMILYASSOCIATED: dictionnary containing different superfamily names possible for a given proteine profile.

	@rtype: string
	@return: Name of the superFamily.
	"""
	####	APPROACH : Add the superFamily names associated to the PROTPROFILE. Then Check if the superFamily name associated is uniq or not
	####	String that will contain the superFamily name
	name = ""
	####	List that will contain the superFamily name associated to PROTPROFILE
	superFamilyAssociated = retrieveSuperFamilyAssociated(PROTPROFILE, SUPERFAMILYASSOCIATED)

	####	If only ONE superFamily name is associated to PROTPROFILE
	if len(superFamilyAssociated) == 1:
		name = superFamilyAssociated[0]
	####	Else MULTIPLE superFamily names are associated to PROTPROFILE, name is undefined with the PROTPROFILE
	else:
		name="undefined"
		for proteine in PROTPROFILE:
			name += str("_" + proteine)

	return name

def compareMultipleProtProfiles(PROTPROFILES, SUPERFAMILYASSOCIATED):
	"""
	Comparison between the blast and the proteines profiles keywords to check if it's in accordance with Wicker classification.
	This classification is implemented in the SUPERFAMILYASSOCIATED dictionnary. Key = blast keyword : Value = protProfile keywords possible

	Keyword arguments:
	@type PROTPROFILES: dictionnary
	@param PROTPROFILES: All the proteines profiles keywords found and their proportion for a given sequence.
	@type SUPERFAMILYASSOCIATED: dictionnary
	@param SUPERFAMILYASSOCIATED: dictionnary containing different superfamily names possible for a given proteine profile.

	@rtype: string
	@return: Name of the superFamily.
	"""
	####	APPROACH : Proceed in 2 steps,
	####	First creates a list were all the PROTPROFILES founded that are not uniq, are added (allow to check if the proteines profiles are associated to a common superFamily)
	####	Then check if the superFamily associated founded is uniq. Finally: Assigns this name as final name else, as undefined with the matches.

	####	String that will contain the superFamily name
	name = ""
	####	List that will contain all the superFamily name associated to PROTPROFILES
	superFamilyAssociated = retrieveSuperFamilyAssociated(PROTPROFILES, SUPERFAMILYASSOCIATED)

	####	Count all the superFamily names associated
	superFamilyCount = Counter(superFamilyAssociated)
	####	Remove all the superFamily names founded just once
	uniqSuperFamilyAssociated = list(set(k for k in superFamilyAssociated if superFamilyCount[k]>1))

	####	If only ONE superFamily name is associated to PROTPROFILES
	if len(uniqSuperFamilyAssociated) == 1 :
		name = uniqSuperFamilyAssociated[0]
	####	Else MULTIPLE superFamily names are associated to PROTPROFILES, name is undefined with all the PROTPROFILES
	else :
		name="undefined"
		for proteine in PROTPROFILES:
			name += str("_" + proteine)
	return name

def retrieveSuperFamilyAssociated(PROTPROFILE, SUPERFAMILYASSOCIATED):
	"""
	Retrieve all the superFamily names associated to the proteine profile.

	Keyword arguments:
	@type PROTPROFILE: dictionnary
	@param PROTPROFILE: Proteine profile keyword found and its proportion for a given sequence.
	@type SUPERFAMILYASSOCIATED: dictionnary
	@param SUPERFAMILYASSOCIATED: dictionnary containing different superfamily names possible for a given proteine profile.

	@rtype: list
	@return: list of all the superFamily names associated to the proteine profile.
	"""
	####	List that will contain the superFamily name associated to PROTPROFILE
	superFamilyAssociated = []

	for proteine in PROTPROFILE:
		####	Look for the superFamily associated with the PROTPROFILE (check in the SUPERFAMILYASSOCIATED)
		if proteine in SUPERFAMILYASSOCIATED:
			####	Add the superFamily names associated to the PROTPROFILE in a list
			for superFamily in SUPERFAMILYASSOCIATED[proteine]:
				superFamilyAssociated.append(superFamily)
	return superFamilyAssociated

def compareWickerClassification(BLAST, PROTPROFILES, SUPERFAMILYASSOCIATED):
	"""
	Comparison between the blast and the proteines profiles keywords to check if it's in accordance with Wicker classification.
	This classification is implemented in the SUPERFAMILYASSOCIATED dictionnary. Key = blast keyword : Value = protProfile keywords possible

	Keyword arguments:
	@type BLAST: dictionnary
	@param BLAST: All the blast keywords found and their proportion for a given sequence.
	@type PROTPROFILES: dictionnary
	@param PROTPROFILES: All the proteines profiles keywords found and their proportion for a given sequence.
	@type SUPERFAMILYASSOCIATED: dictionnary
	@param SUPERFAMILYASSOCIATED: dictionnary containing different superfamily names possible for a given proteine profile.

	@rtype: string
	@return: Name of the superFamily with differents matches found.
	"""
	####	Check the coherence between specific and non specific keywords
	name = ""
	####	Boolean to know if superFamily name matches with BLAST
	matchBlast = False

	####	List that will contain all the superFamily name associated to PROTPROFILES
	superFamilyMatches = []

	####	First parse PROTPROFILES keyowrds founded
	for proteine in PROTPROFILES:
		####	Look for the superFamily associated with the PROTPROFILES (check in the SUPERFAMILYASSOCIATED)
		if proteine in SUPERFAMILYASSOCIATED:
			####	Add the superFamily names associated to the PROTPROFILES in a list
			for superFamily in SUPERFAMILYASSOCIATED[proteine]:
				####	If the superFamily names associated is in BLAST
				if superFamily in BLAST:
					superFamilyMatches.append(superFamily)
					matchBlast = True
	####	Remove the doublons (ex: Copia, Copia, Mariner => Copia, Mariner)
	superFamilyMatches = list(set(superFamilyMatches))

	####	If the superFamily names associated with PROTPROFILES don't match with the BLAST, name is chimeric
	if not matchBlast:
		####	Retrieve the BLAST name
		for superFamily in BLAST:
			name += "potentialChimeric_" + str(superFamily)
		####	Retrieve all the proteines profiles
		for proteine in PROTPROFILES:
			name += str("_" + proteine)
	####	If the superFamily names associated with PROTPROFILES match with the BLAST, name is superFamily
	else :
		name += str(superFamilyMatches[0])

	return name


def percentageCalculation(KEYWORDFOUND):
	"""
	Calculus of the percentage of each keyword (can be either blast or proteines profiles, choose as arg) found for a given sequence.

	Keyword arguments:
	@type KEYWORDFOUND: dictionnary
	@param KEYWORDFOUND: All the keywords found (can be either blast or proteines) and their occurrence for a given sequence.

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
