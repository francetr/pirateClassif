#!/usr/bin/python3
# coding: utf8
from collections import Counter

""" @author: Tristan Frances """

def superFamilyComparison(SUPERFAMILYFOUND, SUPERFAMILYASSOCIATED):
	"""
	Compare the different superFamilies found for one sequence. For this uses the method combinations of the package itertools allowing to combine 2 by 2 the different elements of a list.

	Keyword arguments:
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
	@param PROTPROFILES: All the proteines profiles keywords found and their count for a given sequence.
	@type SUPERFAMILYASSOCIATED: dictionnary
	@param SUPERFAMILYASSOCIATED: dictionnary containing the different superfamily names possible for a given proteine profile.

	@rtype: string
	@return: Name of the superFamily with different matches found.
	"""
	name = ""
	####	For this part, 7 cases have been identified
	####	First THREE CASES, Consider NO BLAST keywords have been founded
	if len(BLAST) == 0:
		####	NO type of PROPROFILES keyword have been founded
		if len(PROTPROFILES) == 0:
			name = "undefined"

		####	ONE or MULTIPLE type of PROPROFILES keyword have been founded
		elif len(PROTPROFILES) >= 1:
			name = compareProtProfiles(PROTPROFILES, SUPERFAMILYASSOCIATED)

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
	print(name)
	return name

def compareProtProfiles(PROTPROFILES, SUPERFAMILYASSOCIATED):
	"""
	Comparison between the blast and the proteines profiles keywords to check if it's in accordance with Wicker classification.
	This classification is implemented in the SUPERFAMILYASSOCIATED dictionnary. Key = blast keyword : Value = protProfile keywords possible

	Keyword arguments:
	@type PROTPROFILES: dictionnary
	@param PROTPROFILES: All the proteines profiles keywords found and their proportion for a given sequence.
	@type SUPERFAMILYASSOCIATED: dictionnary
	@param SUPERFAMILYASSOCIATED: dictionnary containing the different superfamily names possible for a given proteine profile.

	@rtype: string
	@return: Name of the superFamily.
	"""
	####	APPROACH : Proceed in 2 main steps,
	####	First creates a dictionnary were all the superFamily associated founded for PROTPROFILES are added, with their percentage.
	####	Then add to a list only the superFamily with a percentage >= 75%. Finally: Assigns this name as final name else, as undefined with the matches.

	####	String that will contain the superFamily name
	name = ""
	####	Dictionnary that will contain all the superFamily name associated to PROTPROFILES and their percentage for the sequence
	superFamilyAssociated = retrieveSuperFamilyAssociated(PROTPROFILES, SUPERFAMILYASSOCIATED)

	####	Remove all the superFamily names founded with a percentage lesser than 75 %
	uniqSuperFamilyAssociated = list([k for k in superFamilyAssociated if superFamilyAssociated[k] >= 75])

	####	If only ONE superFamily name is associated to PROTPROFILES
	if len(uniqSuperFamilyAssociated) == 1:
		name = uniqSuperFamilyAssociated[0]
	####	Else MULTIPLE superFamily names are associated to PROTPROFILES, name is undefined with all the PROTPROFILES
	else:
		name = "undefined"
		for proteine in PROTPROFILES:
			name += str("_" + proteine)
	return name

def compareWickerClassification(BLAST, PROTPROFILES, SUPERFAMILYASSOCIATED):
	"""
	Comparison between the blast and the proteines profiles keywords to check if it's in accordance with Wicker classification.
	This classification is implemented in the SUPERFAMILYASSOCIATED dictionnary with Key = blast keyword : Value = protProfile keywords possible

	Keyword arguments:
	@type BLAST: dictionnary
	@param BLAST: All the blast keywords found and their proportion for a given sequence.
	@type PROTPROFILES: dictionnary
	@param PROTPROFILES: All the proteines profiles keywords found and their proportion for a given sequence.
	@type SUPERFAMILYASSOCIATED: dictionnary
	@param SUPERFAMILYASSOCIATED: dictionnary containing the different superfamily names possible for a given proteine profile.

	@rtype: string
	@return: Name of the superFamily with differents matches found.
	"""
	####	APPROACH : Proceed in 2 main steps,
	####	First creates a dictionnary were all the superFamily associated founded for PROTPROFILES are added, with their percentage.
	####	Then add to a list only the superFamily that match with the BLAST and with a percentage >= 75%.
	####	Finally: Assigns this name as final name (if only ONE superFamily have matched); else, as potentialChimeric with the matches (BLAST and PROTPROFILES).

	####	String that will contain the superFamily name
	name = ""

	####	List that will contain all the matches between the BLAST and the superFamily associated.
	superFamilyMatches = []
	####	List that will contain all the superFamily name associated to PROTPROFILES
	superFamilyAssociated = retrieveSuperFamilyAssociated(PROTPROFILES, SUPERFAMILYASSOCIATED)
	for superFamily in superFamilyAssociated.keys():
		####	check if the superFamily associated match with the BLAST AND its percentage is greater than 75%
		if superFamily in BLAST and superFamilyAssociated[superFamily] >= 75:
			superFamilyMatches.append(superFamily)

	####	If the number of matches is 1, name is superFamily
	if len(superFamilyMatches) == 1:
		name += str(superFamilyMatches[0])
	####	If the number of match = 0, name is chimeric
	else :
		name = "potentialChimeric_"
		####	Retrieve the BLAST name
		for superFamily in BLAST:
			name += str(superFamily)
		####	Retrieve all the proteines profiles
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
	@param SUPERFAMILYASSOCIATED: dictionnary containing the different superfamily names possible for a given proteine profile.

	@rtype: dictionnary
	@return: dictionnary of all the superFamily names associated to the proteine profile and their percentage.
	"""
	####	List that will contain the superFamily name associated to PROTPROFILE
	finalSuperFamilyAssociated = {}

	####	Retrieve all the PROTPROFILES
	for proteine in PROTPROFILE:
		####	Look for the superFamily associated with the PROTPROFILE (check in the SUPERFAMILYASSOCIATED)
		if proteine in SUPERFAMILYASSOCIATED:
			####	Add the superFamily names associated to the PROTPROFILE in a list
			for superFamily in SUPERFAMILYASSOCIATED[proteine]:
				####	If the superFamily associated is not in the list of finalSuperFamilyAssociated, add this superFamily in the list with its percentage
				if not superFamily in finalSuperFamilyAssociated.keys():
					finalSuperFamilyAssociated[superFamily] = PROTPROFILE[proteine]
				####	If the superFamily associated is already in the list of finalSuperFamilyAssociated, add this superFamily in the list with its percentage
				else:
					finalSuperFamilyAssociated[superFamily] += PROTPROFILE[proteine]

	return finalSuperFamilyAssociated
