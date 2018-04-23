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
	####	Comparison between names found for a sequence and names in the BASELINE
	####	Run through the list of superFamily found
	name = ""
	specificResult = specificComparison(FEATURES, SUPERFAMILYFOUND, BASELINE)
	# print(specificResult)
	nonSpecificComparison(FEATURES, SUPERFAMILYFOUND, BASELINE, specificResult)

	# for superFamily in SUPERFAMILYFOUND:
	# 	####	Run through the reference name of the BASELINE
	# 	for possibleName in BASELINE["specific"].keys():
	#
	# 		####	If the superfamily name is not in the BASELINE AND is not a key of the dictionnary count, its string is declared as value of this dictionnary and has its counter equal to 1
	# 		if superFamily in BASELINE["specific"][possibleName] and (not(superFamily in superFamilyCount.keys())):
	# 			superFamilyCount[superFamily]=1
	# 		####	If the superfamily name is already in the dictionnary and exists in the BASELINE, its counter is incremented by 1
	# 		elif (superFamily in BASELINE["specific"][possibleName]) and (superFamily in superFamilyCount.keys()):
	# 			superFamilyCount[superFamily]+=1
	# print(FEATURES[0], superFamilyCount)

	# name=""
	# max = 0
	# ####	Parse the keys of the superfamily name dictionnary count
	# for superFamily in superFamilyCount.keys():
	# 	####	If the count of superfamily name is superior to max, max is equal to this count and name is the superFamily name
	# 	if superFamilyCount[superFamily] > max:
	# 		max=superFamilyCount[superFamily]
	# 		name=superFamily
	# 	####	Else if the count is equal to the max
	# 	elif superFamilyCount[superFamily] == max:
	# 		max=superFamilyCount[superFamily]
	# 		name="potentialChimeric"
	# #### if the name is ? convert it into unknown
	# if name=="?":
	# 	name="unknown"

	# print(FEATURES[0], name, max)
	return name
	# print("{}, {}, count: {}, len: {}".format(FEATURES[0], SUPERFAMILYFOUND, superFamilyCount, len(SUPERFAMILYFOUND)))

def specificComparison(FEATURES, SUPERFAMILYFOUND, BASELINE):
	"""
	Search if there one or multiple specific keyword(s) founded for one sequence.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type SUPERFAMILYFOUND: dictionnary
	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily.

	@rtype: tuple
	@return: name of the superfamily corresponding to the sequence and number of specific keyword(s) founded for this sequence.
	"""
	#TODO
	####	If one specific keyword has been found
	if len(SUPERFAMILYFOUND["specific"]) == 1:
		for key in SUPERFAMILYFOUND["specific"].keys():
			nameSuperfamily = key
			return {nameSuperfamily: SUPERFAMILYFOUND["specific"][key]}

	####	If multple specific keywords have been found
	elif len(SUPERFAMILYFOUND["specific"]) >= 1:
		return {"potentialChimeric": len(SUPERFAMILYFOUND["specific"])}

def nonSpecificComparison(FEATURES, SUPERFAMILYFOUND, BASELINE, SPECIFICRESULT):
	"""
	Search if there one or multiple non specific keyword(s) founded for one sequence.

	Keyword arguments:
	@type FEATURES: list
	@param FEATURES: names of the features (potentialChimeric, class, order, ...) find in the sequence.
	@type SUPERFAMILYFOUND: dictionnary
	@param SUPERFAMILYFOUND: count the names of the superfamily found for one sequence during the superFamilyDetermination.
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily.
	@type SPECIFICRESULT: TODO
	@param SPECIFICRESULT: TODO):

	@rtype: tuple
	@return: name of the superfamily corresponding to the sequence and number of non specific keyword(s) founded for this sequence.
	"""
	#TODO
	pass
