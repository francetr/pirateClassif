#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

import sys
import argparse
import re
import categorization


def retrieveArguments():
	"""
	Retrieve the arguments from the command line.

	@rtype: argParse
	@return: attributes NameSpace of object argParse that contain the two arguments name.
	"""
	####	 Mananage the 2 arguments (PASTEC and FASTA file name) when the command is launched
	parser = argparse.ArgumentParser(prog="scriptClassif.py", description="This program is a part of the PiRATE project. It aims to automatized the step of TE classification")
	parser.add_argument("classif", type=str, help="classif file that comes from PASTEC")
	parser.add_argument("baseline", type=str, help="baseline file giving the different names possible for a superfamily")
	parser.add_argument("fasta", type=str, help="fasta file providing the sequence")
	args = parser.parse_args()
	# checkArguments(args.classif, args.fasta)
	return args

# def checkArguments(CLASSIFNAME, FASTANAME):
# 	"""
# 	Check the extension of the two arguments of command line. If it is not correct, stop the programm.
#
# 	Keyword argument
# 	@param CLASSIFNAME: name of the argument for the CLASSIF file
# 	@param FASTANAME: name of the argument for the FASTA file
#
# 	@return: None
#
# 	"""
# 	try:
# 		classifName=re.match(r'[\S]*[/]?[\w\-.]+(classif)$',CLASSIFNAME).groups()[0] ####	regex checking is classif file as good extension
# 	except AttributeError as e:
# 		print("####	Wrong extension for the classif file\n####	 Classification aborted")
# 		sys.exit(1)
#
# 	try:
# 		fastaName=re.match(r'[\S]*[/]?[\w\-.]+(fasta)$', FASTANAME).groups()[0] ####	regex checking is fasta file as good extension
# 	except AttributeError as e:
# 		print("####	Wrong extension for the fasta file\n####	Classification aborted")
# 		sys.exit(1)

def readPastec(PASTEC, NONTE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE):
	"""
	Read a PASTEC file as input line by line, and then proceed to the categorization of each sequence.

	Keyword arguments:
	@type PASTEC: string
	@param PASTEC: name of the classif file that will be opened.
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
	try:
		####	Open the classif file
		with open(PASTEC, "r") as f:
			print("Begin of the categorization\n")
			####	Parse every line/sequence of the file
			for line in f:
				sequence=line.replace("\t\t", "\t").strip() # remove the first column double tab and carriage return of the line
				categorization.initCategorization(sequence, NONTE, POTENTIALCHIMERIC, NOCAT, TE, BASELINE)
				# print("Summary, sequences find by type :\nTE : {}\nnonTE : {}\nnoCat : {}\npotentialChimeric : {}\nTotal : {}".format(len(TE),\
				# len(nonTE), len(noCat), len(potentialChimeric), (len(TE)+len(nonTE)+len(noCat)+len(potentialChimeric))))
			print("End of the categorization\n")

	except FileNotFoundError:
		####	Prevent the opening if the file name is incorrect
		print("/!\	Error: No such file as {}\n####	Classification aborted".format(PASTEC))
		sys.exit(1)

def readBaseline(BASELINE):
	"""
	Read the BASELINE file as input line by line. First, search if the beginning of line contains the character > typical of non specific keyword.
	If not, it call the function dedicated to search specific keywords, else, it calls the functions for non specific keywords
	Then the keywords are then added to the dictionnary

	Return a dictionnary which contains the sequence and the id of each sequences which are in the FASTA file.

	Keyword argument:
	@type BASELINE: string
	@param BASELINE: name of the baseline file that will be opened.

	@rtype: dictionnary
	@return: dictionnary containing the different names of possible superFamily.
	{B{Key} = specific : {B{Key} = name of specific keywords for superFamily used later : B{I{Values}} = [list of possible names for this superFamily]},
	B{Key} = nonSpecific :{B{Key} = name of non specific keyword for superFamily used later : B{I{Values}} = [list of possible names for this superFamily]}}.
	"""
	baselineDictionnary={"specific":{},"nonSpecific":{}}
	try:
		with open(BASELINE, "r") as f:
			for line in f:
				####	check the content of the file and look for > character (It indicates if the keyword is specific or nonSpecific)
				keywords=re.search(r'>([^\n]+)', line)
				####	For specific keywords
				if not keywords:
					readSpecificBaselineKeywords(line, baselineDictionnary)
				####	For nonSpecific keyword
				else:
					readNonSpecificBaselineKeywords(keywords.groups()[0], baselineDictionnary)

	except FileNotFoundError:
		print("/!\	Error: No such file as {}\n####	Classification aborted".format(BASELINE))
	return baselineDictionnary

def readSpecificBaselineKeywords(SPECIFICKEYWORDS, BASELINEDICTIONNARY):
	"""
	From one line of the BASELINE file, it add the specific keyword founded as a key of the BASELINEDICTIONNARY and all the possible names associated as the value.

	Keyword argument:
	@type SPECIFICKEYWORDS: string
	@param SPECIFICKEYWORDS: One line read of the BASELINE file, which contains the specific keywords.
	@type BASELINEDICTIONNARY: dictionnary
	@param BASELINEDICTIONNARY: Dictionnary containing all the keywords founded in the BASELINE fileself.
	It contains 2 key, one for the storage of specific keywords and the other for the storage of specific keywords.


	@rtype: None
	"""
	####	List that will contains the different possibles names of a superfamily
	listPossibleNames=[]
	####	Recover the reference name for the superfamily
	superFamilyNames = SPECIFICKEYWORDS.split("\t")
	####	Recover the possible names for the superfamily
	allPossibleNames = superFamilyNames[1].strip().split(":")
	for possibleName in allPossibleNames:
		####	Add the possibles names for a superFamily in the list
		listPossibleNames.append(possibleName)
	####	Complete the dictionnary with : Key = name of superfamily used later; Value = list of possible names for this superfamily
	BASELINEDICTIONNARY["specific"][superFamilyNames[0]]=listPossibleNames

def readNonSpecificBaselineKeywords(NONSPECIFICKEYWORDS, BASELINEDICTIONNARY):
	"""
	From one line of the BASELINE file, it add the non specific keyword founded as a key of the BASELINEDICTIONNARY and all the possible names associated as the value.

	Keyword argument:
	@type NONSPECIFICKEYWORDS: string
	@param NONSPECIFICKEYWORDS: One line read of the BASELINE file (doesn't have the > character), which contains the non specific keywords..
	@type BASELINEDICTIONNARY: dictionnary
	@param BASELINEDICTIONNARY: Dictionnary containing all the keywords founded in the BASELINE fileself.
	It contains 2 key, one for the storage of specific keywords and the other for the storage of specific keywords.


	@rtype: None
	"""
	####	List that will contains the different possibles names of a superfamily
	listPossibleNames=[]
	####	Recover the reference name for the superfamily
	superFamilyNames = NONSPECIFICKEYWORDS.split("\t")
	####	Recover the possible names for the superfamily
	allPossibleNames = superFamilyNames[1].strip().split(":")
	for possibleName in allPossibleNames:
		####	Add the possibles names for a superFamily in the list
		listPossibleNames.append(possibleName)
	####	Complete the dictionnary with : Key = name of superfamily used later; Value = list of possible names for this superfamily
	BASELINEDICTIONNARY["nonSpecific"][superFamilyNames[0]]=listPossibleNames				# ####	List that will contains the different possibles names of a superfamily


def readFasta(FASTA):
	"""
	Read a FASTA file as input and return it as a string. Use the SeqIO from the package Bio of BioPython.

	Return a dictionnary which contains the sequence and the id of each sequences that are in the FASTA file.

	Keyword argument:
	@type FASTA: string
	@param FASTA: name of the FASTA file containing the sequence that will be opened.

	@rtype: dictionnary
	@return: Dictionnary with {B{key} = id of the sequence and value : B{I{value}} = {"seq" : sequence of the FASTA sequence} }
		- id : id of the sequence
		- name : name of the sequence
		- description : description of the sequence
		- number of features : number of features of the sequence
		- seq : sequence concerned
	"""
	seqReturned={}
	try:
		####	 Open the fasta file
		with open(fastaFile, "rU") as handle:
			print(handle)
		####	Parse every line/sequence of the file
			for record in SeqIO.parse(handle, "fasta"):
				####	Save the sequence in a dictionnary
				seqReturned[record.id]={"seq":record.seq}
		return seqReturned
	except (FileNotFoundError, NameError):
		####	prevent the opening if the file name is incorrect
		print("/!\	Error: No such file as {}\n####	Classification aborted".format(FASTA))
		sys.exit(1)
