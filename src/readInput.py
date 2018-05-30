#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

import sys
import argparse
import re
import categorization
from Bio import SeqIO # For the fasta reading


def retrieveArguments():
	"""
	Retrieve the arguments from the command line.

	@rtype: argParse
	@return: attributes NameSpace of object argParse that contain the two arguments name.
	"""
	####	 Mananage the 2 arguments (PASTEC and FASTA file name) when the command is launched
	parser = argparse.ArgumentParser(prog="scriptClassif.py", description="This program is a part of the PiRATE project. It aims to automatized the step of TE classification")
	parser.add_argument("classif", type=str, help="classif file that comes from PASTEC")
	parser.add_argument("fasta", type=str, help="fasta file providing the sequence")
	parser.add_argument("-i", metavar="IDENTITY", type=int, nargs="?", default=100, help="Threshold for considering two sequences as identical, enter an integer from 0 to 100, default is 100.")
	parser.add_argument("--baseline", type=str, nargs="?", default="base_reference.txt", help="baseline file giving the different names possible for a superfamily")
	args = parser.parse_args()
	# checkArguments(args.classif, args.fasta)
	return args


def readPastec(PASTEC, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD):
	"""
	Read a PASTEC file as input line by line, and then proceed to the categorization of each sequence.

	Keyword arguments:
	@type PASTEC: string
	@param PASTEC: name of the classif file that will be opened.
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)
	@type BASELINE: dictionnary
	@param BASELINE: dictionnary containing different superfamily names possible for a given superfamily (usefull for the function superFamilyComparison).
	@type IDENTITYTHRESHOLD: integer
	@param IDENTITYTHRESHOLD: percentage for which the superFamily name of a sequence will be choosen.

	@rtype: None
	"""
	try:
		####	Open the classif file
		with open(PASTEC, "r") as f:
			print("Begin of the categorization\n")
			####	Parse every line/sequence of the file
			for line in f:
				sequence=line.replace("\t\t", "\t").strip() # remove the first column double tab and carriage return of the line
				categorization.initCategorization(sequence, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD)
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
	{B{Key} = blast : {B{Key} = name of blast keywords for superFamily used later : B{I{Values}} = [list of possible names for this superFamily]},
	B{Key} = protProfiles :{B{Key} = name of proteines profiles keyword for superFamily used later : B{I{Values}} = [list of possible names for this superFamily]}}.
	"""
	baselineDictionnary={"blast":{},"protProfiles":{}}
	try:
		with open(BASELINE, "r") as f:
			for line in f:
				####	check the content of the file and look for > character (It indicates if the keyword is specific or nonSpecific)
				keywords=re.search(r'>([^\n]+)', line)
				####	For specific keywords
				if not keywords:
					readBlastBaselineKeywords(line, baselineDictionnary)
				####	For nonSpecific keyword
				else:
					readProtProfilesBaselineKeywords(keywords.groups()[0], baselineDictionnary)
		return baselineDictionnary

	except FileNotFoundError:
		print("/!\	Error: No such file as {}\n####	Classification aborted".format(BASELINE))
		sys.exit(1)


def readBlastBaselineKeywords(SPECIFICKEYWORDS, BASELINEDICTIONNARY):
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
	BASELINEDICTIONNARY["blast"][superFamilyNames[0]]=listPossibleNames

def readProtProfilesBaselineKeywords(NONSPECIFICKEYWORDS, BASELINEDICTIONNARY):
	"""
	From one line of the BASELINE file, it add the non specific keyword founded as a key of the BASELINEDICTIONNARY and all the possible names associated as the value.

	Keyword argument:
	@type NONSPECIFICKEYWORDS: string
	@param NONSPECIFICKEYWORDS: One line read of the BASELINE file (here doesn't have the > character), which contains the proteines profiles keywords..
	@type BASELINEDICTIONNARY: dictionnary
	@param BASELINEDICTIONNARY: Dictionnary containing all the keywords founded in the BASELINE file.
	It contains 2 key, one for the storage of blast keywords and the other for the storage of proteines profiles keywords.


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
	BASELINEDICTIONNARY["protProfiles"][superFamilyNames[0]]=listPossibleNames

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
		with open(FASTA, "rU") as handle:
		####	Parse every line/sequence of the file
			for record in SeqIO.parse(handle, "fasta"):
				####	Save the sequence in a dictionnary
				seqReturned[record.id]={"seq":record.seq}
		return seqReturned
	except (FileNotFoundError, NameError):
		####	prevent the opening if the file name is incorrect
		print("/!\	Error: No such file as {}\n####	Classification aborted".format(FASTA))
		sys.exit(1)

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
