#!/usr/bin/python2
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
	parser.add_argument("-e", metavar="IDENTITY", type=float, nargs="?", default=100.0, help="Threshold for considering two sequences as identical, enter an integer from 0 to 100, default is 100.")
	parser.add_argument("--baseline", type=str, nargs="?", default="base_reference.txt", help="baseline file giving the different names possible for a superfamily")
	parser.add_argument("--config", type=str, nargs="?", default="Filter_pirateClassif.txt", help="config file for the creation of libraries")
	args = parser.parse_args()
	####	If identity threshold > 100 convert it to the value 100
	if args.e >= 100:
		args.e = 100
	####	If identity threshold =< 0 convert it to the value 0
	elif args.e <= 0:
		args.e = 0
	# checkArguments(args.classif, args.fasta)
	return args


def readPastec(PASTEC, SEQCLASSIFIED, BASELINE, IDENTITYTHRESHOLD):
	"""
	Read a PASTEC file as input line by line, and then proceed to the categorization of each sequence.

	Keyword arguments:
	@type PASTEC: string
	@param PASTEC: name of the classif file that will be opened.
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

	except (FileNotFoundError, NameError) as e:
		####	Prevent the opening if the file name is incorrect
		print("/!\	Error: {}\n####	Classification aborted".format(e))
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

	except (FileNotFoundError, NameError) as e:
		print("/!\	Error: {}\n####	Classification aborted".format(e))
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
	It contains 2 keys, one for the storage of blast keywords and the other for the storage of proteines profiles keywords.

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

	Return a dictionnary which contains the sequence, the id and the length of each sequences that are in the FASTA file.

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
		- length : length of the sequence

	"""
	seqReturned={}
	try:
		####	 Open the fasta file
		with open(FASTA, "rU") as handle:
		####	Parse every line/sequence of the file
			for record in SeqIO.parse(handle, "fasta"):
				####	Save the sequence in a dictionnary
				seqReturned[record.id]={"seq":record.seq.upper(), "length":len(record.seq)}
		return seqReturned
	except (FileNotFoundError, NameError) as e:
		####	prevent the opening if the file name is incorrect
		print("/!\	Error: {}\n####	Classification aborted".format(e))
		sys.exit(1)

def readConfig(CONFIG):
	"""
	Read the config file as input and return it as a dictionnary.

	Return a dictionnary which contains the sequence and the id of each sequences that are in the FASTA file.

	Keyword argument:
	@type CONFIG: string
	@param CONFIG: name of the config file containing the filters to apply in order to create the libraries.

	@rtype: dictionnary
	@return: dictionnary with parameters to filter sequences (which will be used to create final libraries) {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
		There are currently 8 keys for this dictionnary:
		- outputName : name of the output file
		- lengthMin : minimal length of the sequence
		- lengthMax : maximum length of the sequence
		- removedTool : name of the tool used on the sequence. Remove the sequences found by this tool
		- onlySelectedtools : name of the tool used on the sequence. Select the sequences found by this tool
		- autonomousLib : boolean (yes/no) indicating if sequence constitute autonomous librarie
		- totalTELib : boolean (yes/no) indicating if sequence constitute total TEs librarie
		- totalRepeatLib : boolean (yes/no) indicating if sequence constitute repeated elements librarie
	"""
	configFilter={}
	try:
		####	Open the classif file
		with open(CONFIG, "r") as f:
			print("Begin of the categorization\n")
			#### Skip the header row
			next(f)
			####	Parse every line/sequence of the file
			for line in f:
				line=line.replace("\t\t", "\t").strip() # remove the first column double tab and carriage return of the line
				line=line.split("\t")
				#### Defintion of the filter onto the config dictionnary
				configFilter[line[0]]={}
				configFilter[line[0]]["outputName"]=line[1]
				configFilter[line[0]]["lengthMin"]=int(line[2].split(":")[0])
				configFilter[line[0]]["lengthMax"]=int(line[2].split(":")[1])
				if line[3].lower() == "na" :
					configFilter[line[0]]["removedTool"]=line[3]
				else:
					configFilter[line[0]]["removedTool"]=[tool for tool in line[3].split(":")]
				#### Take into account if multiple tools areused as filter
				if line[4].lower() == "na" :
					configFilter[line[0]]["onlySelectedtools"]=line[4]
				else:
					configFilter[line[0]]["onlySelectedtools"]=[tool for tool in line[4].split(":")]
				configFilter[line[0]]["autonomousLib"]=line[5]
				configFilter[line[0]]["totalTELib"]=line[6]
				configFilter[line[0]]["totalRepeatLib"]=line[7]

			return configFilter
	except (FileNotFoundError, NameError) as e:
		####	Prevent the opening if the file name is incorrect
		print("/!\	Error: {}\n####	Read config file aborted".format(e))
		sys.exit(1)
