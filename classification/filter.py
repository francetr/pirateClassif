#!/usr/bin/python2
# coding: utf8

import os
import subprocess
import save
import readInput
import re

def initFilters(CONFIG):
	"""
	Save filtered sequences according the config file.

	Keyword arguments:
	@type CONFIG: dictionnary
	@param CONFIG:  dictionnary with parameters to filter sequences (which will be used to create final libraries) {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
		There are currently 8 keys for this dictionnary:
			- finalClassification : final Classification given to the sequence by the classification script
			- lengthMin : minimal length of the sequence
			- lengthMax : maximum length of the sequence
			- removedTools : name of the tool used on the sequence. Remove the sequences found by this tool
			- onlySelectedtools : name of the tool used on the sequence. Select only the sequences found by this tool
			- autonomousLib : boolean (yes/no) indicating if sequence constitute autonomous librarie
			- totalTELib : boolean (yes/no) indicating if sequence constitute total TEs librarie
			- totalRepeatLib : boolean (yes/no) indicating if sequence constitute repeated elements librarie

	@rtype: None
	"""
	#### String containing all the preLibraries file name
	preLibraries = findFile("classification_result/prelibraries/TE", "*.fasta")

	#### String containing all the preLibraries file name
	noCatLibrarie = findFile("classification_result/prelibraries/", "noCat.fasta")

	listPrelibraries = []
	#### dictionnaries that will contains all the id's sequences for concerned libraries
	dicoLibraries={"autonomousLib":[], "totalTELib":[], "totalRepeatLib":[]}

	listPrelibraries.append(noCatLibrarie[0])
	#### Add all the name of prelibraries in listPrelibraries
	for file in preLibraries:
		listPrelibraries.append(file)

	#### Dictionnary that restain the final classification for a given sequence (helpfull for the intermediateLibraries)
	dicoFinalClassif={}
	#### Parse all the prelibrary
	print("####	Apply the filters to create the intermediate libraries")
	createIntermediateLibraries(listPrelibraries, dicoLibraries, CONFIG, dicoFinalClassif)

	#### List containing all the intermediate librarie file name
	intermediateLibraries = findFile("classification_result/intermediateLibraries", "*.fasta")

	print("####	Apply the cd-hit-est on the intermediate libraries")
	applyCDHIT(intermediateLibraries)

	retriveFinalLibrarieSequences(intermediateLibraries, CONFIG, dicoFinalClassif, dicoLibraries)

	print("####	Creation of the three final libraries")
	createFinalLibraries(intermediateLibraries, dicoLibraries)

	print("Number of sequences in autonomousTE : {nbAutonomous}\nNumber of sequences in totalTE : {nbTotalTE}\nNumber of sequences in totalRepeatLib : {nbRepeated}".format(\
	nbAutonomous=len(dicoLibraries["autonomousLib"]), nbTotalTE=len(dicoLibraries["totalTELib"]), nbRepeated=len(dicoLibraries["totalRepeatLib"])))

def findFile(PATH, PATTERN):
	"""
	Find files according a patern. Return a list with the absolute path of all the files name matching with the pattern.

	@type PATH: string
	@param PATH:  string of the absolute path where the file will be searched
	@type PATTERN: string
	@param PATTERN:  string with parameters to filter sequences (which will be used to create final libraries)

	@rtype: list
	@return: list of the absolute path towrd all the file matching with the pattern
	"""
	#### find all the file located in the PATH corresponding to the PATTERN
	findPattern = subprocess.Popen('find {path} -name "{pattern}" 2> /dev/null'.format(path=PATH, pattern=PATTERN), shell=True, stdout=subprocess.PIPE);
	#### String containing all the path toward the file's name matching with the pattern
	return findPattern.stdout.read().decode("utf-8").strip().split("\n")

def createIntermediateLibraries(LISTPRELIBRARIES, DICOLIBRARIES, CONFIG, DICOFINALCLASSIF):
	"""
	Creates the three final Libraries.

	@type LISTPRELIBRARIES: list
	@param LISTPRELIBRARIES: list of the path to the preLibraries files
	@type DICOLIBRARIES: dictionnary
	@param DICOLIBRARIES: dictionnary that will contain the id of the sequences to construct the 3 final libraries {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
		- listAutonomousTE
		- listTotalTE
		- listRepeated
	@type CONFIG: dictionnary
	@param CONFIG:  dictionnary with parameters to filter sequences (which will be used to create final libraries) {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
	@type DICOFINALCLASSIF: dictionnary
	@param DICOFINALCLASSIF: dictionnary that will contain the id of the sequences and its corresponding finalClassification {B{id} : B{I{finalClassification of the sequence}}}

	@rtype: None
	"""
	#### Parse all the intermediate libraries files
	for preLibrary in LISTPRELIBRARIES:
		#### Retrieve the final classification name of the ET from the file name
		finalClassification = os.path.basename(preLibrary).split(".fasta")[0]
		#### Read and store the fasta sequences of the prelibraries
		sequences=readInput.readFasta(preLibrary)
		#### Parse all the sequences
		for id in sequences:
			#### Check the finalClassification of the sequences is in the ID
			if finalClassification.lower() in id.lower():
				DICOFINALCLASSIF[id]=finalClassification
				applyFiltersForIntermediate(id, sequences, finalClassification, CONFIG, DICOLIBRARIES)

def applyFiltersForIntermediate(ID, SEQUENCES, FINALCLASSIFICATION, CONFIG, DICOLIBRARIES):
	"""
	Apply filters from the CONFIG file to creates the final libraries.

	@type ID: string
	@param ID: string of the id of the sequence
	@type SEQUENCES: dictionnary
	@param SEQUENCES: dictionnary containing the fasta sequence of the preLibraries
	@type CONFIG: dictionnary
	@param CONFIG:  dictionnary with parameters to filter sequences (which will be used to create final libraries) {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
	@type FINALCLASSIFICATION: string
	@param FINALCLASSIFICATION: string indicating the final classification of the sequence (usefull to compare SEQUENCES with the CONFIG file)
	@type DICOLIBRARIES: dictionnary
	@param DICOLIBRARIES: dictionnary that will contain the id of the sequences to construct the 3 final libraries {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
		- listAutonomousTE
		- listTotalTE
		- listRepeated

	@rtype: None
	"""
	#### First we check if tools used to detect the TE are from the removedTool or the onlySelectedtools of the CONFIG file
	#### Parse the dictionnary containing CONFIG information
	for sortedOutput in CONFIG:
		#### Look for te finalClassification corresponding to the sequence
		if CONFIG[sortedOutput]["finalClassification"].lower() == FINALCLASSIFICATION.lower():
			#### Find if the tool used to find the TE is a part of the removedTool from the CONFIG file : True it is find, else False
			findRemovedTool=False
			for removedTool in CONFIG[sortedOutput]["removedTools"]:
				if ID.lower().find(removedTool.lower()) == 0:
					findRemovedTool=True

			#### Find if the tool used to find the TE is one of the onlySelectedtools from the CONFIG file : True it is find, else False
			findSelectedTool=False
			for selectedTool in CONFIG[sortedOutput]["onlySelectedtools"]:
				if ID.lower().find(selectedTool.lower()) == 0 or selectedTool.lower() == "na":
					findSelectedTool=True

			#### Control the length of the sequence with the CONFIG
			if SEQUENCES[ID]["length"] >= CONFIG[sortedOutput]["lengthMin"] and SEQUENCES[ID]["length"] < CONFIG[sortedOutput]["lengthMax"]:
				#### Control if the tool used to detect TE is in the removedTool : if not, sequence can be in library autonomous TE
				if not findRemovedTool and findSelectedTool:
						#### Save intermediate libraries
						save.saveIntermadiateLibraries(SEQUENCES, ID, CONFIG, sortedOutput)


def applyCDHIT(INTERMEDIATELIBRARIES):
	"""
	Apply the cd-hit-est on all the intermediate Libraries.

	@type INTERMEDIATELIBRARIES: list
	@param INTERMEDIATELIBRARIES: list of the path to the intermediateLibraries files

	@rtype: None
	"""
	#### Apply cd-hit-est for all the intermediate library
	for file in INTERMEDIATELIBRARIES:
		fileName = os.path.basename(file).split(".fasta")[0]
		os.chdir("classification_result/intermediateLibraries/")
		subprocess.call('cd-hit-est -aS 0.9 -c 0.9 -g 1 -r 1 -i {input}.fasta -o {output}.fasta_tmp'.format(input=fileName, output=fileName), shell=True)
		subprocess.call("mv {input}.fasta_tmp {output}.fasta".format(input=fileName, output=fileName), shell=True)
		os.chdir("../..")

def retriveFinalLibrarieSequences(INTERMEDIATELIBRARIES, CONFIG, DICOFINALCLASSIF, DICOLIBRARIES):
	"""
	Add the id of the sequences onto the DICOLIBRARIES in order to know for which library a sequence will be saved.

	@type INTERMEDIATELIBRARIES: list
	@param INTERMEDIATELIBRARIES: list of the path to the intermediate libraries files
	@type CONFIG: dictionnary
	@param CONFIG:  dictionnary with parameters to filter sequences (which will be used to create final libraries) {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
	@type DICOFINALCLASSIF: dictionnary
	@param DICOFINALCLASSIF: dictionnary that will contain the id of the sequences and its corresponding finalClassification {B{id} : B{I{finalClassification of the sequence}}}
	@type DICOLIBRARIES: dictionnary
	@param DICOLIBRARIES: dictionnary that will contain the id of the sequences to construct the 3 final libraries {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
		- listAutonomousTE
		- listTotalTE
		- listRepeated

	@rtype: None
	"""
	for file in INTERMEDIATELIBRARIES:
		fileName = os.path.basename(file).split(".fasta")[0]
		#### Read and store the fasta sequences of the prelibraries
		sequences=readInput.readFasta(file)
		for id in sequences:
			if fileName in CONFIG:
				applyFiltersForFinales(id, sequences, DICOFINALCLASSIF[id], CONFIG, fileName, DICOLIBRARIES)

def applyFiltersForFinales(ID, SEQUENCES, FINALCLASSIFICATION, CONFIG, SORTEDOUTPUT, DICOLIBRARIES):
	"""
	Apply filters from the CONFIG file to creates the final libraries.

	@type ID: string
	@param ID: string of the id of the sequence
	@type SEQUENCES: dictionnary
	@param SEQUENCES: dictionnary containing the fasta sequence of the preLibraries
	@type CONFIG: dictionnary
	@param CONFIG:  dictionnary with parameters to filter sequences (which will be used to create final libraries) {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
	@type FINALCLASSIFICATION: string
	@param FINALCLASSIFICATION: string indicating the final classification of the sequence (usefull to compare SEQUENCES with the CONFIG file)
	@type SORTEDOUTPUT: string
	@param SORTEDOUTPUT: string indicating the name of the sorted output file
	@type DICOLIBRARIES: dictionnary
	@param DICOLIBRARIES: dictionnary that will contain the id of the sequences to construct the 3 final libraries {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
		- listAutonomousTE
		- listTotalTE
		- listRepeated

	@rtype: None
	"""
	#### First we check if tools used to detect the TE are from the removedTool or the onlySelectedtools of the CONFIG file
	#### Look for te finalClassification corresponding to the sequence
	if CONFIG[SORTEDOUTPUT]["finalClassification"].lower() == FINALCLASSIFICATION.lower():
		#### Find if the tool used to find the TE is a part of the removedTool from the CONFIG file : True it is find, else False
		findRemovedTool=False
		for removedTool in CONFIG[SORTEDOUTPUT]["removedTools"]:
			if ID.lower().find(removedTool.lower()) == 0:
				findRemovedTool=True

		#### Find if the tool used to find the TE is one of the onlySelectedtools from the CONFIG file : True it is find, else False
		findSelectedTool=False
		for selectedTool in CONFIG[SORTEDOUTPUT]["onlySelectedtools"]:
			if ID.lower().find(selectedTool.lower()) == 0 or selectedTool.lower() == "na":
				findSelectedTool=True

		#### Control the length of the sequence with the CONFIG
		if SEQUENCES[ID]["length"] >= CONFIG[SORTEDOUTPUT]["lengthMin"] and SEQUENCES[ID]["length"] < CONFIG[SORTEDOUTPUT]["lengthMax"]:
			#### Control if the tool used to detect TE is in the removedTool : if not, sequence can be in library autonomous TE
			if not findRemovedTool and findSelectedTool:
				#### Control in which libraries the sequence can be added
				#### List of ID for construct librarie of potential autonomous TE
				if CONFIG[SORTEDOUTPUT]["autonomousLib"].lower()=="yes":
					DICOLIBRARIES["autonomousLib"].append(ID)
				else:
					#### Save the id of non selected sequences for this library
					save.saveRemovedSequences(ID, "autonomousLib")
				#### List of ID for construct librarie of total TE
				if CONFIG[SORTEDOUTPUT]["totalTELib"].lower()=="yes":
					DICOLIBRARIES["totalTELib"].append(ID)
				else:
					#### Save the id of non selected sequences for this library
					save.saveRemovedSequences(ID, "totalTELib")
				#### List of ID for construct librarie of repeated Elements
				if CONFIG[SORTEDOUTPUT]["totalRepeatLib"].lower()=="yes":
					DICOLIBRARIES["totalRepeatLib"].append(ID)

def createFinalLibraries(INTERMEDIATELIBRARIES, DICOLIBRARIES):
	"""
	Creates the three final Libraries.

	@type INTERMEDIATELIBRARIES: list
	@param INTERMEDIATELIBRARIES: list of the path to the intermediateLibraries files
	@type DICOLIBRARIES: dictionnary
	@param DICOLIBRARIES: dictionnary that will contain the id of the sequences to construct the 3 final libraries {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
		- listAutonomousTE
		- listTotalTE
		- listRepeated

	@rtype: None
	"""
	#### Parse all the intermediate libraries files
	for file in INTERMEDIATELIBRARIES:
		fileName = os.path.basename(file).split(".fasta")[0]
		#### Read and store the fasta sequences of the prelibraries
		sequences=readInput.readFasta(file)
		#### Save the three finals libraries
		save.saveLibraries(sequences, DICOLIBRARIES)
