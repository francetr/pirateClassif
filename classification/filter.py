#!/usr/bin/python2
# coding: utf8

import os
import subprocess
import save
import readInput

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

	#### Parse all the prelibrary
	print("####	Apply the filters")
	for preLibrary in listPrelibraries:
		#### Retrieve the final classification name of the ET from the file name
		finalClassification = os.path.basename(preLibrary).split(".fasta")[0]
		#### Read and store the fasta sequences of the prelibraries
		sequences=readInput.readFasta(preLibrary)
		#### Parse all the sequences
		for id in sequences:
			#### Check the finalClassification of the sequences is in the ID
			if finalClassification.lower() in id.lower():
				applyFilters(id, sequences, finalClassification, CONFIG, dicoLibraries)

	#### List containing all the intermediate librarie file name
	intermediateLibraries = findFile("classification_result/intermediateLibraries", "*.fasta")
	#### Apply cd-hit-est for all the intermediate library
	for file in intermediateLibraries:
		fileName = os.path.basename(file).split(".fasta")[0]
		os.chdir("classification_result/intermediateLibraries/")
		subprocess.call('cd-hit-est -aS 0.9 -c 0.9 -g 1 -r 1 -i {input}.fasta -o {output}.fasta_tmp'.format(input=fileName, output=fileName), shell=True)
		subprocess.call("mv {input}.fasta_tmp {output}.fasta".format(input=fileName, output=fileName), shell=True)
		os.chdir("../..")

	print("####	Creation of the three final libraries")
	for file in intermediateLibraries:
		#### Read and store the fasta sequences of the prelibraries
		sequences=readInput.readFasta(file)
		#### Save the three finals libraries
		save.saveLibraries(sequences, dicoLibraries)
	print("Number of sequences in autonomousTE : {nbAutonomous}\nNumber of sequences in totalTE : {nbTotalTE}\nNumber of sequences in totalRepeatLib : {nbRepeated}".format(\
	nbAutonomous=len(dicoLibraries["autonomousLib"]), nbTotalTE=len(dicoLibraries["totalTELib"]), nbRepeated=len(dicoLibraries["totalRepeatLib"])))


def applyFilters(ID, SEQUENCES, FINALCLASSIFICATION, CONFIG, DICOLIBRARIES):
	"""
	Apply filters from the CONFIG file to creates the final libraries.

	@type ID: string
	@param ID: string of the id of the sequence
	@type SEQUENCES: dictionnary
	@param SEQUENCES: dictionnary containing the fasta sequence of the preLibraries
	@type FINALCLASSIFICATION: string
	@param FINALCLASSIFICATION: string indicating the final classification of the sequence (usefull to compare SEQUENCES with the CONFIG file)
	@type CONFIG: dictionnary
	@param CONFIG:  dictionnary with parameters to filter sequences (which will be used to create final libraries) {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }

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
					#### Control in which libraries the sequence can be added
					#### List of ID for construct librarie of potential autonomous TE
					if CONFIG[sortedOutput]["autonomousLib"].lower()=="yes":
						DICOLIBRARIES["autonomousLib"].append(ID)
					#### List of ID for construct librarie of total TE
					if CONFIG[sortedOutput]["totalTELib"].lower()=="yes":
						DICOLIBRARIES["totalTELib"].append(ID)
					#### List of ID for construct librarie of repeated Elements
					if CONFIG[sortedOutput]["totalRepeatLib"].lower()=="yes":
						DICOLIBRARIES["totalRepeatLib"].append(ID)

					#### Save intermediate libraries
					save.saveIntermadiateLibraries(SEQUENCES, ID, CONFIG, sortedOutput)
