#!/usr/bin/python2
# coding: utf8

import re
import os
import shutil
import subprocess
import readInput

""" @author: Tristan Frances """

def saveResult(FASTA, SEQCLASSIFIED):
	"""

	Save the categorized sequences in different FASTA files.
	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 8 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 4 : for the results (length, completeness, class and finalDegree), that we'll find for each sequences
		- 3 : for the order, for the predictedSuperFamily and 1 for the proofs. (These 2 dic are just for TE or potentialChimeric)
		- 1 : for the unknown keyword (unknown_keyword)

	@rtype: None
	"""
	#TODO : Save sequences according their class, order and their superFamily
	####	Creat a directory to stor each files
	try:
		#### Remove the existing prelibraries directories if the script has already been launched previously
		shutil.rmtree("classification_result")
	except OSError as e:
		# print(os.strerror(e.errno))
		pass
	finally:
		#### Creates the prelibraries directories
		os.mkdir("classification_result")
		os.mkdir("classification_result/prelibraries/")
		os.mkdir("classification_result/prelibraries/TE")
		os.mkdir("classification_result/prelibraries/TE/class")
		os.mkdir("classification_result/prelibraries/TE/order")
		os.mkdir("classification_result/prelibraries/TE/superFamily")

	#TODO construction automatic of different libraries
	#### List all the prelibraries that will be created
	allSaves = findAllSaves(SEQCLASSIFIED)

	fileSummary = open("classification_result/Classification_summary.txt", "w")
	print("Save classification summary of the sequences into \"%s\" file"%(fileSummary.name))

	fileUnknownKeyword = open("classification_result/unknown_keyword.txt", "w")
	print("Save the unknown keywords into \"%s\" file"%(fileUnknownKeyword.name))

	#### Write the headers of the Classification_summary file
	fileSummary.write("ID\tLENGTH\tSAVE_TYPE\tCOMPLETENESS\tCLASS\tORDER\tPREDICTED_SUPERFAMILY\tBLAST_PROOFS\tPROT_PROOFS\tFINAL_DEGREE\n")

	#### Create files for each superFamily and order
	for seqName in SEQCLASSIFIED:
		#### Save sequences according their saveType (nonTE, or potentialChimeric, or noCat)
		for saveType in allSaves["saveType"]:
			if (SEQCLASSIFIED[seqName]["saveType"] == saveType) and (saveType != "TE"):
				fileType=open("classification_result/prelibraries/{saveType}.fasta".format(saveType=saveType), "a")
				fileType.write(">{seqName}\n{seq}\n".format(seqName=seqName, seq=FASTA[seqName]['seq']))
				fileType.close()

		#### Only for TE sequences
		if SEQCLASSIFIED[seqName]["saveType"] == "TE":
			#### Save sequences with finalDegree equal to a class (I or II)
			for seqClass in allSaves["class"]:
				if SEQCLASSIFIED[seqName]["finalDegree"] == seqClass:
					fileClass=open("classification_result/prelibraries/TE/class/{seqClass}_undefined.fasta".format(seqClass=seqClass), "a")
					fileClass.write(">{seqName}_{seqClass}_undefined_undefined\n{seq}\n".format(seqName=seqName, seqClass=SEQCLASSIFIED[seqName]['class'], seq=FASTA[seqName]['seq']))
					fileClass.close()

			#### Save sequences with finalDegree equal to an order (LTR, DIRS, PLE, ...)
			for order in allSaves["order"]:
				if SEQCLASSIFIED[seqName]["finalDegree"] == str("%s_undefined"%(order)):
					fileOrder=open("classification_result/prelibraries/TE/order/{order}_undefined.fasta".format(order=order), "a")
					fileOrder.write(">{seqName}_{seqClass}_{order}_undefined\n{seq}\n".format(seqName=seqName, seqClass=SEQCLASSIFIED[seqName]['class'], order=SEQCLASSIFIED[seqName]['order'], seq=FASTA[seqName]['seq']))
					fileOrder.close()

			#### Check there is a superFamily for the sequence
			if "superFamily" in SEQCLASSIFIED[seqName]:
				#### Save sequences according their superFamily, if exists, (Copia, Gypsy, Mariner ...)
				for superFamily in allSaves["superFamily"]:
					#### Check the superFamily name is not the same than the order name (like Maverick or Helitron)
					if SEQCLASSIFIED[seqName]["finalDegree"] == superFamily:
						fileSuperFamily=open("classification_result/prelibraries/TE/superFamily/{superFamily}.fasta".format(superFamily=superFamily), "a")
						fileSuperFamily.write(">{seqName}_{seqClass}_{order}_{superFamily}\n{seq}\n".format(seqName=seqName, seqClass=SEQCLASSIFIED[seqName]['class'], order=SEQCLASSIFIED[seqName]['order'], superFamily=SEQCLASSIFIED[seqName]['superFamily'], seq=FASTA[seqName]['seq']))
						fileSuperFamily.close()

		#### Save the unknown keywords of the baseline, if any
		if SEQCLASSIFIED[seqName]["unknown_keyword"] != "\n" :
			print("Write unknown keywords for sequence %s in the file %s"%(seqName, fileUnknownKeyword.name))
			saveUnknownKeyword(fileUnknownKeyword, SEQCLASSIFIED[seqName])

		#### Save ckassification summary of the sequences
		saveSummary(fileSummary, seqName, SEQCLASSIFIED[seqName])
	fileUnknownKeyword.close()
	fileSummary.close()
	print("Number of different saveType files : {saveTypeFile}\nNumber of different class files : {classFile}\nNumber of different order files : {orderFile}\nNumber of different superFamily files : {superFamilyFile}\n".format(\
	saveTypeFile=len(os.listdir("classification_result/prelibraries/")), classFile=len(os.listdir("classification_result/prelibraries/TE/class")), orderFile=len(os.listdir("classification_result/prelibraries/TE/order")), superFamilyFile=len(os.listdir("classification_result/prelibraries/TE/superFamily"))))

def findAllSaves(SEQCLASSIFIED):
	"""
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 8 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 4 : for the results (length, completeness, class and finalDegree), that we'll find for each sequences
		- 3 : for the order, for the predictedSuperFamily and 1 for the proofs. (These 2 dic are just for TE or potentialChimeric)
		- 1 : for the unknown keyword (unknown_keyword)

	@rtype: dictionnary
	@return: dictionnary containing the names of the prelibraries files that will be created.
	"""
	allSaves = {"saveType":[], "class":[], "order":[], "superFamily":[]}
	for sequence in SEQCLASSIFIED.keys():
		#### Check if the saveType of the sequence is not in the saveType list, avoid redundancy
		if not SEQCLASSIFIED[sequence]["saveType"] in allSaves["saveType"]:
			allSaves["saveType"].append(SEQCLASSIFIED[sequence]["saveType"])
		if SEQCLASSIFIED[sequence]["saveType"] == "TE":
			#### Check if there is a defined superFamily for the sequence
			if not SEQCLASSIFIED[sequence]["class"] in allSaves["class"]:
				#### Check if the class of the sequence is not in the saveType list, avoid redundancy
				allSaves["class"].append(SEQCLASSIFIED[sequence]["class"])

			#### Check if there is a defined superFamily for the sequence
			if not SEQCLASSIFIED[sequence]["order"] in allSaves["order"]:
				#### Check if the order of the sequence is not in the saveType list, avoid redundancy
				allSaves["order"].append(SEQCLASSIFIED[sequence]["order"])

			#### Check if there is a defined superFamily for the sequence
			if "superFamily" in SEQCLASSIFIED[sequence]:
				#### Search for possible superFamily name
				name = re.search(r'([^_]+)', SEQCLASSIFIED[sequence]["superFamily"]).groups()[0]
				#### Check if the superFamily of the sequence is not in the saveType list, avoid redundancy
				if not name in allSaves["superFamily"]:
					#### Define the name of the superFamily sequence as FINALSUPERFAMILYNAME
					allSaves["superFamily"].append(name)
	return allSaves

def saveUnknownKeyword(FILEUNKNOWNKEYWORD, UNKNOWNKEYWORD):
	"""
	Save the unknown_keyword if keywords hasn't be found in sequences in an unknown_keyword file.

	Keyword arguments:
	@type FILEUNKNOWNKEYWORD: TextIOWrapper
	@param FILEUNKNOWNKEYWORD: File onto which the unknown keywords for a sequence, during search superFamily name, will be written
	@type UNKNOWNKEYWORD: string
	@param UNKNOWNKEYWORD: name of the FASTA file containing the sequence that will be opened.

	@rtype: None
	"""
	FILEUNKNOWNKEYWORD.write("{unknown_keyword}".format(unknown_keyword=UNKNOWNKEYWORD["unknown_keyword"]))


def saveSummary(FILESUMMARY, SEQNAME, SUMMARY):
	"""
	Save the summary of the classification steps.

	Keyword arguments:
	@type FILESUMMARY: TextIOWrapper
	@param FILESUMMARY: File onto which the summary of the sequences will be written
	@type SEQNAME: string
	@param SEQNAME: Name of the sequence
	@type SUMMARY: string
	@param SUMMARY: Proofs used during the classification step

	@rtype: None
	"""
	FILESUMMARY.write("{id}\t{length}\t{saveType}\t{completeness}\t{seqClass}".format(id=SEQNAME, length=SUMMARY["length"], saveType=SUMMARY["saveType"], completeness=SUMMARY["completeness"], seqClass=SUMMARY["class"]))
	if "order" in SUMMARY:
		FILESUMMARY.write("\t{order}".format(order=SUMMARY["order"]))
	else:
		FILESUMMARY.write("\tNA")
	if "superFamily" in SUMMARY:
		FILESUMMARY.write("\t{superFamily}\t{blastProofs}\t{protProofs}".format(superFamily=SUMMARY["superFamily"], blastProofs=SUMMARY["superFamilyProofs"]["blast"], protProofs=SUMMARY["superFamilyProofs"]["protProfiles"]))
	else:
		FILESUMMARY.write("\tNA\tNA\tNA")
	FILESUMMARY.write("\t{finalDegree}\n".format(finalDegree=SUMMARY["finalDegree"]))

def saveFilteredSequences(CONFIG):
	"""
	Save filtered sequences according the config file.

	Keyword arguments:
	@type CONFIG: dictionnary
	@param CONFIG:  dictionnary with parameters to filter sequences (which will be used to create final libraries) {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }
		There are currently 8 keys for this dictionnary:
		- outputName : name of the output file
		- lengthMin : minimal length of the sequence
		- lengthMax : maximum length of the sequence
		- removedTool : name of the tool used on the sequence. Remove the sequences found by this tool
		- onlySelectedtools : name of the tool used on the sequence. Select only the sequences found by this tool
		- autonomousLib : boolean (yes/no) indicating if sequence constitute autonomous librarie
		- totalTELib : boolean (yes/no) indicating if sequence constitute total TEs librarie
		- totalRepeatLib : boolean (yes/no) indicating if sequence constitute repeated elements librarie

	@rtype: None
	"""
	#### find all the prelibraries
	findTE = subprocess.Popen('find classification_result/prelibraries/TE -name "*.fasta" 2> /dev/null', shell=True, stdout=subprocess.PIPE);
	#### String containing all the preLibraries file name
	preLibraries = findTE.stdout.read().decode("utf-8").strip()
	
	#### find noCat prelibrary
	findNoCat = subprocess.Popen('find classification_result/prelibraries/ -name "noCat.fasta" 2> /dev/null', shell=True, stdout=subprocess.PIPE);
	#### String containing all the preLibraries file name
	noCatLibrarie = findNoCat.stdout.read().decode("utf-8").strip()

	listPrelibraries = []
	#### dictionnaries that will contains all the id's sequences for concerned libraries
	dicoLibraries={"listTEAutonomous":[], "listTotalTE":[], "listRepeated":[]}

	listPrelibraries.append(noCatLibrarie)
	#### Add all the name of prelibraries in listPrelibraries
	for file in preLibraries.split("\n"):
		listPrelibraries.append(file)
	#### Parse all the prelibrary
	for preLibrary in listPrelibraries:
		#### Retrieve the final classification name of the ET from the file name
		finalClassification = os.path.basename(preLibrary).split(".fasta")[0]
		#### Read and store the fasta sequences of the prelibraries
		sequences=readInput.readFasta(preLibrary)
		#### Parse all the sequences
		for id in sequences:
			#### Check the finalClassification of the sequences is in the ID
			if finalClassification in id:
				applyFilters(id, sequences, finalClassification, CONFIG, dicoLibraries)
				pass


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
	@type DICOLIBRARIES: dictionnary
	@param DICOLIBRARIES:  dictionnary that will contain the id of the sequences to construct the 3 final libraries {B{key} = hightest classification given for the TE : B{I{value}} = {"seq" : sequence of the CONFIG sequence} }

	@rtype: None
	"""
	#### First we check if tools used to detect the TE are from the removedTool or the onlySelectedtools of the CONFIG file
	#### Find if the tool used to find the TE is a part of the removedTool from the CONFIG file : True it is find, else False
	findRemovedTool=False
	for removedTool in CONFIG[FINALCLASSIFICATION]["removedTool"]:
		if ID.lower().find(removedTool.lower()) == 0:
			findRemovedTool=True

	#### Find if the tool used to find the TE is one of the onlySelectedtools from the CONFIG file : True it is find, else False
	findSelectedTool=False
	for selectedTool in CONFIG[FINALCLASSIFICATION]["onlySelectedtools"]:
		if ID.lower().find(selectedTool.lower()) == 0:
			findSelectedTool=True
	# print(ID, CONFIG[FINALCLASSIFICATION]["onlySelectedtools"], findSelectedTool)

	#### Control the length of the sequence with the CONFIG
	if SEQUENCES[id]["length"] >= CONFIG[FINALCLASSIFICATION]["lengthMin"] and SEQUENCES[id]["length"] <= CONFIG[FINALCLASSIFICATION]["lengthMax"]:
		pass
