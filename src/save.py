#!/usr/bin/python3
# coding: utf8

import re
import os
import shutil

""" @author: Tristan Frances """

def save(FASTA, SEQCLASSIFIED):
	"""
	Save the categorized sequences in different FASTA files.

	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 8 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 3 : for the results (length, class and finalDegree), that we'll find for each sequences
		- 3 : for the order, for the predictedSuperFamily and 1 for the proofs. (These 2 dic are just for TE or potentialChimeric)
		- 1 : for the unknown keyword (unknown_keyword)

	@rtype: None
	"""
	#TODO : Save sequences according their class, order and their superFamily
	####	Creat a directory to stor each files
	try:
		#### Remove the existing libraries directories if the script has already been launched previously
		shutil.rmtree("libraries")
	except OSError as e:
		# print(os.strerror(e.errno))
		pass
	finally:
		#### Creates the libraries directories
		os.mkdir("libraries/")
		os.mkdir("libraries/class")
		os.mkdir("libraries/order")
		os.mkdir("libraries/superFamily")

	#### List all the libraries that will be created
	#TODO construction automatic of different libraries
	allSaves = {"saveType":[], "class":[], "order":[], "superFamily":[]}
	for sequence in SEQCLASSIFIED.keys():
		allSaves["saveType"].append(SEQCLASSIFIED[sequence]["saveType"])
		if SEQCLASSIFIED[sequence]["saveType"] == "TE":
			allSaves["class"].append(SEQCLASSIFIED[sequence]["class"])
			allSaves["order"].append(SEQCLASSIFIED[sequence]["order"])
			if "superFamily" in SEQCLASSIFIED[sequence]:
				#### Search for possible superFamily name
				name = re.search(r'([^_]+)', SEQCLASSIFIED[sequence]["superFamily"]).groups()[0]
				#### Define the name of the superFamily sequence as FINALSUPERFAMILYNAME
				allSaves["superFamily"].append(name)

	allSaves["saveType"] = list(set(allSaves["saveType"]))
	allSaves["class"] = list(set(allSaves["class"]))
	allSaves["order"] = list(set(allSaves["order"]))
	allSaves["superFamily"] = list(set(allSaves["superFamily"]))

	fileSummary = open("Classification_summary.txt", "w")
	print("Save classification summary of the sequences into \"%s\" file"%(fileSummary.name))

	fileUnknownKeyword = open("unknown_keyword.txt", "w")
	print("Save the unknown keywords into \"%s\" file"%(fileUnknownKeyword.name))

	#### Create files for each superFamily and order
	for seqName in SEQCLASSIFIED:
		#### Save sequences according their saveType (TE or nonTE, or potentialChimeric, or noCat)
		for saveType in allSaves["saveType"]:
			if SEQCLASSIFIED[seqName]["saveType"] == saveType:
				fileType=open("libraries/{saveType}.fasta".format(saveType=saveType), "a")
				fileType.write(">{seqName}\n{seq}\n".format(seqName=seqName, seq=FASTA[seqName]['seq']))
				fileType.close()

		#### Save sequences according their class (I or II)
		for seqClass in allSaves["class"]:
			if SEQCLASSIFIED[seqName]["class"] == seqClass:
				fileClass=open("libraries/class/{seqClass}.fasta".format(seqClass=seqClass), "a")
				fileClass.write(">{seqName}\n{seq}\n".format(seqName=seqName, seq=FASTA[seqName]['seq']))
				fileClass.close()

		#### Save sequences according their order, if exist (LTR, DIRS, PLE, ...)
		if "order" in SEQCLASSIFIED[seqName]:
			for order in allSaves["order"]:
				if SEQCLASSIFIED[seqName]["order"] == order:
					fileOrder=open("libraries/order/{order}.fasta".format(order=order), "a")
					fileOrder.write(">{seqName}\n{seq}\n".format(seqName=seqName, seq=FASTA[seqName]['seq']))
					fileOrder.close()

		#### Save sequences according their superFamily, if exists, (Copia, Gypsy, Mariner ...)
		if "superFamily" in SEQCLASSIFIED[seqName]:
			for superFamily in allSaves["superFamily"]:
				if SEQCLASSIFIED[seqName]["superFamily"] == superFamily:
					fileSuperFamily=open("libraries/superFamily/{superFamily}.fasta".format(superFamily=superFamily), "a")
					fileSuperFamily.write(">{seqName}\n{seq}\n".format(seqName=seqName, seq=FASTA[seqName]['seq']))
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
	saveTypeFile=len(allSaves["saveType"]), classFile=len(allSaves["class"]), orderFile=len(allSaves["order"]), superFamilyFile=len(allSaves["superFamily"])))

# def writeInFile(FILE, FASTA, SEQNAME, SAVETYPE):
# 	####	Case saveType is uncategorized sequences
# 	if SAVETYPE == "noCat":
# 		FILE.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=SAVETYPE["class"], seq=FASTA[SEQNAME]["seq"]))
#
# 	####	Case saveType is potential Chimeric sequences
# 	elif SEQCLASSIFIED[seqName]["saveType"] == "potentialChimeric":
# 		####	Case there is not superFamily determined for this sequence
# 		if not "superFamily" in FILE:
# 			FILE.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=SAVETYPE["class"], seq=FASTA[SEQNAME]["seq"]))
# 		####	Case a superFamily have been determined for this sequence
# 		else:
# 			FILE.write(">{seqName}:{seqClass}:{order}:{superFamily}\n{seq}\n".format(seqName=SEQNAME, seqClass=SAVETYPE["class"], order=SAVETYPE["order"], superFamily=SAVETYPE["superFamily"], seq=FASTA[SEQNAME]["seq"]))
#
# 	####	Case saveType is TE
# 	elif SEQCLASSIFIED[seqName]["saveType"] == "TE":
# 		####	Case there is not superFamily determined for this sequence
# 		if not "superFamily" in FILE:
# 			FILE.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=NOCAT["class"], seq=FASTA[SEQNAME]["seq"]))
# 		####	Case a superFamily have been determined for this sequence
# 		else:
# 			FILE.write(">{seqName}:{seqClass}:{order}:{superFamily}\n{seq}\n".format(seqName=SEQNAME, seqClass=SAVETYPE["class"], order=SAVETYPE["order"], superFamily=SAVETYPE["superFamily"], seq=FASTA[SEQNAME]["seq"]))
#
# 	####	Case saveType is nonTE
# 	elif SEQCLASSIFIED[seqName]["saveType"] == "nonTE":
# 		FILE.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=SAVETYPE["class"], seq=FASTA[SEQNAME]["seq"]))
#
# 	####	Save the unknown keywords (keywords unknown, founded during search process in a sequence), just if unknown keywords have been found
# 	if SEQCLASSIFIED[seqName]["unknown_keyword"] != "\n" :
# 		print("Write unknown keywords for sequence %s in the file %s"%(seqName, fileUnknownKeyword.name))
# 		saveUnknownKeyword(fileUnknownKeyword, SEQCLASSIFIED[seqName])
# 	####	Save summary of the sequences
# 	saveSummary(fileLog, SEQCLASSIFIED[seqName])


#### ONLY THE FUNCTIONS OF SUMMARY AND UNKNOWNKEYWORD ARE USED RIGHT NOW
# def saveNoCat(FILENOCAT, FASTA, SEQNAME, NOCAT):
# 	"""
# 	Save the sequences considered as noCat in a FASTA file.
#
# 	Keyword arguments:
# 	@type FILENOCAT: TextIOWrapper
# 	@param FILENOCAT: File onto which the uncategorized sequences will be written
# 	@type FASTA: dictionnary
# 	@param FASTA: dictionnary with the nucleotide sequence which has been classified
# 	@type SEQNAME: string
# 	@param SEQNAME: name of the saved sequence
# 	@type NOCAT: dictionnary
# 	@param NOCAT: dictionnary for uncategorized sequences.
#
# 	@rtype: None
# 	"""
# 	FILENOCAT.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=NOCAT["class"], seq=FASTA[SEQNAME]["seq"]))
#
#
# def savePotentialChimeric(FILECHIMERIC, FASTA, SEQNAME, POTENTIALCHIMERIC):
# 	"""
# 	Save the sequences considered as potentialChimeric in a FASTA file.
#
# 	Keyword arguments:
# 	@type FILECHIMERIC: TextIOWrapper
# 	@param FILECHIMERIC: File onto which the potential chimeric sequences will be written
# 	@type FASTA: dictionnary
# 	@param FASTA: dictionnary with the nucleotide sequence which has been classified
# 	@type SEQNAME: string
# 	@param SEQNAME: name of the saved sequence
# 	@type POTENTIALCHIMERIC: dictionnary
# 	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
#
# 	@rtype: None
# 	"""
# 	####	Case chimeric have been found first
# 	if not "superFamily" in POTENTIALCHIMERIC:
# 		FILECHIMERIC.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=POTENTIALCHIMERIC["class"], seq=FASTA[SEQNAME]["seq"]))
# 	####	Case chimeric have been found during superFamily determination
# 	else:
# 		FILECHIMERIC.write(">{seqName}:{seqClass}:{order}:{superFamily}\n{seq}\n".format(seqName=SEQNAME, seqClass=POTENTIALCHIMERIC["class"], order=POTENTIALCHIMERIC["order"], superFamily=POTENTIALCHIMERIC["superFamily"], seq=FASTA[SEQNAME]["seq"]))
#
#
# def saveTE(FILETE, FASTA, SEQNAME, TE):
# 	"""
# 	Save the sequences considered as TE in a FASTA file.
#
# 	Keyword arguments:
# 	@type FILETE: TextIOWrapper
# 	@param FILETE: File onto which the TE sequences will be written
# 	@type FASTA: dictionnary
# 	@param FASTA: dictionnary with the nucleotide sequence which has been classified
# 	@type SEQNAME: string
# 	@param SEQNAME: name of the saved sequence
# 	@type TE: dictionnary
# 	@param TE: dictionnary for Transposable Element.
#
# 	@rtype: None
# 	"""
# 	if not "superFamily" in TE:
# 		FILETE.write(">{seqName}:{seqClass}:{order}\n{seq}\n".format(seqName=SEQNAME, seqClass=TE["class"], order=TE["order"], seq=FASTA[SEQNAME]["seq"]))
# 	else:
# 		FILETE.write(">{seqName}:{seqClass}:{order}:{superFamily}\n{seq}\n".format(seqName=SEQNAME, seqClass=TE["class"], order=TE["order"], superFamily=TE["superFamily"], seq=FASTA[SEQNAME]["seq"]))
#
# def saveNonTE(FILENONTE, FASTA, SEQNAME, NONTE):
# 	"""
# 	Save the sequences considered as NONTE in a FASTA file.
#
# 	Keyword arguments:
# 	@type FILENONTE: TextIOWrapper
# 	@param FILENONTE: File onto which the nonTE sequences will be written
# 	@type FASTA: dictionnary
# 	@param FASTA: dictionnary with the nucleotide sequence which has been classified
# 	@type SEQNAME: string
# 	@param SEQNAME: name of the saved sequence
# 	@type NONTE: dictionnary
# 	@param NONTE: dictionnary for non Transposable Element.
#
# 	@rtype: None
# 	"""
# 	FILENONTE.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=NONTE["class"], seq=FASTA[SEQNAME]["seq"]))

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
	FILESUMMARY.write("{id}\t{saveType}\tCLASS:\t{seqClass}\tFINALDEGREE:\t{finalDegree}".format(id=SEQNAME, saveType=SUMMARY["saveType"], seqClass=SUMMARY["class"], finalDegree=SUMMARY["finalDegree"]))
	if "order" in SUMMARY:
		FILESUMMARY.write("\tORDER:\t{order}".format(order=SUMMARY["order"]))
	if "superFamily" in SUMMARY:
		FILESUMMARY.write("\tPREDICTED_SUPERFAMILY:\t{superFamily}\t{proofs}".format(superFamily=SUMMARY["superFamily"], proofs=SUMMARY["superFamilyProofs"]))
	FILESUMMARY.write("\n")
