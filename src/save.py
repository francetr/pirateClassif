#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

def save(FASTA, SEQCLASSIFIED):
	"""
	Save the categorized sequences in different FASTA files.

	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type SEQCLASSIFIED: dictionnary
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 6 dictionnaries :
		- 1 : for the file which saves sequences (saveType: TE; or nonTE; or potentialChimeric; or noCat);
		- 3 : for the results (class, order and superFamily)
		- 1 : for the log (log)
		- 1 : for the error(error)

	@rtype: None
	"""
	#TODO: Organize the way TE are saved
	####	First create the files in which the sequences and log will be written
	fileNoCat = open("noCat.fasta", "w")
	print("Save uncategorized sequences into \"%s\" file"%(fileNoCat.name))

	filePotentialChimeric=open("potentialChimeric.fasta", "w")
	print("Save potentialChimeric sequences into \"%s\" file"%(filePotentialChimeric.name))

	fileTE = open("TE.fasta", "w")
	print("Save TE sequences into \"%s\" file"%(fileTE.name))

	fileNonTE = open("nonTE.fasta", "w")
	print("Save nonTE sequences into \"%s\" file"%(fileNonTE.name))

	fileLog = open("log.txt", "w")
	print("Save log of the sequences into \"%s\" file"%(fileLog.name))

	fileError = open("error.txt", "w")
	print("Save the errors into \"%s\" file"%(fileError.name))

	for seqName in SEQCLASSIFIED:
		####	Save uncategorized sequences
		if SEQCLASSIFIED[seqName]["saveType"] == "noCat":
			saveNoCat(fileNoCat, FASTA, seqName, SEQCLASSIFIED[seqName])
		####	Save potentialChimeric sequences
		elif SEQCLASSIFIED[seqName]["saveType"] == "potentialChimeric":
			savePotentialChimeric(filePotentialChimeric, FASTA, seqName, SEQCLASSIFIED[seqName])
		####	Save TE sequences
		elif SEQCLASSIFIED[seqName]["saveType"] == "TE":
			saveTE(fileTE, FASTA, seqName, SEQCLASSIFIED[seqName])
		####	Save nonTE sequences
		elif SEQCLASSIFIED[seqName]["saveType"] == "nonTE":
			saveNonTE(fileNonTE, FASTA, seqName, SEQCLASSIFIED[seqName])

		####	Save the errors (keywords unknown in a sequence), just if errors have been found
		if SEQCLASSIFIED[seqName]["error"] != "\n":
			print("Write unknown keywords for sequence %s in the file %s"%(seqName, fileError.name))
			saveError(fileError, SEQCLASSIFIED[seqName])
		####	Save log of the sequences
		saveLog(fileLog, SEQCLASSIFIED[seqName])

	####	Close all the saving files
	fileError.close()
	fileLog.close()
	fileNonTE.close()
	fileTE.close()
	filePotentialChimeric.close()
	fileNoCat.close()


def saveNoCat(FILENOCAT, FASTA, SEQNAME, NOCAT):
	"""
	Save the sequences considered as noCat in a FASTA file.

	Keyword arguments:
	@type FILENOCAT: TextIOWrapper
	@param FILENOCAT: File onto which the uncategorized sequences will be written
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been classified
	@type SEQNAME: string
	@param SEQNAME: name of the saved sequence
	@type NOCAT: dictionnary
	@param NOCAT: dictionnary for uncategorized sequences.

	@rtype: None
	"""
	FILENOCAT.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=NOCAT["class"], seq=FASTA[SEQNAME]["seq"]))


def savePotentialChimeric(FILECHIMERIC, FASTA, SEQNAME, POTENTIALCHIMERIC):
	"""
	Save the sequences considered as potentialChimeric in a FASTA file.

	Keyword arguments:
	@type FILECHIMERIC: TextIOWrapper
	@param FILECHIMERIC: File onto which the potential chimeric sequences will be written
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been classified
	@type SEQNAME: string
	@param SEQNAME: name of the saved sequence
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.

	@rtype: None
	"""
	if not "superFamily" in POTENTIALCHIMERIC:
		FILECHIMERIC.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=POTENTIALCHIMERIC["class"], seq=FASTA[SEQNAME]["seq"]))
	else:
		FILECHIMERIC.write(">{seqName}:{seqClass}:{order}:{superFamily}\n{seq}\n".format(seqName=SEQNAME, seqClass=POTENTIALCHIMERIC["class"], order=POTENTIALCHIMERIC["order"], superFamily=POTENTIALCHIMERIC["superFamily"], seq=FASTA[SEQNAME]["seq"]))


def saveTE(FILETE, FASTA, SEQNAME, TE):
	"""
	Save the sequences considered as TE in a FASTA file.

	Keyword arguments:
	@type FILETE: TextIOWrapper
	@param FILETE: File onto which the TE sequences will be written
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been classified
	@type SEQNAME: string
	@param SEQNAME: name of the saved sequence
	@type TE: dictionnary
	@param TE: dictionnary for Transposable Element.

	@rtype: None
	"""
	if not "superFamily" in TE:
		FILETE.write(">{seqName}:{seqClass}:{order}\n{seq}\n".format(seqName=SEQNAME, seqClass=TE["class"], order=TE["order"], seq=FASTA[SEQNAME]["seq"]))
	else:
		FILETE.write(">{seqName}:{seqClass}:{order}:{superFamily}\n{seq}\n".format(seqName=SEQNAME, seqClass=TE["class"], order=TE["order"], superFamily=TE["superFamily"], seq=FASTA[SEQNAME]["seq"]))

def saveNonTE(FILENONTE, FASTA, SEQNAME, NONTE):
	"""
	Save the sequences considered as NONTE in a FASTA file.

	Keyword arguments:
	@type FILENONTE: TextIOWrapper
	@param FILENONTE: File onto which the nonTE sequences will be written
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been classified
	@type SEQNAME: string
	@param SEQNAME: name of the saved sequence
	@type NONTE: dictionnary
	@param NONTE: dictionnary for non Transposable Element.

	@rtype: None
	"""
	FILENONTE.write(">{seqName}:{seqClass}\n{seq}\n".format(seqName=SEQNAME, seqClass=NONTE["class"], seq=FASTA[SEQNAME]["seq"]))

def saveLog(FILELOG, LOG):
	"""
	Save the log of the classification steps.

	Keyword arguments:
	@type FILELOG: TextIOWrapper
	@param FILELOG: File onto which the log of the sequences will be written
	@type LOG: string
	@param LOG: Proofs used during the classification step

	@rtype: None
	"""
	FILELOG.write(LOG["log"])


def saveError(FILEERROR, ERROR):
	"""
	Save the error if keywords hasn't be found in sequences in an error file.

	Keyword arguments:
	@type FILEERROR: TextIOWrapper
	@param FILEERROR: File onto which the unknown keywords for a sequence, during search superFamily name, will be written
	@type ERROR: string
	@param ERROR: name of the FASTA file containing the sequence that will be opened.

	@rtype: None
	"""
	FILEERROR.write("{error}".format(error=ERROR["error"]))
