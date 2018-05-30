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
	@param SEQCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)

	@rtype: None
	"""
	#TODO: Organize the way TE are saved
	saveNoCat(FASTA, SEQCLASSIFIED["noCat"])
	savePotentialChimeric(FASTA, SEQCLASSIFIED["potentialChimeric"])
	saveTE(FASTA, SEQCLASSIFIED["TE"])
	saveNonTE(FASTA, SEQCLASSIFIED["nonTE"])
	saveLog(SEQCLASSIFIED["log"])

def saveNoCat(FASTA, NOCAT):
	"""

	Save the sequences considered as noCat in a FASTA file.

	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type NOCAT: string
	@param NOCAT: name of the FASTA file containing the sequence that will be opened.

	@rtype: None
	"""
	saveFileName = "noCat.fasta"
	print("Save uncategorized sequences into \"%s\" file"%(saveFileName))
	with open(saveFileName, "w") as f:
		for id in NOCAT.keys():
			f.write(">{id}:{seqClass}\n{seq}\n".format(id=id, seqClass=NOCAT[id]["class"], seq=FASTA[id]["seq"]))


def savePotentialChimeric(FASTA, POTENTIALCHIMERIC):
	"""

	Save the sequences considered as potentialChimeric in a FASTA file.

	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.

	@rtype: None
	"""
	saveFileName = "potentialChimeric.fasta"
	print("Save potential chimeric sequences into \"%s\" file"%(saveFileName))
	with open(saveFileName, "w") as f:
		for id in POTENTIALCHIMERIC.keys():
			if not "superFamily" in POTENTIALCHIMERIC[id]:
				f.write(">{id}:{seqClass}\n{seq}\n".format(id=id, seqClass=POTENTIALCHIMERIC[id]["class"], seq=FASTA[id]["seq"]))
			else:
				f.write(">{id}:{seqClass}:{order}:{superFamily}\n{seq}\n".format(id=id, seqClass=POTENTIALCHIMERIC[id]["class"], order=POTENTIALCHIMERIC[id]["order"], superFamily=POTENTIALCHIMERIC[id]["superFamily"], seq=FASTA[id]["seq"]))


def saveTE(FASTA, TE):
	"""

	Save the sequences considered as TE in a FASTA file.

	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type TE: dictionnary
	@param TE: dictionnary for Transposable Element.

	@rtype: None
	"""
	saveFileName = "TE.fasta"
	print("Save Transposable Elements sequences into \"%s\" file"%(saveFileName))
	with open(saveFileName, "w") as f:
		for id in TE.keys():
			if not "superFamily" in TE[id]:
				f.write(">{id}:{seqClass}:{order}\n{seq}\n".format(id=id, seqClass=TE[id]["class"], order=TE[id]["order"], seq=FASTA[id]["seq"]))
			else:
				f.write(">{id}:{seqClass}:{order}:{superFamily}\n{seq}\n".format(id=id, seqClass=TE[id]["class"], order=TE[id]["order"], superFamily=TE[id]["superFamily"], seq=FASTA[id]["seq"]))

def saveNonTE(FASTA, NONTE):
	"""

	Save the sequences considered as NONTE in a FASTA file.

	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type NONTE: dictionnary
	@param NONTE: dictionnary for non Transposable Element.

	@rtype: None
	"""
	saveFileName = "nonTE.fasta"
	print("Save non Transposable Elements sequences into \"%s\" file"%(saveFileName))
	with open(saveFileName, "w") as f:
		for id in NONTE.keys():
			f.write(">{id}:{seqClass}\n{seq}\n".format(id=id, seqClass=NONTE[id]["class"], seq=FASTA[id]["seq"]))

def saveLog(LOG):
	"""

	Save the log of the classification steps.

	Keyword arguments:
	@type LOG: string
	@param LOG: Proofs used during the classification step

	@rtype: None
	"""
	saveFileName = "log.txt"
	print("Save log of proofs used for classification of the sequences into \"%s\" file"%(saveFileName))
	with open(saveFileName, "w") as f:
		f.write("Sequence name\tClass Found\tOrder found\tProofs for superFamily classification\tSuperFamily found\n")
		for sequence in LOG:
			f.write(LOG[sequence])
