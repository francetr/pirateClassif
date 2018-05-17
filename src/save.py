#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

def save(FASTA, SEQUENCESCLASSIFIED):
	"""

	Save the categorized sequences in different FASTA files.

	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type SEQUENCESCLASSIFIED: dictionnary
	@param SEQUENCESCLASSIFIED: dictionnary storing the result of the classification into 4 dictionnaries (TE, nonTE, potentialChimeric and noCat)

	@rtype: None
	"""
	#TODO: Organize the way TE are saved
	saveNoCat(FASTA, SEQUENCESCLASSIFIED["noCat"])
	savePotentialChimeric(FASTA, SEQUENCESCLASSIFIED["potentialChimeric"])
	saveTE(FASTA, SEQUENCESCLASSIFIED["TE"])
	saveNonTE(FASTA, SEQUENCESCLASSIFIED["nonTE"])

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
	with open("noCat.fasta", "w") as f:
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
	with open("potentialChimeric.fasta", "w") as f:
		for id in POTENTIALCHIMERIC.keys():
			f.write(">{id}:{seqClass}\n{seq}\n".format(id=id, seqClass=POTENTIALCHIMERIC[id]["class"], seq=FASTA[id]["seq"]))

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
	with open("TE.fasta", "w") as f:
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
	with open("nonTE.fasta", "w") as f:
		for id in NONTE.keys():
			f.write(">{id}:{seqClass}\n{seq}\n".format(id=id, seqClass=NONTE[id]["class"], seq=FASTA[id]["seq"]))
