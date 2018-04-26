#!/usr/bin/python3
# coding: utf8

""" @author: Tristan Frances """

def save(FASTA, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""

	Save the categorized sequences in different FASTA files.

	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type NONTE: dictionnary
	@param NONTE: dictionnary for non transposable element (nonTE).
	@type POTENTIALCHIMERIC: dictionnary
	@param POTENTIALCHIMERIC: dictionnary for potential chimeric element.
	@type NOCAT: dictionnary
	@param NOCAT: dictionnary for non categorized element (noCat).
	@type TE: dictionnary
	@param TE: dictionnary for transposable element (I or II).

	@rtype: None
	"""
	#TODO: Organize the way TE are saved
	saveNoCat(FASTA, NOCAT)
	savePotentialChimeric(FASTA, POTENTIALCHIMERIC)
	saveTE(FASTA, TE)
	saveNonTE(FASTA, NONTE)

def saveNoCat(FASTA, NOCAT):
	"""

	Save the sequences considered as noCat in a FASTA file.

	Keyword arguments:
	@type FASTA: dictionnary
	@param FASTA: dictionnary with the nucleotide sequence which has been categorized
	@type FASTA: string
	@param FASTA: name of the FASTA file containing the sequence that will be opened.

	@rtype: None
	"""
	with open("noCat.fasta", "w") as f:
		for id in NOCAT.keys():
			f.write(">{}:{}\n{}\n".format(id, NOCAT[id]["class"], FASTA[id]["seq"]))


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
			f.write(">{}:{}\n{}\n".format(id, POTENTIALCHIMERIC[id]["class"], FASTA[id]["seq"]))

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
			f.write(">{}:{}:{}\n{}\n".format(id, TE[id]["class"], TE[id]["order"], FASTA[id]["seq"]))


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
			f.write(">{}:{}\n{}\n".format(id, NONTE[id]["class"], FASTA[id]["seq"]))
