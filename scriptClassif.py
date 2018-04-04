#!/usr/bin/python3.5
# coding: utf8
import sys
from Bio import SeqIO # For the fasta reading
import argparse
import re



############
#   Commandes to launch this script :
#	   python3 <path/toward/this/script> <path/toward/the/classif/file> <path/toward/the/fasta/file>
# example	   python3 code/scriptClassif.py sortie_pastec/TisoRepet1.classif hat.fasta
#	   or  ./<path/toward/the>/scriptClassif.py <path/toward/the/classif/file> <path/toward/the/fasta/file>
############

# @author: Tristan Frances

def retrieveArguments():
	"""
	Retrieve the arguments from the command line

	:return args: attributes NameSpace of object argParse that contain the two arguments name

	"""
	####	 Mananage the 2 arguments (PASTEC and FASTA file name) when the command is launched
	parser = argparse.ArgumentParser(prog="scriptClassif.py", description="This program is a part of the PiRATE project. It aims\
	to automatized the step of TE classification")
	parser.add_argument("classif", type=str, help="classif file that comes from PASTEC")
	parser.add_argument("fasta", type=str, help="fasta file providing the sequence")
	args = parser.parse_args()
	# checkArguments(args.classif, args.fasta)

	return(args)

# def checkArguments(CLASSIFNAME, FASTANAME):
# 	"""
# 	Check the extension of the two arguments of command line. If it is not correct, stop the programm.
#
# 	Keyword argument
# 	:param CLASSIFNAME: name of the argument for the CLASSIF file
# 	:param FASTANAME: name of the argument for the FASTA file
#
# 	:return: None
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



def readPastec(PASTEC, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""
	Read a PASTEC file as input line by line, and then proceed to the categorization of each sequence

	Keyword argument
	:param PASTEC: name of the classif file that will be opened
	:param NONTE: dictionnary for non transposable element (nonTE)
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:parma TE: dictionnary for transposable element (I or II)

	:return: None

	"""
	try:
		####	 open the classif file
		with open(PASTEC, "r") as f:
		####	 parse every line/sequence of the file
			for line in f:
				sequence=line.replace("\t\t", "\t").strip() # remove the first column double tab and carriage return of the line
				categorization(sequence, NONTE, POTENTIALCHIMERIC, NOCAT, TE)
				# print("Summary, sequences find by type :\nTE : {}\nnonTE : {}\nnoCat : {}\npotentialChimeric : {}\nTotal : {}".format(len(TE),\
				# len(nonTE), len(noCat), len(potentialChimeric), (len(TE)+len(nonTE)+len(noCat)+len(potentialChimeric))))

	except FileNotFoundError:
		####	 prevent the opening if the file name is incorrect
		print("####	No such file as {}\n####	Classification aborted".format(PASTEC))
		sys.exit(1)

def readFasta(FASTA):
	"""
	Read a FASTA file as input and return it as a string. Use the SeqIO from the package Bio of BioPython

	Return a dictionnary with keys:
	id : id of the sequence
	name : name of the sequence
	description : description of the sequence
	number of features : number of features of the sequence
	seq : sequence concerned

	Keyword argument
	:param FASTA: name of the FASTA file containing the sequence that will be opened

	:return seqReturned: Dictionnary with key: id of the sequence and value : {"seq" : sequence of the FASTA sequence}

	"""
	seqReturned={}
	try:
		####	 Open the fasta file
		with open(fastaFile, "rU") as handle:
		####	 parse every line/sequence of the file
			for record in SeqIO.parse(handle, "fasta"):
				####	 save the sequence in a dictionnary
				seqReturned[record.id]={"seq":record.seq}
		return(seqReturned)
	except(FileNotFoundError, NameError):
		####	prevent the opening if the file name is incorrect
		print("####	No such file as {}\n####	Classification aborted".format(FASTA))
		sys.exit(1)

def categorization(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""
	Categorize the Transposable Element from the input argument

	Return a dictionnary which contains the different catagories which caracterize the sequence.

	Keyword argument
	:param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed. Usefull values of this list :
	- [0] : id of the sequence
	- [3] : potentiel chimeric
	- [4] : class of the sequence
	- [5] : order of the sequence
	- [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
	:param NONTE: dictionnary for non transposable element (nonTE)
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	:return: None

	"""
	# listOrder=["SSR", "noCat", "LTR", "TRIM", "PotentialHost", "TIR", "DIRS", "LARD", "LINE", "MITE", "Helitron", "Maverick", "SINE"]
	features=SEQUENCE.split("\t")
	#TODO Treatment of the superfamily
	classDetermination(features, nonTE, potentialChimeric,  noCat, TE)


def classDetermination(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	"""
	Determine the class of a sequence.
	For this, complete four dictionnaries, passed onto arguments, that will contain the different catagories that caracterize the sequence.

	Keyword argument
	:param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed
	:param NONTE: dictionnary for non transposable element (nonTE)
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	:return: None

	"""
	# typeOrder=None
	####	 check first if the sequence is chimeric
	if(SEQUENCE[3] != "PotentialChimeric"):
		####	 check the class of the sequence if it is not chimeric : class I
		if(SEQUENCE[4] == "I"):
			TE[SEQUENCE[0]]={"class":"I"}
			# print(SEQUENCE[0])
			orderDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE)
			pass
		####	 classII
		elif(SEQUENCE[4] == "II"):
			TE[SEQUENCE[0]]={"class":"II"}
			orderDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE)
			# print(SEQUENCE[0])
			pass
		####	 NoCat
		elif(SEQUENCE[4] == "noCat"):
			####	 noCat
			NOCAT[SEQUENCE[0]]={"class":"noCat"}
			# print("NoCat : %s %s" % (SEQUENCE[4], SEQUENCE[5]))
			pass
		####	 NonTE
		else:
			NONTE[SEQUENCE[0]]={"class":"nonTE"}
			# print("NA : %s %s" % (SEQUENCE[4], SEQUENCE[5]))
			pass
	else:
		POTENTIALCHIMERIC[SEQUENCE[0]]={"Class":"potentialChimeric"}


def orderDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE):
	"""
	Determine the order of a sequence. Doing so it complete a dictionnary containing
	the transposable element. If the order wasn't found, the sequence is placed in
	NOCAT dictionnary and removed from TE dictionnary

	Keyword argument
	:param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	:return: None

	"""

	TE[SEQUENCE[0]]["order"]=SEQUENCE[5]
	####	 The order of the TE is determined : the super family of the TE will be determined
	if (SEQUENCE[5]!="noCat"):
		superFamilyDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE)

def superFamilyDetermination(SEQUENCE, POTENTIALCHIMERIC, NOCAT, TE):
	#TODO if no coding part, move to POTENTIALCHIMERIC or NOCAT dictionnary?
	"""
	Determine the super family of a sequence. Doing so it complete a dictionnary
	containing the transposable element.
	If the superFamily can't be defined, it will be noCat
	If there is 2 or more superfamily possible, the sequence will be added to the POTENTIALCHIMERIC dictionnary and then removed from TE dictionnary

	:param SEQUENCE: name of the list of strings, which contain the sequence to categorize, that will be parsed
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	:return: None

	"""
	try:
		####	 search if there is a 'coding' part in the 7th value of SEQUENCE => needed to defined the SEQUENCE superfamily
		codingRecord = re.search(r'coding=\(([^\)]+)\)', SEQUENCE[7]).groups()[0]
		# print(codingRecord)
		####	 split the coding part to obtain the different results according to the comparison with different databases (3 possibilties:
		####	 TEBLRtx, TEBLRx and profiles)
		databaseRecords = codingRecord.split(';')
		# matches = {}
		####	 scan the different results obtain for each comparisons
		for dr in databaseRecords:
			####	 split to obtain the name of the database used
			dbName = dr.split(':')[0].strip()
			####	 for the profil hmm database, do nothing
			if dbName=='profiles': continue
			# matches[dbName] = []
			####	 split to obtain each categories (such as class, order, superfamily), and parse them in order to search interesting regex
			# TODO compare the superFamily find for each substr
			for substr in dr.split(','):
				try:
					print(substr)
					searchObj = re.search(r' ([^:]+):Class(I+):([^:]+):([^:]+): ([0-9\.]+)%', substr)
					# matches[dbName].append(searchObj.groups())
					TE[SEQUENCE[0]]["superFamily"]=searchObj.groups()[3]
				except AttributeError:
					print('Issue on : '+substr)
	except AttributeError:
		####	 there is no coding part : the sequence is put in the POTENTIALCHIMERIC dictionnary and removed from th TE dictionnary
		# print(SEQUENCE[7])
		NOCAT[SEQUENCE[0]]=TE[SEQUENCE[0]]
		del TE[SEQUENCE[0]]
		return


def save(FASTA, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
	#TODO
	"""

	Save the sequence in a file

	Return a dictionnary which contains the different catagories which caracterize the sequence.

	Keyword argument
	:param FASTA: name of the list of strings, which contain the nucleotid of the sequences
	:param NONTE: dictionnary for non transposable element (nonTE)
	:param POTENTIALCHIMERIC: dictionnary for potential chimeric element
	:param NOCAT: dictionnary for non categorized element (noCat)
	:param TE: dictionnary for transposable element (I or II)

	"""

	# for sequenceName in TE:
	#	 if FASTA.id!=sequenceName:
	#	 print(FASTA.id, sequenceName)
	# return




####################
#	   MAIN
####################

if __name__ == "__main__":
	# nom="/home/tfrances/Bureau/donnees/sortie_PASTEC/TisoMgescan.classif"

	print("Start of the classification\n")
	####	 Instanciation of dictionnaries that will contain the results
	nonTE, potentialChimeric,  noCat, TE = {}, {}, {}, {}

	####	 retrieve and check the arguments used in command line
	print("####	Check the extensions of the files")
	args = retrieveArguments();
	print("####	Valid extension of the files\n")
	pastecFile=args.classif
	####	Reading of the classif file ####
	try:
		# take the second argument in the terminal command as file name
		print("####	Read of the classif file")
		pastec=readPastec(pastecFile, nonTE, potentialChimeric, noCat, TE)
		# print(TE)
		print("####	End of reading of the classif file\n")
	except IndexError:
		print("####	No PASTEC provided\n####	Classification aborted")
		sys.exit(1)


	print("TE : %d, noCat : %d, nonTE : %d, chimeric: %.2f" %(len(TE), len(noCat), len(nonTE), len(potentialChimeric)))
	####	Procede to the categorisation
	# cat=categorization(pastec)

	# for key in TE.keys():
	#	 if (TE[key]["class"]!="PotentialChimeric"):
	#	 print(TE[key])
	#	 pass
	#
	####	Reading of the fasta file ####
	# try:
	#	 print("####	Read of the FASTA file")
	#	 fastaFile=sys.argv[2] # take the second argument in the terminal command as file name
	#	 fasta=readFasta(fastaFile)
	#	 print("####	End of reading of the FASTA file\n")
	# except IndexError:
	#	 print("####	No Fasta provided\n####	Classification aborted")
	#	 sys.exit(1)
	# for key in fasta.keys():
	#	 print(">%s :\n%s\n" %(key, fasta[key]["seq"]))
	# save(fasta, nonTE, potentialChimeric, noCat, TE)
