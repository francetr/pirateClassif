#!/usr/bin/python3.5
# coding: utf8
import sys
from Bio import SeqIO # For the fasta reading
import argparse
import re

############
#   Commandes to launch this script :
#                   python3 <path/toward/this/script> <path/toward/the/classif/file> <path/toward/the/fasta/file>
# example           python3 code/scriptClassif.py sortie_pastec/TisoRepet1.classif hat.fasta
#               or  ./<path/toward/the>/scriptClassif.py <path/toward/the/classif/file> <path/toward/the/fasta/file>
############

# @author: Tristan Frances
def readPastec(PASTEC):
    """
    Read a PASTEC file as input and return it as a string

    Keyword argument
    PASTEC -- name of the classif file that will be opened

    """
    with open(PASTEC, "r") as f:
        string=f.read().replace("\t\t", "\t").strip() # in the order : read the file; replace the first separator; delete the blank at the end of file
        return(string)

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
    FASTA -- name of the FASTA file containing the sequence that will be opened

    """
    with open(fastaFile, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return(record)


def categorisation(PASTEC):
    """
    Categorize the Transposable Element from the input argument

    Return a dictionnary which contains the different catagories which caracterize the sequence.

    Keyword argument
    PASTEC -- name of the string that will be parsed

    """
    # NB : there seems to be \r\n character due to windows but it don't change the processing of the string if we split with \n or \r\n
    try :
        PASTEC = PASTEC.split("\n") # split the original file into a list with each lines
    except :
        pass
    nonTE, potentialChimeric,  noCat, TE = {}, {}, {}, {}
    listOrder=["SSR", "noCat", "LTR", "TRIM", "PotentialHost", "TIR", "DIRS", "LARD", "LINE", "MITE", "Helitron", "Maverick", "SINE"]
    for line in PASTEC:
        # This loop allow to split each lines of original file in order to retrieve each column and put it into a list
        tmp=line.split("\t")
        #TODO Treatment of the superfamily
        # print(tmp[7])
        classDetermination(tmp, nonTE, potentialChimeric,  noCat, TE)
        # TE[tmp[0]]={"class":tmp[4],"order":tmp[5], "other":tmp[7], "superfamily":"TODO", "autonomous":"TODO", "sequence":"TODO"}
    # print(len(PASTEC))

    print("TE : %d, noCat : %d, nonTE : %d" %(len(TE), len(noCat), len(nonTE)))

    return(TE)


def classDetermination(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
    """
    Categorize the Transposable Element from the input argument.
    For this, complete four dictionnaries, passed onto arguments, that will contain the different catagories that caracterize the sequence.

    Keyword argument
    SEQUENCE -- list of features for one sequence. Usefull values of this list :
        - [0] : id of the sequence
        - [3] : potentiel chimeric
        - [4] : class of the sequence
        - [5] : order of the sequence
        - [7] : superfamily (if determined) of the sequence. Must be extracted with Regex
    NONTE -- dictionnary for non transposable element (nonTE)
    POTENTIALCHIMERIC -- dictionnary for potential chimeric elemenet
    NOCAT -- dictionnary for non categorized elemenet (noCat)
    TE -- dictionnary for transposable elemenet (I or II)
    """
    typeOrder=None
    if(SEQUENCE[4] == "I"):
        # print("I")
        TE[SEQUENCE[0]]={"class":"I"}
        # print(SEQUENCE[0])
        pass
    elif(SEQUENCE[4] == "II"):
        TE[SEQUENCE[0]]={"class":"II"}
        # print(SEQUENCE[0])
        pass
    elif(SEQUENCE[4] == "noCat"):
        NOCAT[SEQUENCE[0]]={"class":"noCat"}
        # print("NoCat : %s %s" % (SEQUENCE[4], SEQUENCE[5]))
        pass
    else:
        NONTE[SEQUENCE[0]]={"class":"nonTE"}
        # print("NA : %s %s" % (SEQUENCE[4], SEQUENCE[5]))
        pass
    if(SEQUENCE[3] == "PotentialChimeric"):
        # print("pot")
        pass



def save(DIC):
    """
    TODO

    Save the sequence in a file

    Return a dictionnary which contains the different catagories which caracterize the sequence.

    Keyword argument
    DIC -- dictionnary containing the different sequences categorisation

    """
    return


"""
parser = argparse.ArgumentParser()
parser.add_argument("classif", type=str, help="classif file providing from PASTEC")
parser.add_argument("fasta", type=str, help="fasta file providing the sequence")

args = parser.parse_args()
answer = args.classif
try:
    classifName=re.match(r'[\S]*[/]?[\w.]+(classif)$', args.classif).groups()[0]
    fastaName=re.match(r'[\S]*[/]?[\w.]+(fasta)$', args.fasta).groups()[0]
    if classifName == "classif" :
        print("the path of {} is {}".format(classifName, answer))
    if fastaName == "fasta":
        print("the path of {} is {}".format(fastaName, args.fasta))
except AttributeError as e:
    print("One of the extension file is incorrect")
    sys.exit
"""


####################
#       MAIN
####################

if __name__ == "__main__":
    # nom="/home/tfrances/Bureau/donnees/sortie_PASTEC/TisoMgescan.classif"
    ####    Reading of the classif file ####
    try:
        pastecFile=sys.argv[1] # take the second argument in the terminal command as file name
        pastec=readPastec(pastecFile)
    except IndexError:
        print("Pas de fichier PASTEC fourni")
        sys.exit(0)

    ####    Procede to the categorisation
    cat=categorisation(pastec)
    for key in cat.keys():
        # print(cat[key]['other'])
        pass

    ####    Reading of the fasta file ####
    try:
        fastaFile=sys.argv[2] # take the second argument in the terminal command as file name
        fasta=readFasta(fastaFile)
    except IndexError:
        print("Pas de fichier Fasta fourni")
        sys.exit(0)
    # print(fasta)
