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
def readPastec(PASTEC, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
    """
    Read a PASTEC file as input line by line, and then proceed to the categorization of each sequence

    Keyword argument
    PASTEC -- name of the classif file that will be opened
    NONTE -- dictionnary for non transposable element (nonTE)
    POTENTIALCHIMERIC -- dictionnary for potential chimeric elemenet
    NOCAT -- dictionnary for non categorized elemenet (noCat)
    TE -- dictionnary for transposable elemenet (I or II)

    """
    try:
        with open(PASTEC, "r") as f:
            for line in f:
                sequence=line.replace("\t\t", "\t").strip()
                categorization(sequence, NONTE, POTENTIALCHIMERIC, NOCAT, TE)
    except FileNotFoundError as e:
        print("No such file as {}".format(PASTEC))
        sys.exit()
        # string=f.readlines().replace("\t\t", "\t").strip() # in the order : read the file; replace the first separator; delete the blank at the end of file

        # return(string)

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
    try:
        with open(fastaFile, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                return(record)
    except(FileNotFoundError, NameError):
        print("No such file as {}".format(FASTA))
        sys.exit()

def categorization(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
    """
    Categorize the Transposable Element from the input argument

    Return a dictionnary which contains the different catagories which caracterize the sequence.

    Keyword argument
    SEQUENCE -- name of the string, which contain the sequence to categorize, that will be parsed
    NONTE -- dictionnary for non transposable element (nonTE)
    POTENTIALCHIMERIC -- dictionnary for potential chimeric elemenet
    NOCAT -- dictionnary for non categorized elemenet (noCat)
    TE -- dictionnary for transposable elemenet (I or II)

    """
    # NB : there seems to be \r\n character due to windows but it don't change the processing of the string if we split with \n or \r\n
    # try :
    #     SEQUENCE = SEQUENCE.split("\n") # split the original file into a list with each lines
    # except :
    #     pass
    listOrder=["SSR", "noCat", "LTR", "TRIM", "PotentialHost", "TIR", "DIRS", "LARD", "LINE", "MITE", "Helitron", "Maverick", "SINE"]
    features=SEQUENCE.split("\t")
    #TODO Treatment of the superfamily
    # print(tmp[7])
    classDetermination(features, nonTE, potentialChimeric,  noCat, TE)
    # TE[tmp[0]]={"class":tmp[4],"order":tmp[5], "other":tmp[7], "superfamily":"TODO", "autonomous":"TODO", "sequence":"TODO"}
    # print(len(PASTEC))



def classDetermination(SEQUENCE, NONTE, POTENTIALCHIMERIC, NOCAT, TE):
    """
    Determine the class of a sequence.
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
    # typeOrder=None
    if(SEQUENCE[4] == "I"):
        # print("I")
        TE[SEQUENCE[0]]={"class":"I"}
        # print(SEQUENCE[0])
        orderDetermination(SEQUENCE, TE)
        pass
    elif(SEQUENCE[4] == "II"):
        TE[SEQUENCE[0]]={"class":"II"}
        orderDetermination(SEQUENCE, TE)
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

def orderDetermination(SEQUENCE, TE):
    """
    Determine the order of a sequence.
    Doing so it complete a dictionnary containing the transposable element

    Keyword argument
    SEQUENCE -- list of features for one sequence.
    TE -- dictionnary for transposable elemenet (I or II)
    """
    TE[SEQUENCE[0]]["order"]=SEQUENCE[5]
    superFamilyDetermination(SEQUENCE, TE)

def superFamilyDetermination(SEQUENCE, TE):
    """
    Categorize the Transposable Element from the input argument.
    For this, complete four dictionnaries, passed onto arguments, that will contain the different catagories that caracterize the sequence.
    """
    pass


def save(DIC):
    #TODO
    """

    Save the sequence in a file

    Return a dictionnary which contains the different catagories which caracterize the sequence.

    Keyword argument
    DIC -- dictionnary containing the different sequences categorisation

    """
    return




####################
#       MAIN
####################

if __name__ == "__main__":
    # nom="/home/tfrances/Bureau/donnees/sortie_PASTEC/TisoMgescan.classif"

    ###     Instanciation of dictionnaries that will contain the results
    nonTE, potentialChimeric,  noCat, TE = {}, {}, {}, {}

    ###     Mananage the 2 arguments (PASTEC and FASTA file name) when the command is launched
    parser = argparse.ArgumentParser()
    parser.add_argument("classif", type=str, help="classif file providing from PASTEC")
    parser.add_argument("fasta", type=str, help="fasta file providing the sequence")
    args = parser.parse_args()
    try:
        classifName=re.match(r'[\S]*[/]?[\w.]+(classif)$', args.classif).groups()[0] ### regex checking is classif file as good extension
        fastaName=re.match(r'[\S]*[/]?[\w.]+(fasta)$', args.fasta).groups()[0] ### regex checking is fasta file as good extension
        if classifName == "classif" :
            # print("the path of {} is {}".format(classifName, answer))
            pass
        if fastaName == "fasta":
            # print("the path of {} is {}".format(fastaName, args.fasta))
            pass
    except AttributeError as e:
        print("One of the extension file is incorrect")
        sys.exit()


    ####    Reading of the classif file ####
    try:
        pastecFile=args.classif # take the second argument in the terminal command as file name
        pastec=readPastec(pastecFile, nonTE, potentialChimeric,  noCat, TE)
    except IndexError:
        print("No PASTEC provided")
        sys.exit()



    print("TE : %d, noCat : %d, nonTE : %d" %(len(TE), len(noCat), len(nonTE)))
    ####    Procede to the categorisation
    # cat=categorization(pastec)

    for key in TE.keys():
        print(TE[key])
    #
    # ####    Reading of the fasta file ####
    # try:
    #     fastaFile=sys.argv[2] # take the second argument in the terminal command as file name
    #     fasta=readFasta(fastaFile)
    # except IndexError:
    #     print("No Fasta provided")
    #     sys.exit()
    # # print(fasta)
