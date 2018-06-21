# pirateClassif

## Aim of this project :

This script is a part of the PiRATE (**Pi**peline to **R**etrieve and **A**nnotate **T**ransposable **E**lements) project (that aim to analyze the Transposable Elements in sequenced genomes). The purpose of this code is to automatize the step of classification of PiRATE.  
You can have acces to the VM of PiRATE to the page:  _http://www.seanoe.org/data/00406/51795/_

For more information, you can read the research paper that lead to the creation of PiRATE :  
_https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4763-1_

## Requirements
This script is written in python, it can be launched either with python2.7 or python3.5. See the official site of Python for more informations about it : _https://www.python.org/_  
However, it requires the ***biopython*** modules (usefull to write the sequences onto fasta files), which you can acces with ***pip***.  
To install Biopython, we invited you to do it via a terminal with the command:
~~~{bash}
Installation of pip:  
For python2:  
sudo apt-get install python-pip  

For python3:  
sudo apt-get install python3-pip  

Then you can install biopython:  
For python2:  
pip install biopython  

For python3:  
pip3 install biopython  
~~~

To be sure biopython is update
~~~
For python2:  
pip install biopython --upgrade  

or

For python3:  
pip3 install biopython --upgrade
~~~

If you got troubles with the installation of biopython or if you need to install it with other OS than Linux please check its official website : _https://biopython.org/wiki/Download_  

## Files used for this script:
* ***pastecFile*** : Output file of the pastec tool. This file must be correctly written, meaning all columns are separated by one or two tabulation characters.  
The script will retrieve the class, order and superfamily for each sequence in this file.
> Example of a line:  
> repeatscout_443	2567	+	PotentialChimeric	II	TIR	incomplete	CI=62; coding=(TE_BLRtx: VANDAL2:ClassII:TIR:MuDR: 16.84%; TE_BLRx: MuDR-1_ALy_3p:ClassII:TIR:MuDR: 96.44%, MuDR-1_ALy_4p:ClassII:TIR:MuDR: 100.00%; profiles: PF02902.14_Peptidase_C48_NA_CYP_20.5: 99.07%(99.07%)); struct=(TElength: >1000bps); other=(SSRCoverage=0.13)


* ***baselineFile*** : File containing the baseline of the superFamilies names (multiple names can define the same superFamily). This is useful to determine if a sequence can be classified as multiple families or if its superFamily is not defined yet. By default it is ***base&#95;reference.txt*** which is used.  
For the classification step, we used ***blast*** and ***protProfiles*** keywords. The first concern keywords that can be find in the ***TE&#95;BLRx*** or ***TE&#95;BLRtx*** part of coding, and the last can be found in the ***profiles*** part.  
***NB:&nbsp;*** It is important that the baseline file is located in the same directory than the current working directory where you launched the script. If it is not the case, you still can specify its path with the argument ***--baseline***.  
This file can be completed, for this, you just need to enter the new family name in a new line followed by a _tabulation_ and the different names (separated by ***:*** ) that can design this superfamily.  
> Example :  
> Mariner&nbsp;&nbsp;&nbsp;&nbsp;Tc1-mariner:mariner  
> &gt;TASE&nbsp;&nbsp;&nbsp;&nbsp;Mariner:hAT:PIF-Harbinger:PiggyBac:Merlin:Transib:CACTA:P


* ***fastaFile*** : File containing the fasta sequences used to run the PASTEC tool. Usefull to retrieve the sequences after the classification has been done.


## Command line to launch the script:
Four arguments can be passed onto the command line but two arguments are mandatory, the third and the fourth are optionnal.
1. The first is the **output** file of PASTEC, usually with an extension ***.txt***
2. The second is the **fasta** file used for PASTEC, usually with an extension ***.fasta***
3. The third is the identity threshold (in percentage) by which the superFamily name will be determined for a sequence. By default it will be ***100&#37;*** but it can be changed with the argument ***-e***.  
***NB:&nbsp;*** By default, if you enter a percentage > 100, it will be 100; and for a perntage < 0, it will be 0.
4. The fourth is the **baseline** file, by default it will be ***base&#95;reference.txt*** but it can be changed with the argument ***--baseline*** followed by the path where the new baseline file is.

Hence multiple command line can be used to launch the script with Python2 or Python3
> python&nbsp;&nbsp;path/toward/this/script&nbsp;&nbsp;path/toward/the/classif/file  path/toward/the/fasta/file   python3&nbsp;&nbsp;path/toward/this/script&nbsp;&nbsp;path/toward/the/classif/file  path/toward/the/fasta/file  
> example:


~~~{bash}
python src/scriptClassif.py ArabiTEdenovo.txt ArabiTEdenovo.fasta
python3 src/scriptClassif.py ArabiTEdenovo.txt ArabiTEdenovo.fasta
python3 src/scriptClassif.py ArabiTEdenovo.txt ArabiTEdenovo.fasta -e 75
python3 src/scriptClassif.py ArabiTEdenovo.txt ArabiTEdenovo.fasta --baseline new_base_reference.txt
python3 src/scriptClassif.py ArabiTEdenovo.txt ArabiTEdenovo.fasta -e 75 --baseline new_base_reference.txt
~~~


</br>
It can also be launched by executing itself (be sure the permission for execution are granted for the script)

> ./path/toward/the/scriptClassif.py&nbsp;&nbsp;path/toward/the/classif/file path/toward/the/fasta/file  
> example:

~~~{bash}
./src/scriptClassif.py ArabiTEdenovo.txt ArabiTEdenovo.fasta
~~~

</br>
For help, you can type the command

~~~{bash}
python3 scriptClassif.py -h
~~~


## What this script do?
This script will retrieve interesting informations from the output file of PASTEC (namely the class, order and superFamily of the sequence).
Before these informations are determined, the type of the sequence is first established, four outcomes are possible:
1. **noCat** : the sequence can't be categorized
2. **nonTE** : the sequence is not a Transposable Element (can be PotentialHostGene)
3. **PotentialChimeric** : the sequence can be categorized as 2 or more superFamilies
4. **TE** : the superFamily's sequence has been established

According the type of the sequence, the script will determine the class, order, and superFamily of the pastec output, and saves the sequence in accordance with th its outcome corresponding (names of the outcomes possibles are cited above).  
In addition with these four outcomes, two additionnals files are available:
* **Classification_summary.txt** : contains for every sequences the proofs founded to classifie the class, order and superFamily of this sequence.  
Here is an example of a line that can be found in this file :  
> ltrharvest_941 &emsp; 10243 &emsp; TE &emsp; I &emsp; LTR &emsp; BLAST : {'Gypsy': &nbsp; [6, 100.0]} &emsp; PROTPROFILES : &nbsp; {'RH': &nbsp; [1, 20.0], 'AP': &nbsp; [1, 20.0], 'RT': &nbsp; [2, 40.0], 'GAG': &nbsp; [1, 20.0]} &emsp; PREDICTED_SUPERFAMILY: &emsp; Gypsy  
> From right to left we got :  
name of the sequence &emsp; length &emsp; type of the sequence(TE, nonTE, potentialChimeric or noCat) &emsp; class &emsp; order &emsp; Blast keywords: &nbsp; name: [nb of occurrence, &nbsp; percentage], &nbsp; Proteine Profile keywords: name:[nb of occurrence, percentage] &emsp; name predicted of the superFamily

* **unknown_keywords.txt** : contains the string which doens't contains keywords that matches with the baseline (during the superFamily determination).

## How the superFamily determination work:
The diagram above show the diffent possibilities which can be founded during the keywords comparison used for the superFamily determination.
![Comparison](./diagrammes/Diagramme_comparaison_superfamille_en.jpeg)
