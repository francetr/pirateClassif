# pirateClassif

## Aim of this project :

This script is a part of the PiRATE project (that aim to analyze the Transposable Element in sequences).
The purpose of this code is to automatize the step of classification of PiRATE.

For more information, you can go read the research paper that lead to the creation of PiRATE :
_http://www.seanoe.org/data/00406/51795/_

## Files used for this script:
* _pastecFile_ : Output file of the pastec tool. The script will retrieve the class, order and superfamily for each sequence in this file

* _baselineFile_ : File containing the baseline of the superFamilies names (multiple names can define the same superFamily). This is useful to determine if a sequence can be classified as multiple families or if her superFamily is not defined yet.

* _fastaFile_ : Original file containing the sequences used. Usefull to retrieve the sequences after the classification has be done


## Command line to launch the script:

>_python3 path/toward/this/script path/toward/the/classif/file path/toward/the/baseline/file path/toward/the/fasta/file_

example:
~~~
python3 code/scriptClassif.py TisoRepet1.classif base_reference.txt hat.fasta
~~~

> _./path/toward/the/scriptClassif.py path/toward/the/classif/file path/toward/the/baseline/file path/toward/the/fasta/file_

example:
~~~
./code/scriptClassif.py TisoRepet1.classif base_reference.txt hat.fasta
~~~
