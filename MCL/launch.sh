#!/bin/bash

cd /home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK/&&

export REPET_HOST=localhost &&
export REPET_USER=root &&
export REPET_PW=jeremy07 &&
export REPET_DB=MCL &&
export REPET_PORT=3306 &&

export REPET_PATH=/home/jeremy/Pipeline/REPET_2.5 &&
export PYTHONPATH=/home/jeremy/Pipeline/REPET_2.5 &&
export REPET_JOBS=MySQL &&
export REPET_JOB_MANAGER=SGE &&
export REPET_QUEUE=SGE &&
export PATH=/home/jeremy/Pipeline/REPET_2.5/bin:/bin:/home/jeremy/perl5/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games &&

#### OLD VERSION
#### Find all the fasta sequences of TE prelibraries and copy them into MCL directory
FILESPATH=`find /home/jeremy/galaxy/tools/Pipeline/pirateClassif/classification/classification_result/prelibraries/TE -name "*.fasta"` &&

#### Loop to copy the FASTA prelibraries, created after the scriptClassif, into the MCL directory.
for FILE in $FILESPATH
do
  #### retrieve the FILE name without the extension
	SUFFIXNAME=`basename $FILE .fasta` &&
  #### Create a repository for each sequence where the MCL result will be store
  mkdir -p "/home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK/$SUFFIXNAME/" &&
  #### Go to the directory where the MCL will be apply
	cd "/home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK/$SUFFIXNAME"  &&
	#### First remove short sequences with cd-hit-est
	cd-hit-est -aS 1 -c 0.9 -g 1 -r 1 -i $FILE -o tmp.fasta
	#### Rename the prelibraries files as it should be
	mv tmp.fasta `basename $FILE`
	# cp $FILE /home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK/$SUFFIXNAME/
done

#### retrieve all the sequences that will be used for MCL
MCLFILESPATH=`find /home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK -name "*.fasta"` &&

#### Loop to apply the MCL to every file containing FASTA sequences
for FILE in $MCLFILESPATH
do
	#### retrieve the FILE name without the extension
	SUFFIXNAME=`basename $FILE .fasta` &&
  #### Apply fasta_formatter to format the max length of all the files containing FASTA sequences to 60 nt (else MCL will not work)
  fasta_formatter -i $FILE -o ${FILE}_new -w 60 &&
  #### Rename the fasta FILE into the normal name
  mv ${FILE}_new $FILE &&

  #### Change the working directory into the directory containing the FASTA sequence, else MCL won't work
	cd /home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK/$SUFFIXNAME
	#### Apply the MCL for each files containing FASTA sequence
	/home/jeremy/Pipeline/REPET_2.5/bin/PostAnalyzeTELib.py -a 1 -i `basename $FILE` -M MCL -v 5 &&
	perl -ne 'print if s/^>(.*)_(MCL)(\d+)(.*)/$2\t$3\t$1$4/' ${SUFFIXNAME}_MCL.fa |sort -n -k2 |perl -ane 'print $F[0].$F[1]."\t".$F[2]."\n"' > Clustering_${SUFFIXNAME}_MCL.txt
done

# /home/jeremy/Pipeline/REPET_2.5/bin/PostAnalyzeTELib.py -a 1 -i "genome.fasta" -M MCL -v 5 &&
# perl -ne 'print if s/^>(.*)_(MCL)(\d+)(.*)/$2\t$3\t$1$4/' genome_MCL.fa |sort -n -k2 |perl -ane 'print $F[0].$F[1]."\t".$F[2]."\n"' > Clustering_MCL.txt
