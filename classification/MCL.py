#!/usr/bin/python2
# coding: utf8

""" @author: Tristan piratClassif """

import os
import glob
import subprocess
import re

def launchMCL():
	"""

	Launch the MCL on all the prelibraries created by the save module in order to obtain superfamilies. For this, it execute bash script that will do the work.

	@rtype: None
	"""
	#### This part is dedicated to the family classification using the MCL of REPET package
	#### First made the script executable
	subprocess.call(["chmod", "0744", "/home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/database.sh"])
	#### Then execute the script
	subprocess.Popen("/home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/database.sh", shell=True).wait()

	#### First made the script executable
	subprocess.call(["chmod", "0744", "/home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/launch.sh"])
	#### Then execute the script
	subprocess.Popen("/home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/launch.sh", shell=True).wait()

def saveLibraries():
	"""

	Creates the three libraries, according the research paper PiRATE :
		- 1 for autonomous TE
		- 1 for all TE
		- 1 for all repeated elements
	For this, it uses the sequences generated by the MCL script.

	@rtype: None
	"""
	#TODO
	listFiles=[]
	listClassification=[]
	listNonAutonomous=["TRIM","LARD","SINE","MITE","7SL","5SL","tRNA"]

	#### Use the find bash command to retrieve all the classification name of the libraries
	bashCommand = subprocess.Popen('find /home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK/ -name "*_MCL.fa" 2> /dev/null', shell=True, stdout=subprocess.PIPE);
	MCLsequences = bashCommand.stdout.read().decode("utf-8")
	#### Add the classification name in a list
	for file in MCLsequences.strip().split("\n"):
		listClassification.append(os.path.basename(file).split("_MCL.fa")[0])
	for file in listClassification:
		if file in listNonAutonomous or file in ["%s_undefined"%(na) for na in listNonAutonomous]:
			print(file)
