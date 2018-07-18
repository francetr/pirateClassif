#!/usr/bin/python2
# coding: utf8

""" @author: Tristan Frances """

import os
import glob
import subprocess

def launchMCL():
	"""

	Launch the MCL on all the prelibraries created by the save module in order to obtain superfamilies. For this, it execute bash script that will do the work.

	@rtype: None
	"""
	#### This part is dedicated to the family classification using the MCL of REPET package
	#### First made the script executable
	subprocess.call(["chmod", "0744", "/home/jeremy/galaxy/tools/Pipeline/classification/MCL/database.sh"])
	#### Then execute the script
	subprocess.Popen("/home/jeremy/galaxy/tools/Pipeline/classification/MCL/database.sh", shell=True).wait()

	#### First made the script executable
	subprocess.call(["chmod", "0744", "/home/jeremy/galaxy/tools/Pipeline/classification/MCL/launch.sh"])
	#### Then execute the script
	subprocess.Popen("/home/jeremy/galaxy/tools/Pipeline/classification/MCL/launch.sh", shell=True).wait()

def saveLibraries():
    listFiles=[]
    for (rep, subRep, files) in os.walk("/home/jeremy/galaxy/tools/Pipeline/classification/MCL/WORK"):
        listFiles.append(files)
