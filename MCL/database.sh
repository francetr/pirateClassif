#!/bin/bash

mysql -u root --password=jeremy07 -e "drop database MCL; create database MCL"

qdel -u jeremy

if [ -d "/home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK" ];then
rm -Rf /home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK &&
mkdir /home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK
else
mkdir /home/jeremy/galaxy/tools/Pipeline/pirateClassif/MCL/WORK
fi
