#!/bin/bash
  
returncode=0

# test output files were created
# lis
if [ -e canout.lis ] 
then
 printf "canout.lis file successfully created\n"
else
 printf "canout.lis file not created!\n"
 returncode=1
fi
# check number of files created by Canal
wordcount=$(ls|grep canout|wc -l)
if [[ "$wordcount" == "      45" ]] 
then
 printf "correct number of files created by Canal\n"
else
 printf "some or all files created by canal missing!\n"
 returncode=1
fi

exit $returncode
