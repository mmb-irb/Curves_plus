#!/bin/bash
  
returncode=0

# lis
if [ -e canionout.lis ]
then
 printf "canionout.lis file successfully created\n"
else
 printf "canionout.lis file not created!\n"
 returncode=1
fi
# check number of files created by Canion
wordcount=$(ls | grep canionout | wc -l)
if [[ "$wordcount" == "       9" ]] 
then
 printf "correct number of files created by Canion\n"
else
 printf "some or all files created by Canion missing!\n"
 returncode=1
fi

exit $returncode
