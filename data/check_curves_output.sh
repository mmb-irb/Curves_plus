#!/bin/bash
  
returncode=0

# check output files were created
# cda
if [ -e curout_ions.cda ]
then
 printf "curout_ions.cda file successfully created\n"
else
 printf "curout_ions.cda file not created!\n"
 returncode=1
fi
# cdi
if [ -e curout_ions.cdi ]
then
 printf "curout_ions.cdi file successfully created\n"
else
 printf "curout_ions.cdi file not created!\n"
 returncode=1
fi
# lis file
if [ -e curout_ions.lis ]
then
 printf "curout_ions.lis file successfully created\n"
else
 printf "curout_ions.lis file not created!\n"
 returncode=1
fi

exit $returncode

