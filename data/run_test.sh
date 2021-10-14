#!/bin/bash


bash exec_curves.sh
bash check_curves_output.sh
bash exec_canal.sh
bash check_canal_output.sh
bash exec_canion.sh
bash check_canion_output.sh

# cleanup 
rm canout* curout_ions.* canionout.*
