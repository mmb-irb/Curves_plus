## Curves+

---

[Official Homepage](https://bisi.ibcp.fr/tools/curves_plus/index.html)

_Authors: R. Lavery, M. Moakher, J. H. Maddocks, D. Petkeviciute, K. Zakrzewska_

**Curves+** is a revised version of the Curves approaech for analysing the structure of nucleic acids. It respects the international conventions for nucleic acid analysis, runs much faster and provides new data.


As well as treating single nucleic acid structures, Curves+, with the help of **Canal**, can analyse molecular dynamics trajectories and generate time series, time averaged properties and search for correlations between variables.


In conjunction with **Canion**, Curves+ can now be used to analyze the distribution of ions or molecules around nucleic acids in curvilinear helicoidal coordinates.

---

Usage:

    # compile Curves+, Canal and Canion using gfortran as compiler
    make all FC=gfortran 

    # compile only Curves+
    make Cur+ FC=gfortran
    
    # compile only Cabal (Curves+ needs to be compiled first)
    make Canal FC=gfortran

    # clean 
    make clean

    # test Cur+ execution on 1bna.pdb, should return no output if successful
    bash test/test.sh
