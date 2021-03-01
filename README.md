## Curves+ v3.0

---

[Official Homepage](https://bisi.ibcp.fr/tools/curves_plus/index.html)
[BSC Curves+ Server Homepage](http://curvesplus.bsc.es/misc)

_Authors: R. Lavery, M. Moakher, J. H. Maddocks, D. Petkeviciute, K. Zakrzewska_

**Curves+**
Conformational analysis of single nucleic acid structures or of molecular dynamics (MD) trajectories

**Canal**
Statistical analysis of nucleic conformational data from MD trajectories produced by Curves+

**Canion**
Analysis of ion/water/ligand distribution data from MD trajectories produced by Curves+
---

Usage:

    # compile Curves+, Canal and Canion using gfortran as compiler
    make all FC=gfortran 

    # compile only Curves+
    make Cur+ FC=gfortran
    
    # compile only Canal (Curves+ needs to be compiled first). Same procedure for Canion
    make Canal FC=gfortran

    # clean 
    make clean
