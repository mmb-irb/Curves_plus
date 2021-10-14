#!/bin/bash
  

export DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib

# -f Curves+ executable
Cur+ <<!
&inp 
 file=abcstride.trj,
 ftop=abctop_nowat.prmtop,
 lis=curout_ions,
 lib=$CONDA_PREFIX/.curvesplus/standard,
 line=.t.,
 fit=.t.,
 test=.t.,
 ions=.t.,
&end        
2 1 -1 0 0                                               
1:18 
36:19     
! 
