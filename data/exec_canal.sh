#!/bin/bash
  
# test Canal executable
Canal <<! 
&inp
  lis=canout,
  series=.t.,
&end
curout_ions.cda GCAACGTGCTATGGAAGC       
!

