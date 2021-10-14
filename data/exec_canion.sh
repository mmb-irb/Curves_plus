#!/bin/bash

# test Canion executable
Canion <<!
&inp
 lis=canionout,
 axfrm=curout_avg,
 dat=curout_ions.cdi,
 solute=abcavg,
 type=K,
&end
!

