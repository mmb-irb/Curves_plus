#!/bin/bash

rm r+bdna*.*
cp test/1bna.pdb .
./Cur+ <<!
 &inp file=1bna, lis=r+bdna,
 lib=standard,
 &end
 2 1 -1 0 0
 1:12
 24:13
!
rm r+bdna*.*
rm 1bna.pdb
