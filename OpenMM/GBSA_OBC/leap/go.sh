#!/bin/bash

export AMBERHOME=/netscratch/mw529/code/AMBER/amber12

$AMBERHOME/bin/tleap -f leap.bat &> leap.log
$AMBERHOME/bin/sander -O
#$AMBERHOME/bin/ambpdb -p prmtop < restrt > out.pdb
