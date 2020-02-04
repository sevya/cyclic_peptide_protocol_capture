#!/bin/bash
INPUT=$1
tag=`basename $INPUT .pdb`
cat > tleap.in << EOF
source leaprc.protein.ff14SB
source leaprc.water.tip4pew
mol = loadpdb $INPUT
saveAmberParm mol dry.${tag}.prmtop dry.${tag}.inpcrd
savepdb mol dry.${tag}.pdb
solvateBox mol TIP4PEWBOX 12
charge mol
addions2 mol Cl- 0
addions2 mol Na+ 0
saveAmberParm mol ${tag}.prmtop ${tag}.inpcrd
EOF
tleap -f tleap.in

