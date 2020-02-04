#!/bin/bash
INPUT=$1
tag=`basename $INPUT .prmtop`
cat > parmed_check.in << EOF
CheckValidity
outparm ${tag}.valid.prmtop
EOF
parmed -p $INPUT -i parmed_check.in
mv $INPUT ${tag}.backup.prmtop
mv ${tag}.valid.prmtop ${tag}.prmtop
