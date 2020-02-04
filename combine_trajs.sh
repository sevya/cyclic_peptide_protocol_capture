#!/bin/bash
PRMTOP=$1
tag=`basename $PRMTOP .prmtop`
cat >> combine_trajs.in < EOF
parm $PRMTOP
trajin ${tag}_prod.*.nc
autoimage
rms @CA
strip :WAT,Na+,Cl- outprefix dry.${tag}.pdb
trajout ${tag}.all.nc
run
quit
EOF
cpptraj -i combine_trajs.in

