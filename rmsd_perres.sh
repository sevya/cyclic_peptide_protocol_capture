#!/bin/bash
# Input variables
PRMTOP=$1
TRAJ=$2
REF=$3
# Derived variables
tag=`basename $TRAJ .nc`
cat > cpptraj_rmsd_perres.in << EOF
parm $PRMTOP
trajin $TRAJ
autoimage
reference $REF
rms reference @CA,C,O,N perres perresout ${tag}_heavyatom.dat perresavg ${tag}_perres_heavyatom.dat
run
quit
EOF
cpptraj -i cpptraj_rmsd_perres.in

