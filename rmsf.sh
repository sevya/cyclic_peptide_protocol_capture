#!/bin/bash
# Input variables
PRMTOP=$1
TRAJ=$2
# Derived variables
tag=`basename $PRMTOP .prmtop`
cat > cpptraj_rmsf.in << EOF
parm $PRMTOP
trajin $TRAJ
autoimage
rms @CA,C,O,N
average crdset AvgStructReceptor
createcrd MyTrajFramesReceptor
run
rms ref AvgStructReceptor @CA,C,O,N
atomicfluct out ${tag}_heavyatom_rmsf.dat @CA,C,O,N byres
run
quit
EOF
cpptraj -i cpptraj_rmsf.in

