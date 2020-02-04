#!/bin/bash

# Input variables
PRMTOP=$1 # Validated parameter file from Leap and ParmEd
NRES=$2 # Number of residues in protein

# Derived variables
tag=`basename $PRMTOP .prmtop`
wationstart=$((NRES + 1))

# Minimization script with protein restraints
rm –f minsolv.in
cat > minsolv.in << EOF
Buffer minimization
 &cntrl
  ntx=1,irest=0,ntrx=1,ntxo=1,
  ntpr=50,ntwx=0,ntwv=0,ntwe=0,

  ntf=1,ntb=1,
  es_cutoff=8.0,
  vdw_cutoff=8.0,

  ibelly=0,ntr=1,

  imin=1,
  maxcyc=1500,
  ncyc=500,
  ntmin=1,dx0=0.1,drms=0.0001,
  ntc=1,tol=0.00001,
 &end
Hold protein fixed
10.0
RES   1 $NRES
END
END
EOF

# Minimization with buffer restraints
rm –f minpro.in
cat > minpro.in << EOF
Protein minimization
&cntrl
  ntx=1,irest=0,ntrx=1,ntxo=1,
  ntpr=50,ntwx=0,ntwv=0,ntwe=0,

  ntf=1,ntb=1,
  es_cutoff=8.0,
  vdw_cutoff=8.0,

  ibelly=0,ntr=1,

  imin=1,
  maxcyc=1500,
  ncyc=500,
  ntmin=1,dx0=0.1,drms=0.0001,
  ntc=1,tol=0.00001,
 &end
Hold buffer fixed
5.0
RES   $wationstart 9999999
END
END
EOF

# Minimization with no restraints
rm –f minall.in
cat > minall.in << EOF
All minimization
&cntrl
  ntx=1,irest=0,ntrx=1,ntxo=1,
  ntpr=50,ntwx=0,ntwv=0,ntwe=0,

  ntf=1,ntb=1,
  es_cutoff=8.0,
  vdw_cutoff=8.0,

  ibelly=0,ntr=0,

  imin=1,
  maxcyc=1500,
  ncyc=500,
  ntmin=1,dx0=0.1,drms=0.0001,
  ntc=1,tol=0.00001,
 &end
EOF
# Minimize the solvent
pmemd.cuda -O \
        -i minsolv.in \
        -o minsolv.out \
        -p $PRMTOP \
        -c ${tag}.inpcrd \
        -ref ${tag}.inpcrd \
        -r minsolv_${tag}.crd

# Minimize the protein
pmemd.cuda -O \
        -i minpro.in \
        -o minpro.out \
        -p $PRMTOP \
        -c minsolv_${tag}.crd \
        -ref minsolv_${tag}.crd \
        -r minpro_${tag}.crd

# Minimize everything
pmemd.cuda -O \
        -i minall.in \
        -o minall.out \
        -p $PRMTOP \
        -c minpro_${tag}.crd \
        -r minall_${tag}.crd

