#!/bin/bash

# Input variables
PRMTOP=$1 # Validated parameter file from Leap and ParmEd

# Derived variables
tag=`basename $PRMTOP .prmtop`

# Generate production file
cat > prod.in << EOF
 Production, 300K 10ns
 &cntrl
  imin=0,          ! Molecular dynamics
  ntx=5,           ! Positions and velocities read NetCDF
  irest=1,         ! Restart calculation
  ntc=2,           ! SHAKE on for bonds with hydrogen
  ntf=2,           ! No force evaluation for bonds with hydrogen
  tol=0.0000001,   ! SHAKE tolerance
  nstlim=5000000,  ! Number of MD steps - 10 ns at 2 fs timestep
  ntt=3,           ! Langevin dynamics
  gamma_ln=2.0,    ! Collision frequency for Langevin dyn.
  temp0=300.0,     ! Simulation temperature (K)
  ntpr=2500,       ! Print to mdout every ntpr steps
  ntwr=2500,       ! Write a restart file every ntwr steps
  ntwx=2500,       ! Write to trajectory file every ntwx steps
  ntwprt=282,      ! Atoms (starting at 1) in output trajectory
  dt=0.002,        ! Timestep (ps) is 2 fs
  ig=-1,           ! Random seed for Langevin dynamics
  ntb=2,           ! Constant pressure PBCs
  ntp=1,           ! Isotropic pressure coupling
  cut=8.0,         ! Nonbonded cutoff (Angstroms)
  ioutfm=1,        ! Write binary NetCDF trajectory
  ntxo=2,          ! Write binary restart file
  barostat=2,      ! Use Monte Carlo barostat (Amber 14+)
  iwrap=1,
 &end
EOF

# Start production with output from heating run
cp heat300K_${tag}.crd ${tag}_prod.0000.crd
# Split trajectories up into 10 ns pieces
for run in `seq 1 400`; do
        pcnt=$((run - 1))
        prev=`echo $pcnt | awk '{printf( "%04d", $0)}'`
        run=`echo $run | awk '{printf( "%04d", $0)}'`
        if [[ -e ${tag}_prod.${prev}.crd ]]; then
                echo -n "Starting pmemd.cuda at "
                date
                pmemd.cuda -O \
                        -i prod.in \
                        -p $PRMTOP \
                        -c ${tag}_prod.${prev}.crd \
                        -o ${tag}_prod.${run}.out \
                        -r ${tag}_prod.${run}.crd \
                        -x ${tag}_prod.${run}.nc \
                        -inf ${tag}_prod.${run}.mdinfo
                echo -n "\npmemd.cuda ended at "
                date
        fi
done

