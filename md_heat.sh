#!/bin/bash

# Input variables
PRMTOP=$1 # Validated parameter file from Leap and ParmEd

# Derived variables
tag=`basename $PRMTOP .prmtop`

cat > heat100K.in << EOF
Heating from 0 to 100 K in NVT
 &cntrl
  imin=0,         ! Molecular dynamics
  ntx=1,          ! Coordinates with no initial velocities
  irest=0,        ! No restart velocities assigned
  ntc=2,          ! SHAKE algorithm bonds with hydrogen
  ntf=2,          ! No force evaluation bonds with hydrogen
  tol=0.0000001,  ! SHAKE tolerance
  nstlim=50000,   ! Number of MD steps
  ntt=3,          ! Langevin dynamics
  gamma_ln=1.0,   ! Collision frequency Langevin dynamics
  ig=-1,          ! Random seed Langevin dynamics
  ntpr=500,
  ntwr=500,
  ntwx=500,       ! Write to trajectory file every ntwx steps
  dt=0.001,       ! Timestep (ps)
  ntb=1,          ! Constant volume PBCs
  ntp=0,
  cut=8.0,

  tempi=10.0,     ! Starting temperature
  TEMP0=100.0,    ! Reference temperature (held here)
 &end
END
END
EOF

cat > heat300K.in << EOF
Heating from 100 to 300 K in NPT
 &cntrl
  imin=0,         ! Molecular dynamics
  ntx=5,          ! Positions and velocities read NetCDF
  irest=1,        ! Restart calculation
  ntc=2,          ! SHAKE algorithm bonds with hydrogen
  ntf=2,          ! No force evaluation bonds with hydrogen
  tol=0.0000001,  ! SHAKE tolerance
  nstlim=500000,  ! Number of MD steps
  ntt=3,          ! Langevin dynamics
  gamma_ln=1.0,   ! Collision frequency Langevin dynamics
  ig=-1,          ! Random seed Langevin dynamics
  ntpr=2500,
  ntwr=2500,
  ntwx=2500,      ! Write to trajectory file every ntwx steps
  dt=0.001,       ! Timestep (ps)
  ntb=2,          ! Constant pressure PBCs
  ntp=1,
  cut=8.0,
 &end
 &wt
  type='TEMP0',   ! Varies the target temperature TEMP0
  istep1=50001,   ! Initial step
  istep2=250000,  ! Final step
  value1=100.0,   ! Initial temp0 (K)
  value2=300.0,   ! final temp0 (K)
 &end
 &wt
  type='TEMP0',   ! Varies the target temperature TEMP0
  istep1=250001,  ! Initial step
  istep2=500000,  ! Final step
  value1=300.0,   ! Initial temp0 (K)
  value2=300.0,   ! final temp0 (K)
 &end
 &wt
   type='END',
 &end
END
END
EOF
# Heat the system from 0 to 100k slowly over 50 ps
pmemd.cuda -O \
        -p $PRMTOP \
        -c minall_${tag}.crd \
        -i heat100K.in \
        -o heat100K_${tag}.out \
        -r heat100K_${tag}.crd \
        -x heat100K_${tag}.nc

# Heat the system from 100K to 300k slowly over 250 ps
pmemd.cuda -O \
        -p $PRMTOP \
        -c heat100K_${tag}.crd \
        -i heat300K.in \
        -o heat300K_${tag}.out \
        -r heat300K_${tag}.crd \
        -x heat300K_${tag}.nc

