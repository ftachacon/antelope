#!/bin/bash
dirname=$1 # Pass a simple name to this scrip
# pe request
#$ -pe mpich 23

# our Job name 
#$ -N antelope

#$ -S /bin/bash

#$ -q dque_ib # $ -V

#$ -cwd

# needs in 
#   $NSLOTS          
#       the number of tasks to be used
#   $TMPDIR/machines 
#       a valid machiche file to be passed to mpirun 
#   enables $TMPDIR/rsh to catch rsh calls if available

echo "Got $NSLOTS slots."
cat $TMPDIR/machines



#######################################################
### mpich 1.2.7p1 (w/ Intel-12.1 compiler)
#######################################################
#
# MPI_HOME=/opt/mpi/intel-12.1/mpich-1.2.7p1
# MPI_EXEC=$MPI_HOME/bin/mpirun
#
# cd $SGE_O_WORKDIR
#
# $MPI_EXEC -machinefile $TMPDIR/machines -np $NSLOTS $SGE_O_WORKDIR/exec_file



#######################################################
### mpich2 1.4.1p1 (w/ Intel-12.1 compiler)
#######################################################
#
# MPI_HOME=/opt/mpi/intel-12.1/mpich2-1.4.1p1
# MPI_EXEC=$MPI_HOME/bin/mpirun

# cd $SGE_O_WORKDIR

# $MPI_EXEC -machinefile $TMPDIR/machines -n $NSLOTS $SGE_O_WORKDIR/exec_file



#######################################################
### openmpi 1.4.4 (w/ Intel-12.1 compiler)
#######################################################
#
# MPI_HOME=/opt/mpi/intel-12.1/openmpi-1.4.4
# MPI_EXEC=$MPI_HOME/bin/mpirun
#
# cd $SGE_O_WORKDIR
#
# $MPI_EXEC -machinefile $TMPDIR/machines -n $NSLOTS /opt/vasp/vasp5.3.3/vasp


###############################
## Some laser field inputs ##
E0=0.0030 #07778 #0.0012963 #0.0019445        #0.007778   #0.003889         # 0.003           # Electric field strength (a.u.)
Ncycles=12          # No. of Opt. cycles
ellip=+0            # Ellipticity


#Dynamical-Time Variables or Params 
dt=0.02             # Time-Steps (a.u.)
dephasing=83. #220.      #1440. #220      # Dephasing (a.u.)

#k-Space variables or params 
Nx=513 #403 #401   #201 #201              # No. of points along kx
Ny=575 #483 #465   #233 #233             # No. of points along ky, ratio Nx/Ny=1.72 or, dpending on box 1.16
ksfactor=1 #BZ scale-factor


#########################################
### SOME INPUTs PARAMETERS
## Set of Haldane M. Parameters

phi0=0.06 #0.0 #1.16           # Magnetic flux or phase of the complex 2nd hopping (rad.)
Mt2=2.54 #0.103309 #       # Ratio of on-site potential and 2nd hopping



gauge=1        # gauge parameter, can be 0, 1 and 2
gbox_ky_down=0.01          #Botton Boundary  (lower) ky-BZ point for gauge modification
gbox_ky_up=0.70            #Upper  Boundary  (higher) ky-BZ point for gauge modification
ky_shift_down=-0.35        #Lowest BZ ky shift in a.u.
fbz_shift=1                 #Controlling BZ shift. it has to be 1 for yes, 0 for none
diag=1                      #diagnostic

###################################
mkdir ${dirname} #creating test directory


#rm exec_hhg_mpi*


####################################
make #compiling code and generating executable 



####################################
#Copying exe* and other files to Test dir.
#rm ./${dirname}/exec_hhg_mpi



cp exec_hhg_mpi ./${dirname}/
cp AlsisPythonShort/An* ./${dirname}/
cp AlsisPythonShort/postprocessingBForVis1.py ./${dirname}/



#######################################################
### openmpi 1.10.1 (w/ Intel-14.0.2 compiler)
#######################################################
#
 #MPI_HOME=/opt/mpi/intel-14.0/openmpi-1.10.1
 #MPI_HOME=/opt/mpi/intel-13.1/openmpi-1.10.1
MPI_HOME=/opt/mpi/gcc-4.9.2/openmpi-1.10.1
MPI_EXEC=$MPI_HOME/bin/mpirun
#
cd $SGE_O_WORKDIR
#
date 
$MPI_EXEC -machinefile $TMPDIR/machines -n $NSLOTS ./exec_hhg_mpi ${phi0}  ${Mt2}  ${Nx}  ${Ny}  ${E0}  ${Ncycles}  ${dt}  ${dephasing} ${ellip} ${eps} ${gauge} ${rflag} ${gbox_ky_down} ${gbox_ky_up} ${fbz_shift} ${ky_shift_down} ${diag} ${ksfactor} > stdout.dat
date 
