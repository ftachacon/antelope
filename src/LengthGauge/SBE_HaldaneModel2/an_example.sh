#!/bin/bash


####################################################
# This script runs antelope code for a basic example.
# You need to give a single "parameter" (directory name)
# via terminal and this scrip will compute:
# 1) evolution of the Semiconductor Bloch Eqs.
#     for Haldane Model
# 2) intra and inter band currents as a function of time



dirname=$1 # Pass a simple name to this scrip




##########################################
## You can modify the No. of MPI tasks

N_MPI_TASKs=2   #No. of MPI task or cores for parallelization




#########################################
### SOME INPUTs PARAMETERS
## Set of Haldane M. Parameters

phi0=0.00       # Magnetic flux or phase of the complex 2nd hopping (rad.)
Mt2=2.54 #1. #2.54        # Ratio of on-site potential and 2nd hopping


gauge=2         # gauge parameter, can be 0, 1 and 2
rflag=0         # flag =0, No-reg.; flag=1, Taylor-Reg.; flag=2, Local-Gauss-Reg.

eps=0.1         # Reg. parameter for a fixed gauge, i.e. gauge=1


Nx=100 #120 #250  #90   #170  #  #160  #90        # No. of points along kx
Ny=58  #70 #146  #52   #70  #52  #94  58 #          # No. of points along ky, ratio Nx/Ny=1.72



###############################
## Some laser field inputs ##
E0=0.0050 #0.007            # Electric field strength (a.u.)
Ncycles=3          # No. of Opt. cycles
ellip=1.           # Ellipticity

dt=0.80             # Time-Steps (a.u.)
dephasing=110.  #1440. #220      # Dephasing (a.u.)





###################################
mkdir ${dirname} #creating test directory


rm exec*
####################################
make #compiling code and generating executable 



####################################
#Copying exe* and other files to Test dir.
cp exec_hhg_mpi ./${dirname}/
cp AlsisPythonShort/An* ./${dirname}/
cp AlsisPythonShort/postprocessingBForVis1.py ./${dirname}/



###################################
#Running test example
cd ./${dirname}


date
ls -lh exe*
time -p mpirun -n ${N_MPI_TASKs} exec_hhg_mpi ${phi0}  ${Mt2}  ${Nx}  ${Ny}  ${E0}  ${Ncycles}  ${dt}  ${dephasing} ${ellip} ${eps} ${gauge} ${rflag}

date 


#############################################
## Next line performe remain momentum integral
## along ky-direction #
##

echo "--------------------"
echo ""
echo "Integration along ky direction via python script"

time -p python postprocessingBForVis1.py / 0 ${N_MPI_TASKs}

echo ""
echo "End of the integration along ky-direction"


##removing temporary files
rm full_*


###################################
##Analysis-first-step 
#cd  ./AlsisPythonShort/



######################################
##Convering data to python analysis codee
#./script_to_analysis.sh



cd ../
echo ""
echo "" 
