
#MKL path
MKLROOT=/opt/intel
COMPOSERLIB=$(MKLROOT)

CFLAGS= -O3 -openmp
INTEL= icpc
MPI= /usr/local/bin/mpicxx
LIBS= -lpthread -lm -ldl -liomp5 -lmkl_rt


#My library a. chacon
MyPath=.
OWN_LIB=${MyPath}/include


ProgramMPI=hhg_rk5_hm_mpi_full_vel_nc.cpp 
exnamempi=exec_hhg_mpi


comp: compileMPI 



compileMPI: ${ProgramMPI}
	${MPI} ${ProgramMPI} -w -I${OWN_LIB} -o ${exnamempi}


clean: ${exnamempi} 
	rm ${exnamempi}

#ROOT_PATH=/home/achacon
#LD_FFTW=$(ROOT_PATH)/libraries/fftw-3.3.4/LFFTW3


