
#MKL path
#MKLROOT=/opt/intel
#COMPOSERLIB=$(MKLROOT)

CXXFLAGS= -O2 -std=c++14
LIBS= -lconfig++ -llapacke -llapack -lblas -lm
#INTEL= icpc
MPI= /usr/bin/mpicxx
#LIBS= -lpthread -lm -ldl -liomp5 -lmkl_rt


#My library a. chacon
MyPath=.
OWN_LIB=${MyPath}/include

Materials=${MyPath}/Materials

ProgramMPI=antelope_wannier.cpp

exnamempi=exec_hhg_mpi
exname=exec_hhg

comp: compileMPI 
#comp: compileSq compileMPI 


compileSq: ${Program0} 
	$(MPI)  ${Program0} -w -I${OWN_LIB} -o ${exname}


compileMPI: ${ProgramMPI}
	${MPI} ${CXXFLAGS} ${ProgramMPI} -w -I${OWN_LIB} -I${Materials} -o ${exnamempi} ${LIBS}


clean: 
	rm a.out

#ROOT_PATH=/home/achacon
#LD_FFTW=$(ROOT_PATH)/libraries/fftw-3.3.4/LFFTW3


