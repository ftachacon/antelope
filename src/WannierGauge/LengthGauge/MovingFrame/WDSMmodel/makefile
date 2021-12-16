
#Makefile :  mpsd
#MKL path
#MKLROOT=/opt/intel
#COMPOSERLIB=$(MKLROOT)

CXXFLAGS= -O3 -xCORE-AVX2 -std=c++17  -no-multibyte-chars
INC=-I${MKL_HOME}/include
LIBS=-L${MKL_HOME}/lib/intel64 -lmkl_intel_lp64 -mkl=cluster -lmkl_core 
MPI= mpiicpc

LIBS_=-L/cobra/u/geonda/bin/lib 
INC_=-I/cobra/u/geonda/bin/include

Materials=${MyPath}/Materials

ProgramMPI=antelope.cpp

exnamempi=exec_hhg_mpi

all: compileMPI

compileMPI: ${ProgramMPI}
	${MPI} ${CXXFLAGS} ${ProgramMPI} -w   -I${Materials} ${INC_} ${INC} -o ${exnamempi} ${LIBS} ${LIBS_}  /u/geonda/bin/lib/libconfig++.a 

	cp ./exec_hhg_mpi ../../../../../bin/
clean:
	rm a.out

