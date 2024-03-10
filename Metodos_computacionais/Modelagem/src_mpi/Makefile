#
# Pressupoe a existencia do diretorio /home/user/bin
# os arquivos executaveis sao movidos para este diretorio
#
BINDIR=$(HOME)/bin
#
#
FFTWROOT=/opt/fftw
FFTWINC=$(FFTWROOT)/include
FFTWLIB=$(FFTWROOT)/lib
FFTWFLAGS=-I$(FFTWINC) -L$(FFTWLIB) -lfftw3
#
MPIROOT=/opt/mpich3
MPIF90=$(MPIROOT)/bin/mpif90
MPIINC=$(MPIROOT)/include
#compilador GNU
F90=gfortran
#FCOPT=-O3 -fopenmp -fcheck=all -ffast-math -std=f2008 -I$(MPIINC)
FCOPT=-Og -fcheck=all -std=f2003 -I$(MPIINC)
#
AFDOBJ=su_types_mod.o \
       su_io_mod.o \
       assertion_mod.o \
       afd2d_mod.o \
       afd2d_mpi.o

.SUFIXES: .f90 .o .mod


HELP:
	@echo "Introducao a modelagem e imageamento com equacao da onda"
	@echo "Prototipo para modelagem acustica 2d usando MPI"
	@echo "Para gerar todos os executaveis: make all"
	@echo " "
	@echo "executavel afd2d_mpi.x  sao criados no diretorio /home/user/bin"
all:
	make afd_mpi.x
	make clean

afd_mpi.x: $(AFDOBJ)
	$(MPIF90) $(FCOPT) $(AFDOBJ) -o $@
	ln -sf $(PWD)/$@  $(BINDIR)/$@

%.o:%.f90
	$(MPIF90) -c $(FCOPT)  $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~

