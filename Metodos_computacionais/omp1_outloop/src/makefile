#
# Pressupoe a existencia do diretorio /home/user/bin
# os arquivos executaveis sao movidos para este diretorio
#
BINDIR=$(HOME)/bin
#
# prototipo de modelagem acustica 2d
OBJ=su_types_mod.o \
    su_io_mod.o \
    afd_2d_omp.o

.SUFIXES: .f90 .o .mod

#compilador GNU
DBG=-Og  -fopenmp -fcheck=all -std=f2003
RLS=-O3 -fopenmp  -std=f2003
F90=gfortran
FCOPT=$(RLS) 

afd_2d.x: $(OBJ)
	$(F90) $(FCOPT) $(OBJ) -o $(BINDIR)/afd_2d.x

%.o:%.f90
	$(F90) -c $(FCOPT) $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~

