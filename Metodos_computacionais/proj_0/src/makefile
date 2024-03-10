#
# Pressupoe a existencia do diretorio /home/user/bin
# os arquivos executaveis sao movidos para este diretorio
#
BINDIR=$(HOME)/bin
#
# prototipo de modelagem acustica 2d
OBJMOD=su_types_mod.o \
       su_io_mod.o \
       afd_2d.o

.SUFIXES: .f90 .o .mod

#compilador GNU
F90=gfortran
FCOPT=-O3 -fopenmp  -std=f2003 

afd_2d.x: $(OBJMOD)
	$(F90) $(FCOPT) $(OBJMOD) -o $(BINDIR)/afd_2d.x

%.o:%.f90
	$(F90) -c $(FCOPT) $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~

