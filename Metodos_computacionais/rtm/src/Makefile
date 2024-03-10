#
# Pressupoe a existencia do diretorio /home/user/bin
# os arquivos executaveis sao movidos para este diretorio
#
BIN=/home/${USER}/bin
# prototipo de modelagem acustica 2d
OBJMOD=su_types_mod.o \
       su_io_mod.o \
       interp_trace.o \
       number_types.o \
       rtm_mod.o \
       afd_2d.o



# compilador intel:
#F90=ifort
#FCOPT=-O3 -assume byterecl -openmp -parallel -fpp -ipo

#compilador GNU
F90=mpif90
DEBUG=-g -assume byterecl -fpp 
SERIAL=-O3 -assume byterecl
FLAGS=$(SERIAL)
#FCOPT=-O3 -ffast-math -std=f2008 

.SUFIXES: .f90 .o .mod


afd_2d.x: $(OBJMOD)
	$(F90) $(FLAGS) $(OBJMOD) -o $(BIN)/afd_2d.x

%.o:%.f90
	$(F90) -c $(FLAGS) $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~

