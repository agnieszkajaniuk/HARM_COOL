#/***********************************************************************************
#  HARM NEUTRINO version 4.85   (released December, 2019)
#***********************************************************************************/
#### set USEICC to 0 if you want gcc compiler options, else set to 1 to use icc
########  gcc generally used for debugging with -g option so we can use gdb 
USEICC = 0
USEMPI = 1


ifeq ($(USEICC),1)
CC       = icc
CCFLAGS  = -O2 -axW -tpp7
#CCFLAGS  = -O2
endif

ifeq ($(USEICC),0)
CC       = g++
CCFLAGS  = -O2
endif

SPECTRUM_FC = gfortran
SPECTRUM_OPT   = -O2 -c -w 

EXTRALIBS = -lm -lgfortran


ifeq ($(USEMPI),1)
CC       = mpic++
CCFLAGS  = -O2 -Wno-unused-result -DMPI_USED -mfpmath=387 -I/usr/include/hdf5/openmpi
SPECTRUM_FC = mpif90
EXTRALIBS = -lm -lgfortran -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5 -ldl -lz -lgsl -lgslcblas
endif



CC_COMPILE  = $(CC) $(CCFLAGS) -c 
CC_LOAD     = $(CC) $(CCFLAGS) 

.SUFFIXES : .o .f90 .c

.c.o:
	$(CC_COMPILE) -o $@ $<

.f90.o:
	${SPECTRUM_FC} ${SPECTRUM_OPT} -o $@ $<

EXE = harm
#all: $(EXE) image_interp
#all: $(EXE) h5_slice
all: $(EXE)


COBJS = \
bounds.o coord.o diag.o dump.o fixup.o \
image.o interp.o main.o metric.o lu.o \
phys.o ranc.o restart.o step_ch.o u2p_util.o\
utoprim_1dfix1.o utoprim_1dvsq2fix1.o utoprim_2d.o \
cooleos.o cache4d.o akima.o random.o utoprim_5d.o mpi.o tracers.o


COBJS += init.o #original init; change here if you have other

INCS = decs.h defs.h u2p_defs.h u2p_util.h mnewt.h random.h


FOBJ = const_mod.o gauss_integral.o bisection.o mnewtf.o number_den.o\
neu_rate.o eos.o accretion.o chem_eq.o abundance.o 


$(EXE): $(COBJS) $(INCS) $(FOBJ)
	$(CC_LOAD) $(COBJS) $(FOBJ) $(EXTRALIBS) -o ../exe/$(EXE)

h5_slice: h5_slice.o
	$(CC_LOAD) h5_slice.o $(EXTRALIBS) -o ../exe/h5_slice

image_interp: image_interp.c
	$(CC_LOAD) $(EXTRALIBS) image_interp.c -o ../exe/image_interp

clean:
	/bin/rm -f *.o *.il *.mod
	/bin/rm -f ../exe/$(EXE) ../exe/image_interp ../exe/h5_slice


newrun:
	/bin/rm -rf ../exe/dumps ../exe/images ../exe/ener.out

