# some system dependent settings...

ifndef SYSTEM
SYSTEM = linux
endif

ifeq ($(SYSTEM),linux)
DBG=
OPENMP= --openmp
FFLAGS=$(DBG) -O3 $(OPENMP) -fPIC
BLAS = -lblas
LAPACK = -llapack
PREO = 
POSTO = $(BLAS) $(LAPACK)
LDFLAGS = -shared
SHELL = /bin/sh
OBJSUF=o
MODSUF=mod
FC = gfortran
MEX=mex
MWRAP=mwrap
MWFLAGS=-c99complex
endif

ifeq ($(SYSTEM),mac)
DBG=
OPENMP=
FFLAGS=$(DBG) -O3 $(OPENMP) -fPIC
PREO = -framework accelerate
POSTO =
LDFLAGS = -shared
SHELL = /bin/sh
OBJSUF=o
MODSUF=mod
FC = gfortran
MEX=mex
MWRAP=mwrap
MWFLAGS=-c99complex
endif

# locations

SRC_DIR = src
BIN_DIR = bin
TEST_DIR = ye-olde-fortran-test
TMP_DIR = tmp
MDIR = matlab
MWRAP_DIR = mwrap

VPATH = $(SRC_DIR):$(BIN_DIR)

# names

LIBBASE = bhdeig
LIBNAME  = lib$(LIBBASE).so
LIBLINK = -l$(LIBBASE)
#FSOURCES = $(shell echo $(SRC_DIR)/*.f)
#F90SOURCES = $(shell echo $(SRC_DIR)/*.f90)
#OBJS = $(patsubst $(SRC_DIR)/%.f,%.o,$(FSOURCES)) $(patsubst $(SRC_DIR)/%.f90,%.o,$(F90SOURCES))

GATEWAY = $(LIBBASE)gateway
MWRAPFILE = bhdeig

GATEWAY2 = chunksgateway
MWRAPFILE2 = chunks

GATEWAY3 = legeexpsgateway
MWRAPFILE3 = legeexps

MODS =
OBJS = hank103.o \
	cdjseval2d.o \
	helmbhrouts.o hbhrouts2d.o\
	helmrouts2d.o \
	prini.o \
	chunks_helmstokes_kerndefs.o \
	chunks_helmstokes_matbuild.o \
	chunks_helmstokes_kerndefs_fark.o \
	chunks_ders.o \
	dlaran.o \
	chunks.o \
	chunks_multi.o \
	chunks_multi_extras.o\
	hkrand.o \
	legeexps.o \
	pplot.o \
	cmat_build.o cmat_build_vec.o\
	chunks_quads.o \
	corners.o \
	adapgaus.o \
	cgmres.o \
	qerrfun.o \
	gammanew_eval.o \
	csvdpiv.o \
	qleigen_trid.o


LOBJS = $(patsubst %.o,../$(BIN_DIR)/%.o,$(OBJS))
# targets

all: lib

.PHONY : all lib profile release \
  install install-strip uninstall clean distclean setup_dir\
	deepclean

setup_dir: 
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(TMP_DIR)

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

lib: setup_dir $(MODS) $(OBJS) 
	cd $(BIN_DIR); $(FC) $(LDFLAGS) -o $(LIBNAME) $(OBJS)

clean:
	cd $(BIN_DIR); rm -f *
	cd $(MDIR); rm -f *.m $(GATEWAY).c *.mex* *~

deepclean:
	cd $(BIN_DIR); rm -f *
	cd $(TMP_DIR); rm -f *

distclean: clean
	cd $(BIN_DIR); rm -f $(LIBNAME)

printflags: 
	@echo $(DBG) $(OPENMP)

## tests

# note that, without install, we need to point to the library
# at run-time as well... (see the -Wl,-rpath part)
test%: setup_dir lib
	cd $(TMP_DIR); $(FC) $(FFLAGS) -o $@ ../$(TEST_DIR)/$@.f $(PREO) ../$(BIN_DIR)/$(LIBNAME) $(POSTO) -Wl,-rpath=../$(BIN_DIR)
	cd $(TMP_DIR); ./$@

## matlab

$(GATEWAY).c: $(MWRAP_DIR)/$(MWRAPFILE).mw Makefile
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw


# all mex-ing

mexfiles: mexfile mexfile2 mexfile3

# easier just to link against object files than library,
# if it's not being installed (I think)
mexfile: $(GATEWAY).c $(OBJS) Makefile
	cd $(MWRAP_DIR); $(MEX) $(GATEWAY).c $(LOBJS) -largeArrayDims -lgfortran -lmwblas -lgomp -lm

$(GATEWAY2).c: $(MWRAP_DIR)/$(MWRAPFILE2).mw Makefile
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY2) -mb $(MWRAPFILE2).mw
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -mex $(GATEWAY2) -c $(GATEWAY2).c $(MWRAPFILE2).mw


# easier just to link against object files than library,
# if it's not being installed (I think)
mexfile2: $(GATEWAY2).c $(OBJS) Makefile
	cd $(MWRAP_DIR); $(MEX) $(GATEWAY2).c $(LOBJS) -largeArrayDims -lgfortran -lmwblas -lgomp -lm


$(GATEWAY3).c: $(MWRAP_DIR)/$(MWRAPFILE3).mw Makefile
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY3) -mb $(MWRAPFILE3).mw
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -mex $(GATEWAY3) -c $(GATEWAY3).c $(MWRAPFILE3).mw


mexfile3: $(GATEWAY3).c $(OBJS) Makefile
	cd $(MWRAP_DIR); $(MEX) $(GATEWAY3).c $(LOBJS) -largeArrayDims -lgfortran -lmwblas -lgomp -lm
