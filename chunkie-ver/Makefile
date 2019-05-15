# some system dependent settings...

ifndef SYSTEM
SYSTEM = linux
endif

ifeq ($(SYSTEM),linux)
DBG=
OPENMP= --openmp
FFLAGS=$(DBG) -O3 $(OPENMP) -fPIC
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
SHELL = /bin/sh
OBJSUF=o
MODSUF=mod
FC = gfortran
MEX=mex
MWRAP=mwrap
MWFLAGS=-c99complex
endif

# locations

SRC_DIR = fsrc
BIN_DIR = bin
MWRAP_DIR = mwrap

VPATH = $(SRC_DIR):$(BIN_DIR)

# names

LIBBASE = obh2d
GATEWAY = $(LIBBASE)gateway
MWRAPFILE = obh2d


OBJS = hank103.o \
	cdjseval2d.o \
	obhgreen2d.o ostokgreen2d.o \
	helmrouts2d.o 

LOBJS = $(patsubst %.o,$(BIN_DIR)/%.o,$(OBJS))
# targets

all: matlab

.PHONY : all lib profile release \
  install install-strip uninstall clean distclean \
	deepclean

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

clean:
	cd $(BIN_DIR); rm -f *

deepclean:
	cd $(BIN_DIR); rm -f *

distclean: clean
	cd $(BIN_DIR); rm -f $(LIBNAME)

printflags: 
	@echo $(DBG) $(OPENMP)

## matlab

# all mex-ing

# easier just to link against object files than library,
# if it's not being installed (I think)

$(GATEWAY).c: $(MWRAP_DIR)/$(MWRAPFILE).mw Makefile
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAP_DIR)/$(MWRAPFILE).mw
	cd $(MWRAP_DIR); $(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw

matlab: $(GATEWAY).c $(OBJS) Makefile
	$(MEX) $(MWRAP_DIR)/$(GATEWAY).c $(LOBJS) -largeArrayDims -lgfortran -lmwblas -lgomp -lm
