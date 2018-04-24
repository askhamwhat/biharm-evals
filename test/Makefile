
ifndef SYSTEM
SYSTEM = linux
endif

FC = gfortran
FFLAGS = -g

ifeq ($(SYSTEM),linux)
BLAS = -lblas
LAPACK = -llapack
PREO = 
POSTO = $(BLAS) $(LAPACK)
endif

ifeq ($(SYSTEM),mac)
PREO = -framework accelerate
POSTO = 
endif

TMP_DIR = ../tmp
BIN_DIR = ../bin
SRC_DIR = ../src

OBJS = hank103.o \
	cdjseval2d.o \
	helmbhrouts.o hbhrouts2d.o\
	helmrouts2d.o \
	prini.o \
	chunks_helmstokes_kerndefs.o \
	dlaran.o \
	chunks.o \
	hkrand.o \
	legeexps.o \
	pplot.o \
	cmat_build.o cmat_build_vec.o\
	chunks_quads.o

LOBJS = $(patsubst %.o,$(BIN_DIR)/%.o,$(OBJS))

VPATH = $(SRC_DIR):$(BIN_DIR)

%.o: %.f
	$(FC) $(DBGFL) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

%.o: %.f90
	$(FC) $(DBGFL) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

testhelmbhrouts: $(OBJS)
	mkdir -p $(TMP_DIR)
	$(FC) $(FFLAGS) -o $(TMP_DIR)/$@ $@.f $(PREO) $(LOBJS) $(POSTO) 
	cd $(TMP_DIR); ./$@

testhelmbhders: $(OBJS)
	mkdir -p $(TMP_DIR)
	$(FC) $(FFLAGS) -o $(TMP_DIR)/$@ $@.f $(PREO) $(LOBJS) $(POSTO) 
	cd $(TMP_DIR); ./$@

testhbhstokesgreenid: $(OBJS)
	mkdir -p $(TMP_DIR)
	$(FC) $(FFLAGS) -o $(TMP_DIR)/$@ $@.f $(PREO) $(LOBJS) $(POSTO) 
	cd $(TMP_DIR); ./$@

testhbhstokesmatform: $(OBJS)
	mkdir -p $(TMP_DIR)
	$(FC) $(FFLAGS) -o $(TMP_DIR)/$@ $@.f $(PREO) $(LOBJS) $(POSTO) 
	cd $(TMP_DIR); ./$@

clean: 
	cd $(BIN_DIR); rm -f *
	cd $(TMP_DIR); rm -f *
