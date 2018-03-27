
FC = gfortran
FFLAGS = -g

TMP_DIR = ../tmp
BIN_DIR = ../bin
SRC_DIR = ../src

OBJS = hank103.o \
	cdjseval2d.o \
	helmbhrouts.o \
	helmrouts2d.o \
	prini.o

LOBJS = $(patsubst %.o,$(BIN_DIR)/%.o,$(OBJS))

VPATH = $(SRC_DIR):$(BIN_DIR)

%.o: %.f
	$(FC) $(DBGFL) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

%.o: %.f90
	$(FC) $(DBGFL) $(FFLAGS) -c $< -o $(BIN_DIR)/$@

testhelmbhrouts: $(OBJS)
	mkdir -p $(TMP_DIR)
	$(FC) $(FFLAGS) -o $(TMP_DIR)/$@ $@.f $(LOBJS)
	cd $(TMP_DIR); ./$@

testhelmbhders: $(OBJS)
	mkdir -p $(TMP_DIR)
	$(FC) $(FFLAGS) -o $(TMP_DIR)/$@ $@.f $(LOBJS)
	cd $(TMP_DIR); ./$@

clean: 
	cd $(BIN_DIR); rm -f *
	cd $(TMP_DIR); rm -f *
