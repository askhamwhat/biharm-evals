PROJECT = int2

ifndef SYSTEM
SYSTEM = mac
endif

FC = gfortran -c -w
FFLAGS = -O3

ifeq ($(SYSTEM),mac)
FLINK = gfortran -w -o $(PROJECT) -framework accelerate
endif
ifeq ($(SYSTEM),linux)
FLINK = gfortran -w -o $(PROJECT) 
FEND = -lblas -llapack
endif

EFOL = test7Data
SRC = ../src

.PHONY: all clean list


SOURCES =  test7.f \
  $(SRC)/legeexps.f \
  $(SRC)/prini.f \
  $(SRC)/dlaran.f \
  $(SRC)/corrand.f \
  $(SRC)/hkrand.f \
  $(SRC)/pplot.f \
  $(SRC)/chebexps.f \
  $(SRC)/hank103.f \
  $(SRC)/cdjseval2d.f \
  $(SRC)/helmbhrouts.f \
  $(SRC)/helmrouts2d.f \
  $(SRC)/hbhrouts2d.f \
  $(SRC)/chunks_helmstokes_kerndefs.f90 \
  $(SRC)/chunks_helmstokes_matbuild.f90 \
  $(SRC)/chunks_helmstokes_kerndefs_fark.f90 \
  $(SRC)/chunks_ders.f \
  $(SRC)/chunks.f \
  $(SRC)/chunks_multi.f \
  $(SRC)/chunks_multi_extras.f \
  $(SRC)/cmat_build.f90 \
  $(SRC)/cmat_build_vec.f90 \
  $(SRC)/chunks_quads.f \
  $(SRC)/corners.f \
  $(SRC)/adapgaus.f \
  $(SRC)/cgmres.f \
  $(SRC)/qerrfun.f \
  $(SRC)/gammanew_eval.f \
  $(SRC)/csvdpiv.f \
  $(SRC)/qleigen_trid.f \


OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(PROJECT)
	$(FLINK) $(OBJECTS) $(FEND)
	mv $(PROJECT) $(EFOL)/
	cd $(EFOL) && ./$(PROJECT)

clean:
	rm -f $(OBJECTS)
	rm -f $(PROJECT)

list: $(SOURCES)
	$(warning Requires:  $^)





  
