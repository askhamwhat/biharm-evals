PROJECT = int2

FC = gfortran -c -w
FFLAGS = -O3
FLINK = gfortran -w -o $(PROJECT) -framework accelerate

EFOL = test2Data
SRC = ../src

.PHONY: all clean list


SOURCES =  test2.f \
  $(SRC)/l2dterms.f \
  $(SRC)/laprouts2d.f \
  $(SRC)/legeexps.f \
  $(SRC)/prini.f \
  $(SRC)/cdjseval2d.f \
  $(SRC)/chunks_helmstokes_kerndefs_fark.f90 \
  $(SRC)/dlaran.f \
  $(SRC)/corrand.f \
  $(SRC)/hank103.f \
  $(SRC)/hbhrouts2d.f \
  $(SRC)/helmbhrouts.f \
  $(SRC)/helmrouts2d.f \
  $(SRC)/hkrand.f \
  $(SRC)/cgmres.f \
  $(SRC)/cmat_build.f90 \
  $(SRC)/chunks.f \
  $(SRC)/chunks_multi.f \
  $(SRC)/chunks_multi_extras.f \
  $(SRC)/chunks_ders.f \
  $(SRC)/pplot.f \
  $(SRC)/chunks_quads.f \
  $(SRC)/corners.f \
  $(SRC)/adapgaus.f \
  $(SRC)/gammanew_eval.f \
  $(SRC)/qerrfun.f \


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





  
