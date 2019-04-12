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

EFOL = test5Data
SRC = ../src

.PHONY: all clean list


SOURCES =  test5.f \
  $(SRC)/legeexps.f \
  $(SRC)/prini.f \
  $(SRC)/dlaran.f \
  $(SRC)/corrand.f \
  $(SRC)/hkrand.f \
  $(SRC)/pplot.f \
  $(SRC)/chebexps.f \


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





  
