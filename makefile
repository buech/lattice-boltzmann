CFLAGS += -std=c99 -Wall -Wextra -pedantic
FFLAGS += -std=f2003 -Wall -Wextra -pedantic -cpp
LDLIBS += -lhdf5 -lm
FLIBS += -lhdf5_fortran
HDF5_PATH ?= /usr

N ?= 128
M ?= 48
T ?= 20000
INTERVAL ?= 100
RE ?= 160.0
ULB ?= 0.04

CPPFLAGS += -DN=$(N) -DM=$(M) -DT=$(T) -DINTERVAL=$(INTERVAL) -DRE=$(RE) -DULB=$(ULB)

OMP ?= 0
ifeq ($(OMP), 1)
	CFLAGS += -fopenmp
	FFLAGS += -fopenmp
endif

all: lb lbf

lb: lb.c

lbf: lbf.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) -I$(HDF5_PATH)/include $< $(FLIBS) -o $@

clean:
	rm -f lb lbf constants.mod *.o

.PHONY: all clean
