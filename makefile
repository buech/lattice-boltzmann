FC = gfortran
CFLAGS += -std=c99 -Wall -Wextra -pedantic
FFLAGS += -std=f2003 -Wall -Wextra -pedantic -cpp
LDLIBS = -lm

N ?= 128
M ?= 48
T ?= 20000
RE ?= 160.0
ULB ?= 0.04

CPPFLAGS += -DN=$(N) -DM=$(M) -DT=$(T) -DRE=$(RE) -DULB=$(ULB)

OMP ?= 0
ifeq ($(OMP), 1)
	CFLAGS += -fopenmp
	FFLAGS += -fopenmp
endif

all: lb lbf

lb: lb.c

lbf: lb.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) $< -o $@

clean:
	rm -f lb lbf constants.mod

.PHONY: all clean
