FC = gfortran
CFLAGS += -std=c99
FFLAGS += -Wall -Wextra -pedantic -cpp
LDLIBS = -lm

OMP ?= 0
ifeq ($(OMP), 1)
	CFLAGS += -fopenmp
	FFLAGS += -fopenmp
endif

all: lb lbf

lb: lb.c

lbf: lb.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) $< -o $@

clean:
	rm -f lb lbf constants.mod

.PHONY: all clean
