CFLAGS += -std=c99
CFLAGS += -Wall -Wextra -pedantic
FFLAGS += -Wall -Wextra -pedantic

LDLIBS = -lm

FC = gfortran

all: lb lbf

lb: lb.c

lbf: lb.f90
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -f lb lbf constants.mod

.PHONY: all clean
