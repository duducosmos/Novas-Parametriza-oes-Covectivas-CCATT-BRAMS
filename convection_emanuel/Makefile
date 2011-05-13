# ============================================================================
# Name        : Makefile
# Author      : Eduardo S. Pereira
# Version     : 1.0
# Copyright   : GLP
# Description : Makefile for Hello World in Fortran
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler
FORTRAN_COMPILER = gfortran

all: src/convect1.f90
	$(FORTRAN_COMPILER) -O2 -g \
		-o bin/Convec \
		src/Convec.f90

clean:
	rm -f bin/Convec *.mod
