# Makefile for mandelbrot area code

#
# Fortran compiler and options
#
.SUFFIXES: .f90 

FC=	ifort  -O3 -qno-opt-matmul
CC=     icc

#
# Object files
#
OBJ=	 tile.o matmul.o second_wall.o 

#
# Compile
#
 tile:	$(OBJ)
	$(FC) -o $@ $(OBJ) 

.f90.o:
	$(FC) -c $< 

#
# Clean out object files and the executable.
#
clean:
	rm *.o  tile
