# Makefile for mandelbrot area code

.SUFFIXES: .f90 
#
# Fortran compiler and options
#
FC=	ifort
FFLAGS= -O3 -unroll0 -no-vec
CC=     icc

#
# Object files
#
OBJ=	 unroll.o sum.o second_wall.o 

#
# Compile
#
 unroll:	$(OBJ)
	$(FC) -o $@ $(OBJ) 

.f90.o:
	$(FC) $(FFLAGS) -c $< 

#
# Clean out object files and the executable.
#
clean:
	rm *.o  unroll
