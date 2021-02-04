# Makefile for mandelbrot area code

#
# Fortran compiler and options
#
.SUFFIXES:.f90 
FC=	ifort  -O0 -fno-inline-functions 
CC=     icc

#
# Object files
#
OBJ=	 simplify.o fred.o second_wall.o 

#
# Compile
#
 simplify:	$(OBJ)
	$(FC) -o $@ $(OBJ) 

.f90.o:
	$(FC) -c $< 

#
# Clean out object files and the executable.
#
clean:
	rm *.o  simplify
