#
# C compiler and options
#
CC=   icc   -O0 -fno-inline-functions

#
# Object files
#
OBJ=	 simplify.o second_wall.o 

#
# Compile
#
simplify:	$(OBJ)
	$(CC) -o $@ $(OBJ) 

.c.o:
	$(CC) -c $< 

#
# Clean out object files and the executable.
#
clean:
	rm *.o simplify
