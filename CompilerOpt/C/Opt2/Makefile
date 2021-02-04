#
# C compiler and options
#
CC=     icc
CFLAGS= -O3 -unroll=0 -no-vec

#
# Object files
#
OBJ=	 unroll.o second_wall.o 

#
# Compile
#
unroll:	$(OBJ)
	$(CC) -o $@ $(OBJ) 

.c.o:
	$(CC) $(CFLAGS) -c $< 

#
# Clean out object files and the executable.
#
clean:
	rm *.o unroll
