#
# C compiler and options
#
CC=     icc  -O3 

#
# Object files
#
OBJ=	 tile.o second_wall.o 

#
# Compile
#
tile:	$(OBJ)
	$(CC) -o $@ $(OBJ) 

.c.o:
	$(CC) -c $< 

#
# Clean out object files and the executable.
#
clean:
	rm *.o tile
