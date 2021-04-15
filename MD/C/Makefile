SRC=MD.c control.c util.c 
OBJ=$(SRC:.c=.o)
CC=icc 
#CFLAGS= -g -O0 -check=uninit -check-pointers:rw -no-vec
#CFLAGS= -O2 -vec -qopt-report=5 -qopt-report-phase=all -align -xCORE-AVX2 -restrict -parallel -qopenmp
CFLAGS= -O2 -vec -qopt-report=5 -qopt-report-phase=all -align -xCORE-AVX2 -restrict -ipo  -qopenmp -parallel


all: MD

MD: $(OBJ)
	$(CC) $(CFLAGS)  -o $@  $(OBJ) -lm


output.dat: MD input.dat
	./MD


clean:
	rm -f MD $(OBJ) 

$(OBJ) : coord.h Makefile


