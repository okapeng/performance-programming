# makefile for matrix power code

#OBJ= matpower.o

FC	= ifort -pg -O0 
CC	= icc

FFLAGS = -c -C -mkl
LDFLAGS = -mkl


.f.o:
	$(FC) $(FFLAGS) $<
.c.o:
	$(CC) $(CFLAGS) $<

all: matpower naive blas_matpower

matpower:  matpower.o
	$(FC) $(LDFLAGS) $(LFLAGS) -o $@ matpower.o $(LIBS)

blas_matpower:  blas_matpower.o
	$(FC) $(LDFLAGS) $(LFLAGS) -o $@ blas_matpower.o $(LIBS)  


naive: matpower_naive.o
	$(FC) $(LDFLAGS) $(LFLAGS) -o $@ matpower_naive.o  $(LIBS)


clean:
	rm -f matpower blas_matpower naive gmon.out *.o *~ Matpower.[oe]*


tar:
	tar cvf matpower.tar matpower.f blas_matpower.f matpower_naive.f README Makefile matpower.pbs input*
