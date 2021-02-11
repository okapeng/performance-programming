CC = icc
CFLAGS = -mcmodel=medium -m64 -qno-openmp -O2 -DSTREAM_ARRAY_SIZE=100000000


FF = ifort
FFLAGS = -mcmodel=medium -m64 -qno-openmp -O2 -DSTREAM_ARRAY_SIZE=100000000

all: stream_f.exe stream_c.exe

stream_f.exe: stream.f mysecond.o
	$(CC) $(CFLAGS) -c mysecond.c
	$(FF) $(FFLAGS) -c stream.f
	$(FF) $(FFLAGS) stream.o mysecond.o -o stream_f.exe

stream_c.exe: stream.c
	$(CC) $(CFLAGS) stream.c -o stream_c.exe

clean:
	rm -f stream_f.exe stream_c.exe *.o

# an example of a more complex build line for the Intel icc compiler
