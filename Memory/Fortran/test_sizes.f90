      program test_sizes
      integer i,j,k
      double precision d
      real r
      type fred
      integer a
      integer*2 b
      double precision d
      integer*2 b2
      double precision d2
      end type
      type(fred) var      
      double precision a(30,30)
10    FORMAT(A10,I8,Z16.0,I4)

      write(*,10) 'i',sizeof(i),LOC(i)

      end program
