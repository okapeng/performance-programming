      PROGRAM matrixpower

      IMPLICIT NONE

      INTEGER*8 Nside, Nnum,Two

      PARAMETER(Nside = 512, Nnum=Nside*Nside,Two=2)

      DOUBLE PRECISION mat(Nside,Nside)

      INTEGER(KIND=8) flops


      DOUBLE PRECISION mflop

      INTEGER i,j, count

      integer :: start, finish, rate
      real :: time

c  this variable will count the total number of 
c  floating point ops executed by the code

      flops = 0

c  read the power to raise the matrix to

      OPEN(12,FILE="input")
      READ(12,*) count

c  matrix to unit matrix
      mat = 0.0
      DO i=1,Nside
       mat(i,i) =1.0
      END DO

c      start = etime(tarray)

      call system_clock(count=start, count_rate=rate)

      CALL matpower(count,mat)

      call system_clock(count=finish)
      time=float(finish-start)/rate

      write(6,10)time
   10 format('Time = ', g11.4)

      mflop = DBLE(flops)/(time*1.0D6)

      WRITE(unit=6, fmt=*) 'Flop count is ', flops
      WRITE(unit=6, fmt=*) 'Flop rate is ', mflop, ' Mflops'

c write non-zero elements

      DO j=1,Nside
        DO i=1,Nside
          IF( mat(i,j).ne.0.0D0) THEN
            WRITE(*,'(2I4,F16.8)') i,j,mat(i,j)
          END IF
        END DO
      END DO

      CONTAINS

C
C This subroutine raises a matrix to a +ve power
C
      Subroutine matpower(n, mat)
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION mat(Nside,Nside), work(Nside,Nside,0:1)
      INTEGER toggle
      INTEGER i,j

      IF ( n .le. 0 ) STOP

      IF( n .eq. 1 ) RETURN

      toggle=0
C  start with the identity
        DO j=1,Nside
          DO i=1,Nside
            IF( i.ne.j)THEN
              work(i,j,toggle) = 0.0D0
            ELSE
              work(i,j,toggle) = 1.0D0
            END IF
          END DO
        END DO

C multiply by mat n times
      DO i=1,n
        WRITE(*,*) i
        CALL prod(work(1,1,(1-toggle)),mat,work(1,1,toggle))
        toggle= 1-toggle
      END DO

C copy the final result back to mat
      DO j=1,Nside
        DO i=1,Nside
          mat(i,j) = work(i,j,toggle)
        END DO
      END DO

      RETURN
      END Subroutine


      Subroutine prod( res, mat1, mat2 )
      IMPLICIT NONE
      DOUBLE PRECISION res(Nside,Nside)
      DOUBLE PRECISION mat1(Nside,Nside), mat2(Nside,Nside), tmp
      INTEGER i,j,k
C multiply the matrices
C  NFOPS = Nnum*Nside*2
      DO j=1,Nside
        DO i=1,Nside
          tmp = 0.0D0
          DO k=1,Nside
            tmp = tmp + mat1(i,k) * mat2(k,j)
          END DO
          res(i,j) = tmp
        END DO
      END DO      

      flops = flops +  Nnum*Nside*Two
  
      RETURN
      END Subroutine

      END


