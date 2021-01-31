      PROGRAM matrixpower

      IMPLICIT NONE

      INTEGER*8 Nside, Nnum, Two

      PARAMETER(Nside = 512, Nnum=Nside*Nside, Two=2)

      DOUBLE PRECISION mat(Nside,Nside)

      INTEGER*8 flops

      DOUBLE PRECISION mflop

      INTEGER i,j, count

      integer :: start, finish, rate
      real :: time
  
c  this variable will count the total number of floating 
c  point operations executed by the code

      flops = 0

c  read the power to raise the matrix to

      OPEN(12,FILE="input")
      READ(12,*) count

c  matrix to unit matrix
      mat = 0.0
      DO i=1,Nside
       mat(i,i) =1.0
      END DO

      call system_clock(count=start, count_rate=rate)

      CALL matpower(count,mat)

      call system_clock(count=finish)
      time=float(finish-start)/rate

      write(6,10)time
   10 format('Time = ', g11.4)


      mflop = DBLE(flops)/(time*1.0D6)
      WRITE(*,*) 'flop count ', flops
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
      DO j=1,Nside
        DO i=1,Nside
          work(i,j,toggle) = mat(i,j)
        END DO
      END DO

      IF( MOD(n,2) .EQ. 0)THEN
C  set mat to identity
        DO j=1,Nside
          DO i=1,Nside
            IF( i.ne.j)THEN
              mat(i,j) = 0.0D0
            ELSE
              mat(i,j) = 1.0D0
            END IF
          END DO
        END DO
      END IF
      n = n/2           

      DO WHILE( n .gt. 0 )
        WRITE(*,*) n
        IF( MOD(n,2) .EQ. 0)THEN
          CALL square(work(1,1,toggle), work(1,1,(1-toggle)))
        ELSE
          CALL square_mult(work(1,1,toggle), work(1,1,(1-toggle)),mat)
        END IF
        toggle = 1-toggle
        n = n/2
      END DO

      RETURN
      END SUBROUTINE


      Subroutine square( mat1, mat2 )
      IMPLICIT NONE
      DOUBLE PRECISION mat1(Nside,Nside), mat2(Nside,Nside), tmp
      INTEGER i,j,k
C square the matrix
C  NFOPS = Nnum*Nside*2
      DO j=1,Nside
        DO i=1,Nside
          tmp = 0.0D0
          DO k=1,Nside
            tmp = tmp + mat1(i,k) * mat1(k,j)
          END DO
          mat2(i,j) = tmp
        END DO
      END DO      

      flops = flops +  Nnum*Nside*Two

  
      RETURN
      END SUBROUTINE

      Subroutine square_mult( mat1,mat2, mat )
      IMPLICIT NONE
      DOUBLE PRECISION mat1(Nside,Nside), mat2(Nside,Nside), tmp
      DOUBLE PRECISION mat(Nside,Nside)
      INTEGER i,j, k

C square the matrix and multiply last argument by result
C OK to overwrite mat1
C  NFOPS = Nnum*Nside*2*2
      DO j=1,Nside
        DO i=1,Nside
          tmp = 0.0D0
          DO k=1,Nside
            tmp = tmp + mat1(i,k) * mat1(k,j)
          END DO
          mat2(i,j) = tmp
        END DO
      END DO      

      DO j=1,Nside
        DO i=1,Nside
          tmp = 0.0D0
          DO k=1,Nside
            tmp = tmp + mat2(i,k) * mat(k,j)
          END DO
          mat1(i,j) = tmp
        END DO
      END DO      

      DO j=1,Nside
        DO i=1,Nside
          mat(i,j) = mat1(i,j)
        END DO
      END DO
  
      flops = flops +  Nnum*Nside*Two*Two
      RETURN
      END SUBROUTINE
      END



