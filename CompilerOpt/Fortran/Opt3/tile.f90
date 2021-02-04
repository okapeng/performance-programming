program tile

integer, parameter :: n = 729, reps=1
real :: a(n,n), b(n,n), c(n,n)
real(kind=8) :: time1,second
 
a=7.0
b=6.0
c=23.0 

time1=second()
do i=1,reps 
   call matmul(a,b,c,n)
end do 
time1=second()-time1

print *, "c(1,1) = ",c(1,1)
print *,"Time =      ",time1
print *,"Mflops = ",real(2*n*n)*real(n*reps)/(1000000.*time1)

end

