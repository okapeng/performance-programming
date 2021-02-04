program unroll

integer, parameter :: n = 2187, reps=200000
real :: a(n), b(n), c(n), d
real(kind=8) :: time1,second
 
a=7.0
b=6.0
c=23.0 
d=5.0

time1=second()
do i=1,reps 
   call sum(a,b,c,d,n)
end do 
time1=second()-time1

print *,"a(n) = ",a(n)
print *,"Time =  ",time1
print *,"Mflops = ",real(2*n*reps)/(1000000.*time1)

end

