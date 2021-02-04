program simplify

integer, parameter :: n=2187, reps=50000
integer :: a
real(kind=8) :: time1,second
 
a=7

time1=second()
do i=1,reps 
   call fred(a,n)
end do 
time1=second()-time1

print *,"a = ",a
print *,"Time =      ",time1
print *,"Iterations per second = ",real(n*reps)/(time1)

end

