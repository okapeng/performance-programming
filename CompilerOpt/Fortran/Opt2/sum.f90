subroutine sum(a,b,c,d,n)

integer :: i, n
real :: a(n), b(n), c(n), d
 

do i=1,n 
   a(i) = b(i) + d * c(i) 
end do 

return 
end 
