subroutine matmul(a,b,c,n)

integer :: i, n
real :: a(n,n), b(n,n), c(n,n)
 

do k=1,n 
   do i=1,n
      do j=1,n 
         c(i,j) = c(i,j) + a(i,k) * b(k,j)
      end do
   end do
end do 

return 
end 
