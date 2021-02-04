subroutine fred(a,n)

integer :: i, j, n
integer :: a, b, c, d, e
 
i=0
e=4
c=23
do j=0,n
   e = c - j - j
   b = c - j + c - j + c - j
   d = b + e
   i = i + 2
   d = i
   a = a + d + d + d
end do

return 
end 
