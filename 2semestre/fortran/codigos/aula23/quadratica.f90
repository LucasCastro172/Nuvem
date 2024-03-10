program testecomplex
implicit none
real :: a,b,c
complex :: raiz_1,raiz_2

write(*,*) 'A,B,C:'
read(*,*) A,b, c

call eq_quadratica(a,b,c,raiz_1,raiz_2)

write(*,*) 'Raizes:'
write(*,100) 'R1:', real(raiz_1), ' + ', aimag(raiz_1),'i'
write(*,100) 'R2:', real(raiz_2), ' + ', aimag(raiz_2),'i'

100 format(A4,f10.3,a3,f10.3,a2)

contains

subroutine eq_quadratica(a,b,c,r1,r2)
!ax^2 + bx + c = 0
implicit none
real, intent(in) :: a,b,c
complex, intent(out) :: r1,r2
complex :: raiz_delta,delta

delta = cmplx(b**2 - 4*a*c,0.)
raiz_delta = sqrt(delta)

r1 = (-B +raiz_delta)/ (2*A)

r2 = (-B -raiz_delta)/ (2*A)



end subroutine 

end program
