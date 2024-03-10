!Criar a funcao chapeu utilizando as condicionais com o if
!Aluno: Lucas de Castro Costa
!Data:29/04/2021
!pressupoe que o "a" será escolhido um numero positivo
Real function teste (x, a)  !x e a variavel e o a e constante
implicit none
real, intent(in) :: x, a
!condicoes para a funcao chapeu
  if ( x <= -a) then
  	teste = 0.
  else if ( x <=  0. ) then
	teste = (a + x)/a
  else if ( x <=  a) then
  	teste = (a - x)/a
  else
  	teste= 0.
  end if


end function teste

!!!!!!!!
Program testechapeu
implicit none
real :: x,a
real, external :: teste

a = 1.

write(*,*) 'x?'
read(*,*) x

write(*,*) 'f(x)=', teste (x,a)


end program
