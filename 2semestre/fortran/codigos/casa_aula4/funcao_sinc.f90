!Criar a funcao sinc para todos os pontos utilizando o if
!Aluno: Lucas de Castro Costa
!Data:04/05/2021
!Criar a funcao
real function sinc (x)  !funcao sinc
implicit none
real , intent(in) :: x !variavel

	if ( abs(x) > 1.e-6) then 
		sinc = sin(x)/x
	else
		sinc = 1
	end if

end function sinc !Fim da funcao

!!!!!!Criando o programa

Program testesinc

real :: x
real, external :: sinc

write(*,*) 'Entre com o valor de x:'
read(*,*) x

write(*,*) 'O valor da sinc é', sinc(x)


end program
