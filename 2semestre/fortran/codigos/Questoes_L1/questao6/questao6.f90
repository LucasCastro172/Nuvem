!Programa para criar uma função para o logaritmo
!Aluno: Lucas de Castro Costa
!Data: 14/05/2021
real function log_ax(x, a, b)
implicit none
real, intent(in) :: x, a, b
!if ((x > 0) .and. (a > 0) .and.  (b > 0)) then
	log_ax = (log(x)/log(b))/(log(a)/log(b))
       ! else
       ! write(*,*) 'Não é possível'

!end if
  !log2 = log(x) / log(2.)
end function log_ax

Program Calculo_log
real, external :: log_ax
real  :: x
real :: a, b

write(*,*) 'Entre com o valor do logaritmando'
read(*,*) x

write(*,*) 'Entre com os valores da base'
read(*,*) a,b


write(*,*) 'O valor do logaritmo com logaritmando x com base a é', log_ax(x,a,b)








end program 
!Program llog
!real :: x = 20.
!real :: y, log2

!y = log(x)

!log2 = log(x) / log(2.)

!write(*,*) 'log', y

!write(*,*) 'log2', log2




!end program
