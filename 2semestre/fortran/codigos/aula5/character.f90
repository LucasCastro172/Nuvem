!!!!!!!!!!!!!!! funcao para media
!Usar o character
!Criando funcao
Character(3) function Conceito (nota)
implicit none
real, intent(in) :: nota
!criando as condiconais
!Pressupoe q a nota está entre 0 e 10
if (Nota < 0. ) then
	Conceito = 'ERR'
else if (nota < 5.) then
	Conceito = "INS"
else if (nota < 7.) then
	Conceito = "REG"
else if (nota < 9.) then
	Conceito = "BOM"
else if (nota <= 10.) then
	Conceito = "EXC"
else
	Conceito = "ERR"
end if



end function Conceito

!!!!!!

Program teste_media
implicit none
character(3),external :: Conceito
real :: nota
character(3) :: con

write(*,*) 'Qual a nota?'
read(*,*) nota

con = Conceito(nota)

if (con == 'ERR') then
	write(*,*) 'Nota invalida'
else
	write(*,*) 'O seu conceito é ', con
end if


end program




