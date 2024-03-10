!!!!!!!!!!
!Programa para calcular o ponto de impacto no solo após um lançamento de um projetil
!Aluno: Lucas de Castro Costa
!Data: 25/05/2021
program Primeira_avaliacao
implicit none
real :: v0, angulo, h0 !Entrada do programa: angulo, velocidade inicial e altura inicial. Unidades -> Velocidade inicial: m/s e 
!altura: metros
real :: t1, t2 !Saída do programa: raízes da equação de segundo grau que determinam o tempo de impacto. Unidade: segundos
integer :: nt !Termo para definir qual equação será utilizada da subroutina
real :: x ! Saída do programa: distância até o ponto de impacto. Unidade: metros
real, parameter :: pi = 3.14159265359 !Definição de parametro: pi

!deslocamento =  distancia percorrida

!Entradas do programa
write(*,*) 'Entre com a velocidade inicial(m/s):'
read(*,*) v0 

write(*,*) 'Entre com o angulo em graus:'
read(*,*) angulo

write(*,*) 'Entre com a altura inicial(m):'
read(*,*) h0

if (h0  < 0.) then !garantir que não pode colocar uma altura negativa.
	write(*,*) 'Erro. Entre com uma altura não negativa.'
	stop
end if

if (v0 < 0.) then !garantir que a velocidade será positiva. 
	write(*,*) 'Erro. Entre com uma velocidade inicial não negativa.'
	stop
end if

call lancamento_projetil( v0, angulo, h0, nt, t1, t2) !Chamando a subroutina da resolução da equação de 2 grau

if (nt == 2) then
	if (t1 > 0.) then !limitação para valores acima de 0 para não gerar distancias com valor negativo
		if (angulo < 90. .and. angulo > -90.) then
			write(*,*) 't1:', t1, 'segundo(s)' !lendo a saída
			x = v0*cos(angulo * (pi/180))*t1 !Calculando a distancia até o ponto de impacto
			Write(*,*) 'O deslocamento é de', x, 'metro(s).' !lendo a saída
		else if (angulo == 90. .or. angulo ==  -90.) then
			write(*,*) 't1:', t1, 'segundo(s)' !lendo a saída
			Write(*,*) 'A distância percorrida é nula.' !lendo a saída
		end if		
	else  
		stop
	end if
	
	if (t2 > 0.) then
		if (angulo < 90. .and. angulo > -90.) then
			write(*,*) 't2:', t2, 'segundo(s)' !lendo a saída
			x = v0*cos(angulo * (pi/180))*t2 !Calculando a distancia até o ponto de impacto
			Write(*,*) 'O deslocamento é de', x, 'metro(s).' !lendo a saída
		else if (angulo == 90. .or. angulo ==  -90.) then
			write(*,*) 't2:', t2, 'segundo(s)'!lendo a saída
			Write(*,*) 'A distância percorrida é nula.' !lendo a saída
		end if	
	else 
		stop
	end if

	
else if (nt == 1) then
		if (t1 > 0.) then
			if (angulo < 90. .and. angulo > -90.) then
				write(*,*) 't1:', t1, 'segundo(s)'!lendo a saída
				x = v0*cos(angulo * (pi/180))*t1 !Calculando a distancia até o ponto de impacto
				Write(*,*) 'O deslocamento é de', x, 'metro(s).' !lendo a saída
			else if (angulo == 90. .or. angulo ==  -90.) then
				write(*,*) 't1:', t1, 'segundo(s)' !lendo a saída
				Write(*,*) 'A distância percorrida é nula.' !lendo a saída
			end if	
else 
			write(*,*) 't1:', t1, 'segundo(s)' !lendo a saída
			write(*,*) 'Esse tempo não pode ser utilizado para determinar o deslocamento.'
		end if
	
end if



CONTAINS

!Criando a subroutina

subroutine lancamento_projetil( v0, angulo, h0, nt, t1, t2) !subroutina para solucionar a equaçao de 2 grau
!calculo das raizes da equacao gt^2 - 2v0sin(angulo)*t -2 h0 =0
implicit none
real, intent(in) :: v0, angulo, h0 !Valores de entrada
real :: graus !fator que está transformado o angulo de radianos em graus
integer, intent(out) :: nt !termo para chamar no programa, definindo o formato das raízes
real, parameter :: g = 9.8 !gravidade
real, parameter :: pi = 3.14159265359 ! pi 
real, intent(out) ::  t1, t2 !raizes
real :: delta !termo que guarda o calculo do delta

if (angulo <= 90 .and. angulo >= -90) then !limitando a solução entre -90 e 90 graus
	graus =  angulo * (pi/180)
	Delta = (-2.*v0*sin(graus))**2 - 4.*g*(-2.*h0)

	if (Delta > 0.) then !resolução com duas raizes reais e diferentes
		nt = 2
		t1 =  (2.*v0*sin(graus) + sqrt(delta))/(2.*g)
		t2 =  (2.*v0*sin(graus) - sqrt(delta))/(2.*g)
	else if (delta == 0.) then !resolução com duas raizes reais e iguais
		nt = 1
		t1 = (v0*sin(graus))/g
		t2 = t1
	end if
else 
	write(*,*) 'Utilize um angulo entre -90 e 90 graus' !erro para angulos fora desse intervalo
	
end if


end subroutine lancamento_projetil !finalizando a subroutina
	
!Perguntas
!Sua rotina irá funcionar para ângulos menores que zero? Se sim, também
!haverá alguma restrição de ângulo para este caso?

!Sim, está funcionando. A restrição ocorre para valores menores que -90.
!Isto ocorre pois o cosseno é positivo no primeiro e no quarto quadrante, fazendo com que gere um deslocamento
!positivo entre -90 e 90. Nos angulos -90 e 90 devem dar distancia percorrida igual a 0, então isto foi definido 
!no programa, pois, por problema de precisão, o cosseno de -90 e de 90 davam valores muito proximos de 0, só que
!negativos.



end program	
