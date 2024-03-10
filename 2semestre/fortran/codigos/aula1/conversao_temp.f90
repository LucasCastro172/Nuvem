!primeiro codigo
Program Fahhreint
!Programa para calcular a conversao de temperatura
!Aluno: Lucas de Castro Costa
!Data: 15/04/2021
Implicit none
!parte de cima e essencial para todos os programas
!declaracao de variaveis
Real:: F !Fahreinht
Real:: C !Celsius
!!!Escrever
Write(*,*) 'Qual a temperatura em F?' !primeiro asterisco para escrever na tela e segundo o formato
!tudo q estiver entre apostrofe será reproduzido no terminal
!!!Leitura
Read(*,*) F !primeiro asterisco para ler e segundo o formato
 C = 5.*(F-32.)/9.  !formula de conversao de celsius para fahrenheit
Write(*,*) 'A temperatura em Celsius e', C !Resposta da conversao em celsius
 



	
	
	
End Program

