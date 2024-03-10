!!!!!!!!!!
!Programa para modelar um levantamento um levantamento gravimétrico sobre um alvo de formato cilín-
!drico em um meio encaixante uniforme.
!Aluno: Lucas de Castro Costa
!Data: 06/07/2021
module Levantamento_gravimetrico
implicit none
Integer,parameter :: pr = kind(1.d0) !variável de precisão das variáveis
Real(pr)          :: R, p, L, x !parametros globais 


Contains
!-----------------------------------
Real(pr) function cilindro(xlinha) !escrevendo a função que será integrada
implicit none
Real(pr), intent(in) :: xlinha         !entrada da função
Real(pr)             :: Termo1, Termo2

Termo1 =((p+L)**2+(x-xlinha)**2)*(sqrt(R**2-xlinha**2)+sqrt(p**2+R**2+ x**2-2*xlinha*x))**2

Termo2 =(p**2+(x-xlinha)**2)*(sqrt(R**2-xlinha**2)+sqrt((p+L)**2+R**2+ x**2-2*xlinha*x))**2

 cilindro=log(Termo1/Termo2)    


End Function
!-----------------------------------
subroutine trapezios(FUN,a,b,S ) !subroutina para calcular as integrais numericas
implicit none
real(pr),external       :: Fun  !função de entrada
real(pr),intent(in)     :: a,b  !limites
real(pr),intent(out)    ::  S   !resultado
real(pr)                :: soma,sex,h,x  
real(pr)                :: S1,S2
real(pr)                :: eps
integer                 :: i,N,k
        
eps = 1.d-6
N = 4
h = (b-a)/N
sex = FUN(a) + FUN(b)
soma = 0.
x = a+h
do i = 1, N-1
	soma = soma + FUN(x)
        x = x + h
end do
S1 = h * (sex + 2*soma)/2.d0
        
k = 0
do 
        x = a - h/2.d0
        do i = 1, N
        	x = x + h
        	soma = soma + FUN(x)
        end do
        h = h / 2.d0
        s2 = h * (sex + 2*soma )/2.d0
        S = (1.d0/3.d0) * (4*S2 - S1)
        if (abs(S-S1) > eps) then
        	N = 2 * N
        	S1 = S
        else
                exit
        end if
        k = k + 1
end do
        
    
end subroutine trapezios!----------------------------------- 
End module !final do modulo
!Inicio do programa
program Gravimetria_cilindro
use Levantamento_gravimetrico
implicit none
Real(pr), Allocatable   :: S(:)     !RESULTADOS DAS INTEGRAIS
Real(pr), Allocatable   :: xlinha(:) !valores de x
Integer                 :: i  !contador
Integer                 :: n  !numeros de pontos

open(11,file='dados.dat',status='old',action='read')  !guardando os valores dos parametros globais, numero de pontos e valores dos pontos
open(22,file='levantamento.dat',status='replace',action='write') !Arquivo para guardar valores do eixo x com relação ao resultado da integral

Read(11,*) R, L, p   !valores dos parametros globais


Read(11,*) n         !numero de pontos

Allocate(xlinha(n),S(n)) !fazendo a alocação dos vetores


Do i=1,n
	Read(11,*) x       
	xlinha(i)= x
	Call trapezios(cilindro,-R,R,S(i))   !aplicando a integração   
End Do    

Do i=1,n
   write(22,*)xlinha(i),S(i) !guardando os valores em um arquivo de saída
End do

End program
