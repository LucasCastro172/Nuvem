function sinc(x)
!
! definicao da funcao sinc(x) = sin(pi x ) / ( pi x )
!
implicit none
real,intent(in) :: x
real :: sinc
!
real,parameter::pi=3.14159265358979323846264338328d0

if ( abs(x) > 0.0 ) then
   sinc = sin(pi*x) / (pi*x)
else
   sinc = 1.0
end if

end function sinc

function interp_trace(t,nt,t0,dt,trace)
!
! Interpolacao sinc 
!
implicit none
integer,intent(in) :: nt         ! nt = numero de amostras no arranjo trace(:);
real,intent(in)    :: t,t0,dt    ! t0 = amostra inicial;
                                 ! dt = intervalo de amostragem original;
                                 ! t  = tempo de amostragem  
real,intent(in)    :: trace(nt)  ! arranjo de dados 
real :: interp_trace             ! valor interpolado trace(t)
!
real,external :: sinc
!
real,parameter :: tol = 1.0e-3
integer,parameter :: nlag = 10 
integer :: it,itb,ite
real :: aux0,aux1

aux0 = (t - t0) / dt
it  = floor(aux0)+1
!
if ( abs(aux0-(it-1)) < tol ) then
     interp_trace = trace(it)
end if
!
interp_trace = 0.0
if ( t0 <= t .and. t <= t0+(nt-1)*dt ) then
   itb = max(it-nlag,1)
   ite = min(it+nlag,nt)
   do it=itb,ite
      aux1 = aux0 - real(it-1)
      interp_trace = interp_trace + trace(it)*sinc(aux1)
   end do
end if
end function interp_trace


