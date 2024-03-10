program afd_2d 
!#################################################################### 
! 
! Prototipo para modelagum acustica FD 2D
! 
! Curso: Introducao a modelagem e imageamento com equacao da onda
! Instrutor: jesse costa (UFPA, INCT-GP)
!
! Semana de Inverno de Geofisica
! UNICAMP jul/2012
!
!
!####################################################################
!
use omp_lib
use su_types_mod
use su_io_mod
!
implicit none
!
! precisao dos tipos ponto flutuante
! ==================================
!
! Dupla precisao (dp): calculos.
integer,parameter::sp=kind(1.0)
! Precisao simples (sp): armazanamento em disco.
integer,parameter::dp=kind(1.d0)
!
logical,parameter :: free_surf = .false.
!
! Parametros: 
! =================
! nframes: numero de frames gravados em disco
! nborda : numero de pontos utilizados nas bordas de absorcao
!          Fronteira absorvente aplicada nas bordas da malha
!
integer,parameter  :: nframes = 100
integer,parameter  :: nborda  = 75    !  comp. de onda maximo
real(dp),parameter :: R0pml   = 1.0d-4    ! R_0
real(dp) :: R0log
!
real(dp),allocatable:: gama_x(:), gama_z(:) ! fronteiras de absorcao
real(dp),allocatable:: vel(:,:)           ! modelo de velocidade
real(dp),allocatable:: srcmask(:,:)            ! campo de onda
real(dp),allocatable,target:: pfield1(:,:)            ! campo de onda
real(dp),allocatable,target:: pfield2(:,:)            ! campo de onda
real(dp),pointer:: fswap(:,:)          ! auxiliar pointer
real(dp),pointer:: p1(:,:)            ! campo de onda
real(dp),pointer:: p2(:,:)            ! campo de onda
!
! 
type(su_trace) :: vtrace     ! traco instantaneo do campo de onda
type(su_trace) :: trace
type(su_trace),allocatable :: shotgather(:)  
real(dp) :: dswap
!
! Modelo de velocidade
! ========
! header com geometria da malha e arquivo de dados
! velhdr =     nz,nx,z0,x0,dz,dx 
!              fname.bin 
!
!             | nz,nx dimensao da malha       
!             | z0,x0     
!             | dz,dx intervalo de amostragem
!             | fname.bin: arquivo binario com modelo de velocidade
!
Character(200):: velhdr,veldat
Character(200):: sumovie 
! unidades de arquivo:
integer,parameter :: MDID=72
integer,parameter :: RCID=75
integer,parameter :: SUID=90
integer,parameter :: SUMVID=69
!
!Parametros
!
real(dp),parameter :: pi      = 3.14159265358979323846264338328_dp
real(dp),parameter :: twopi   = 2.0_dp*pi
real(dp),parameter :: cfl     = 0.25d0
real,parameter     :: sec2mcs = 1.0e6
real(dp),parameter :: dtrec   = 0.004d0  ! ms

integer,parameter  :: DRV1LEN = 1
real(dp),parameter :: deriv1(DRV1LEN)=(/1.0d0/)
!
! operadores de diferencas finitas:
!integer,parameter  :: DRVLEN=5
!real(dp),parameter :: deriv2(DRVLEN)=(/-2.84722d0,  &
!                                         1.60000d0,  &
!                                        -2.00000d-1, &
!                                         2.53968d-2, &
!                                        -1.785714d-3/)

integer,parameter  :: DRVLEN=7
real(dp),parameter :: deriv2(DRVLEN)=(/-3.171097979234166, &
                                        1.883703563645937, &
                                       -0.392013173326410, &
                                        0.126419594523541, &
                                       -0.043436978346425, &
                                        0.013633723622292, &
                                       -0.002757740501852 /)

real(dp)::drv2_x(DRVLEN),drv2_z(DRVLEN)

!
!integer,parameter  :: DRVLEN=9
!real(dp),parameter :: deriv2(DRVLEN)=(/-3.182139061273932d0, & 
!                                        1.894741056579879d0, & 
!                                       -0.401774032196622d0, & 
!                                        0.134525605392713d0, & 
!                                       -0.049857793344935d0, & 
!                                        0.017889220775041d0, & 
!                                       -0.005634464140223d0, & 
!                                        0.001370433179210d0, & 
!                                       -0.000190495608096d0 /)

real(dp) :: zf,xf,zr,xr,alpha,beta,delt,gama,invpgama,mgama,laplacian,fonte
real(dp) :: scalx,scalz,scalt
real(dp) :: ttotal,t,vmin,vmax,dxmax
real(sp) :: swap
real(dp) :: x0,z0,dx,dz,dt,freq
real(dp) :: fdtbeg,fdtend,fd_update_time
real(dp) :: swptbeg,swptend,swp_time
integer  :: nx,nz,nt,it 
integer  :: ixf,izf,ixr,izr,ix,iz,it0,iconv
integer  :: ixb,ixe,izb,ize,nxx,nzz,ndt
integer  :: ierr,nsnap,iframe,irec,itrc,itrcmv,idir,ntrc,ns
integer  :: chunksize,nthreads,tid
!
!
! entrada dos parametros da modelagem
! =============================================================
!
! entrada do modelo de velcidades 
!
write(*,'("Entre vel. model header file: ")')
read(*,'(A)') velhdr
write(*,'("Entre com a posicao da fonte (zs,xs) : ")')
read(*,*) zf,xf
write(*,'("Entre com a posicao do geofone (zs,xs) : ")')
read(*,*) zr,xr
write(*,'("Entre freq. pico do pulso fonte: ")')
read(*,*) freq
write(*,'("Entre tempo total de modelagem: ")')
read(*,*) ttotal
write(*,'("Entre arquivo de saida (.su): ")')
read(*,*) sumovie
!
write(*,'("Entre vel. model header file: ",a)') trim(velhdr) 
write(*,'("Entre com a posicao da fonte (zs,xs) : ",f12.6,f12.6)') zf,xf
write(*,'("Entre com a posicao do geofone (zs,xs) : ",f12.6,f12.6)') zr,xr
write(*,'("Entre freq. pico do pulso fonte: ",f12.6)') freq
write(*,'("Entre tempo total de modelagem: ",f12.6)')ttotal
write(*,'("Entre arquivo de saida (.su): ",a)') trim(sumovie)
!
open(30, file=trim(velhdr),iostat=ierr) 
if ( ierr /= 0 ) then
     write(*,'("Erro ao abrir arquivo: ",A)') trim(velhdr)
     stop 
end if
!
read(30,*) nz,nx,z0,x0,dz,dx
read(30,'(a)') veldat
close(30)

! condicoes de absorcao nos 4 lados da malha:
   nxx=nx+2*nborda
   ixb=nborda+1
   ixe=nborda+nx
   nzz=nz+2*nborda
   izb=nborda+1
   ize=nborda+nz
!
allocate ( vel(nzz,nxx),pfield1(nzz,nxx),pfield2(nzz,nxx),        &   ! alocacao de memoria na malha incluindo bordas de absorcao
           srcmask(nzz,nxx),gama_x(nxx),gama_z(nzz), &   ! 
           stat=ierr)


!
if (ierr/=0) then                      ! informa falha na alocacao de memoria
    stop 'falha de alocaçao da malha'
end if

p1=>pfield1(:,:)
p2=>pfield2(:,:)


allocate(vtrace%trcdata(nz),stat = ierr)
if (ierr/=0) then                      
        stop 'falha de alocaçao do traco'
end if
!
! header do traco para instataneo do campo de onda
vtrace%trcheader%ns=nz
vtrace%trcheader%ntr=nx
vtrace%trcheader%d1=dz
vtrace%trcheader%d2=dx
vtrace%trcheader%f1=0.0
vtrace%trcheader%f2=0.0

! acessa o arquivo de velocidades
idir=max(index(velhdr,'\',.true.), index(velhdr,'/',.true.))
veldat=velhdr(1:idir)//trim(veldat) ! adiciona caminho do diretorio
open(MDID,file=trim(veldat), form='unformatted', access='direct', recl= SP, iostat=ierr)
if ( ierr /= 0 ) then
     write(*,'("Erro ao abrir arquivo: ",A)') trim(veldat)
     stop 
end if
!
  irec = 0
  ! leitura do modelo de velocidade
  do ix=ixb,ixe
     do iz=izb,ize
        irec = irec + 1
        read(MDID,rec=irec) swap
        vel(iz,ix) = swap 
     end do
  end do
close(MDID) !  
!
! calcula os vmax e vmin 
vmax=vel(izb,ixb)                      !  
vmin=vel(izb,ixb)                      !
do ix=ixb,ixe                          ! 
   do iz=izb,ize                       ! 
      vmin=min(vmin,vel(iz,ix))        !   
      vmax=max(vmax,vel(iz,ix))        !  
   end do
end do
!
! preenche o modelo de velocidades nas bordas 
!
do ix=1,nborda
   do iz=izb,ize
      vel(iz,ix)     = vel(iz,ixb)
      vel(iz,ixe+ix) = vel(iz,ixe)
   end do
end do

do ix=1,nxx
   do iz=1,izb-1
      vel(iz,ix)=vel(izb,ix)
   end do
   do iz=ize+1,nzz
      vel(iz,ix)=vel(ize,ix)
   end do
end do
!
! posicao da fonte na malha
! 
ixf=int((xf-x0)/dx)+ixb
izf=int((zf-z0)/dx)+izb
!
! posicao da fonte na malha
! 
ixr=int((xr-x0)/dx)+ixb
izr=int((zr-z0)/dx)+izb

!
!srcmask(:,:)     = 0.d0
!srcmask(izf,ixf) = 1.d0
!
ns     = floor(ttotal/dtrec)+1
!
! verifica condicoes de  estabilidade e dispersao numerica
!
dxmax = vmin / (6.0_dp*freq)
dt    = cfl *sqrt(2.0/sum(abs(deriv2(:))))*min(dx,dz)/ vmax
!dt    = cfl * dx   / vmax
ndt   = floor(dtrec/dt+0.5) + 1      
!
! intervalo de tempo de registro multiplo inteiro
! do intervalo de modelagem  dt
!
if ( ndt > 1) then
   dt = dtrec/real(ndt)
else
   dt = dtrec
end if

nt =ceiling(ttotal/dt)          
if( max(dx,dz) > dxmax ) then
   stop 'freq muito alta'
end if
!
! velocidade de propagacao normalizada: 
!
do ix=1,nxx
   do iz=1,nzz
      vel(iz,ix)=vel(iz,ix)**2
   end do
end do
scalx=(dt/dx)**2
scalz=(dt/dz)**2
scalt=dt*dt/(dx*dz)

!drv2_z(:)=scalz*deriv2(:)
!drv2_x(:)=scalx*deriv2(:)


! mascara de absorcao:
gama_x(:) = 0.0_dp
gama_z(:) = 0.0_dp

R0log = log(R0pml) !/log(10.d0)
alpha = -1.5d0*R0log*vmax/(real(Nborda,dp)*dx)
print *,'alpha = ',alpha
alpha = dt*alpha

beta= pi*freq*dt
   do ix=1,Nborda
      gama = abs(real(ix,dp)/real(Nborda,dp))
      gama = alpha*(gama-sin(twopi*gama)/twopi)
      gama_x(ixb-ix) = gama
      gama_x(ixe+ix) = gama
      gama_z(izb-ix) = gama
      gama_z(ize+ix) = gama
   end do
!
! inervalo entre instantaneos do campo de onda
nsnap  = max(1,nt/nframes)

!
! distribuicao espacial da fonte
!
!$omp parallel                                & 
!$omp shared(nzz,nxx,srcmask,izf,ixf) &
!$omp private(ix,iz)
!$omp do schedule(static) collapse(2)
do ix=1,nxx
   do iz=1,nzz
      srcmask(iz,ix) = 4.d0*scalt*exp(-2.0_dp*((real(iz-izf,dp)**2+real(ix-ixf,dp)**2)))/twopi
   end do
end do
!$omp end do 
!$omp end parallel 

print *, '# of time steps  = ', nt
print *, 'time interval dt = ',dt
!
! Instataneous do campo de onda sao armazenados
! no arquivo fdmovie.su
!
ierr=OpenSU(SUMVID,trim(sumovie),'w',nz)
if ( ierr /= 0 ) then
     write(*,'("Erro ao abrir arquivo: ",A)') trim(sumovie)
     stop 
end if
!
open(RCID,file='trace_fd.dat')
!
!condiçao inicial: 
! p1(:,:) pressure field at time t-dt
! p2(:,:,) pressure field at time t
! set to zero
!$omp workshare
p1(:,:) = 0.0d0
!$omp end workshare

!$omp workshare
p2(:,:) = 0.0d0
!$omp end workshare


!
! medidas de tempo usado omp_lib
fd_update_time=0.0_dp
swp_time=0.0_dp
!
iframe = 0
itrcmv = 0
!
! pulso fonte causal 
!
it0 = floor(sqrt(2.d0)/(freq*dt))+1 ! pico pulso Ricker            
!
!
!$omp parallel default(none)                                                     &
!$omp shared(nzz,nxx,nz,p1,p2,vel,fswap,gama_x,gama_z,fonte,srcmask)             &
!$omp shared(fd_update_time,swp_time,nt,it0,dx,ixb,ixe,izb,ize,vtrace)           &
!$omp shared(nsnap,iframe,itrcmv,beta,ndt,izr,ixr,scalx,scalz,drv2_z,drv2_x)     &
!$omp private(ierr,iz,ix,iconv)                                                  &
!$omp private(gama,invpgama,mgama,dswap,laplacian,fdtbeg,fdtend,swptbeg,swptend) &
!$omp private (it,delt,tid)
loop_t :do it=1,nt                    ! Loop de evolucao do campo de pressao
!
! Pulso Ricker:
        delt  = (beta*real(it-it0,dp))**2
        fonte = exp(-(delt))*(1.0_dp-2.0_dp*delt) ! pulso Ricker

!
!******************************************************************************************************************
! update PML auxiliary fields
fdtbeg=omp_get_wtime()
!******************************************************************************************************************
! atualizacao do campo de onda:
!
!$omp do schedule(static) collapse(2)

!
! localidade: permutar a ordem dos loops
do ix=DRVLEN,nxx-DRVLEN
   do iz=DRVLEN,nzz-DRVLEN

   gama      = gama_x(ix)+gama_z(iz)
   !invpgama  = 1.0_dp/(1.0_dp+gama)
   !mgama     = 1.0_dp-gama
   !
   ! calculo do laplaciano:
   laplacian= (scalx+scalz)*deriv2(1)*p2(iz,ix)
   do iconv=1,DRVLEN-1
      laplacian =  laplacian +  & 
                   deriv2(iconv+1)*(scalz*(p2(iz-iconv,ix)+p2(iz+iconv,ix)) + &
                                    scalx*(p2(iz,ix-iconv)+p2(iz,ix+iconv)) ) 
   end do
   !   
   !laplacian= (drv2_z(1)+drv2_x(1))*p2(iz,ix)
   !do iconv=1,DRVLEN-1
   !   laplacian =  laplacian +  & 
   !                drv2_z(iconv+1)*(p2(iz-iconv,ix)+p2(iz+iconv,ix)) + &
   !                drv2_x(iconv+1)*(p2(iz,ix-iconv)+p2(iz,ix+iconv)) 
   !end do


   ! FD update:
   !p1(iz,ix) = invpgama*(                                                                              & 
   !                       2.d0*p2(iz,ix)-mgama*p1(iz,ix) +           & 
   !                       vel(iz,ix)*( laplacian + fonte * srcmask(iz,ix) ) & 
   !                     )
   p1(iz,ix) = (2.d0-gama)*p2(iz,ix)-(1.d0-gama)*p1(iz,ix) +           & 
                          vel(iz,ix)*( laplacian + fonte * srcmask(iz,ix) )  
                        


    end do
end do
!$omp end do
fdtend=omp_get_wtime()
fd_update_time=fd_update_time+(fdtend-fdtbeg)
!

!$omp single
swptbeg=omp_get_wtime()
!
! permutacao de arranjos no time:
fswap=>p2
p2=>p1
p1=>fswap
swptend=omp_get_wtime()
swp_time=swp_time+(swptend-swptbeg)
!  
! armazena instantaneos do campo de onda
if (mod(it-1,nsnap)==0.0) then
    iframe=iframe+1
    if ( mod(iframe,10) == 0 ) write(*,'("writing snapshot # ",i5.5)') iframe
    do ix=ixb,ixe
       itrcmv=itrcmv+1
       vtrace%trcheader%tracl=itrcmv
       vtrace%trcdata(1:nz) = real(p2(izb:ize,ix),sp)
       ierr=PutSUTrace(SUMVID,itrcmv,vtrace)
       flush(SUMVID)
    end do
end if
!$omp end single 

end do loop_t ! fim da evolucao
!$omp end parallel
close(RCID)
close(SUMVID)
!
print *,"fd_time = ",fd_update_time
print *,"swap_time = ",swp_time
!
! free allocated memory
deallocate(pfield1,pfield2,vel,srcmask,gama_x,gama_z)
!                     
end program afd_2d
