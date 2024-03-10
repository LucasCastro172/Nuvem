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
logical,parameter :: free_surf = .true.
!
! Parametros: 
! =================
! nframes: numero de frames gravados em disco
! nborda : numero de pontos utilizados nas bordas de absorcao
!          Fronteira absorvente aplicada nas bordas da malha
!
integer,parameter  :: nframes = 100
integer,parameter  :: nborda  = 75    !  comp. de onda maximo 
!
real(dp),allocatable:: vel(:,:)             ! modelo de velocidade
real(dp),allocatable:: p(:,:,:)             ! campo de onda
real(dp),allocatable:: sourcemask(:,:)      ! fonte
real(dp),allocatable:: gama_x(:), gama_z(:) ! fronteiras de absorcao
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
real(dp),parameter :: cfl     = 0.25d0
real,parameter     :: sec2mcs = 1.0e6
real(dp),parameter :: dtrec   = 0.004d0  ! ms
!
! operadores de diferencas finitas:
!integer,parameter  :: DRVLEN=5
!real(dp),parameter :: deriv2(DRVLEN)=(/-2.84722d0,  &
!                                         1.60000d0,  &
!                                        -2.00000d-1, &
!                                         2.53968d-2, &
!                                        -1.785714d-3/)

!integer,parameter  :: DRVLEN=6
!real(dp),parameter :: deriv2(DRVLEN)=(/-2.643735850583032, &
!                                        1.457463853662507, &
!                                       -0.149099855401351, &
!                                        0.008388073784959, &
!                                        0.007161491917418, &
!                                       -0.002045638672017 /)

!integer,parameter  :: DRVLEN=7
!real(dp),parameter :: deriv2(DRVLEN)=(/-3.171097979234166, &
!                                        1.883703563645937, &
!                                       -0.392013173326410, &
!                                        0.126419594523541, &
!                                       -0.043436978346425, &
!                                        0.013633723622292, &
!                                       -0.002757740501852 /)

!integer,parameter  :: DRVLEN=9
!real(dp),parameter :: deriv2(DRVLEN)=(/-3.197665217827583d0, &
!                                       1.909538627823793d0, &
!                                      -0.414569081291830d0, &
!                                       0.144471891998828d0, &
!                                      -0.056701894105249d0, &
!                                       0.021948574612036d0, &
!                                      -0.007605843676202d0, &
!                                       0.002080291018278d0, &
!                                      -0.000329957465863d0 /)
!
integer,parameter  :: DRVLEN=9
real(dp),parameter :: deriv2(DRVLEN)=(/-3.182139061273932d0, & 
                                        1.894741056579879d0, & 
                                       -0.401774032196622d0, & 
                                        0.134525605392713d0, & 
                                       -0.049857793344935d0, & 
                                        0.017889220775041d0, & 
                                       -0.005634464140223d0, & 
                                        0.001370433179210d0, & 
                                       -0.000190495608096d0 /)

logical, parameter :: homog   = .false.
real(dp),parameter :: c_homog = 3000.d0

real(dp) :: zf,xf,zr,xr,beta,delt,gama,invpgama,mgama,laplacian,fonte
real(sp) :: ttotal,t,vmin,vmax,dxmax,swap
real(dp) :: x0,z0,dx,dz,dt,freq
integer  :: nx,nz,nt,it 
integer  :: ixf,izf,ixr,izr,ix,iz,it0,iconv
integer  :: ixb,ixe,izb,ize,nxx,nzz,ndt
integer  :: ierr,nsnap,iframe,irec,itrc,itrcmv,idir,ntrc,ns
!
!
! entrada dos parametros da modelagem
! =============================================================
!
! entrada do modelo de velcidades 
!
write(*,'("Entre vel. model header file: ",$)')
read(*,'(A)') velhdr
write(*,'("Entre com a posicao da fonte (zs,xs) : ",$)')
read(*,*) zf,xf
write(*,'("Entre com a posicao do geofone (zs,xs) : ",$)')
read(*,*) zr,xr
write(*,'("Entre freq. pico do pulso fonte: ",$)')
read(*,*) freq
write(*,'("Entre tempo total de modelagem: ",$)')
read(*,*) ttotal
write(*,'("Entre arquivo de saida (.su): ",$)')
read(*,*) sumovie
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
   if ( free_surf ) then
   nzz=nz+nborda
   izb=1
   ize=nz
   else
   nzz=nz+2*nborda
   izb=nborda+1
   ize=nborda+nz
   end if
! 
allocate (p(nzz,nxx,2),            &   ! alocacao de memoria na malha incluindo bordas de absorcao
          vel(nzz,nxx),            &   ! 
          sourcemask(nzz,nxx),     &   !  
          gama_x(nxx),gama_z(nzz), &   ! 
          stat=ierr)
!
if (ierr/=0) then                      ! informa falha na alocacao de memoria
    stop 'falha de alocaçao da malha'
end if

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
!  teste meio homogeneo 
if ( homog ) then
  vel(:,:) = c_homog
end if
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
! distribuicao espacial da fonte
!
!do ix=1,nxx
!   do iz=1,nzz
!      sourcemask(iz,ix) = exp(-4.0_dp*((real(iz-izf,dp)**2+real(ix-ixf,dp)**2)))
!   end do
!end do
!
sourcemask(:,:)     = 0.d0
sourcemask(izf,ixf) = 1.d0
!
!
!
ns     = floor(ttotal/dtrec)+1
!
! verifica condicoes de  estabilidade e dispersao numerica
!
dxmax = vmin / (6.0_dp*freq)
dt    = cfl *sqrt(2.0/sum(abs(deriv2(:)))) * dx   / vmax
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
if( dx > dxmax ) then
   stop 'freq muito alta'
end if
!
! velocidade de propagacao normalizada: 
!
do ix=1,nxx
   do iz=1,nzz
      vel(iz,ix)=(vel(iz,ix)*dt/dx)**2
   end do
end do
! mascara de absorcao:
gama_x(:) = 0.0_dp
gama_z(:) = 0.0_dp

beta= pi*freq*dt
if ( free_surf ) then
   do ix=1,Nborda
      gama = beta*(real(ix,dp)/real(Nborda,dp))**2
      gama_x(ixb-ix) = gama
      gama_x(ixe+ix) = gama
      gama_z(ize+ix) = gama
   end do
else
   do ix=1,Nborda
      gama = beta*(real(ix,dp)/real(Nborda,dp))**2
      gama_x(ixb-ix) = gama
      gama_x(ixe+ix) = gama
      gama_z(izb-ix) = gama
      gama_z(ize+ix) = gama
   end do
end if
!
! inervalo entre instantaneos do campo de onda
nsnap  = max(1,nt/nframes)

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
! p(:,:,1) pressure field at time t-dt
! p(:,:,2) pressure field at time t
!
p(:,:,:)= 0.0_dp
!
iframe = 0
itrcmv = 0
!
! pulso fonte causal 
!
it0 = floor(sqrt(3.d0/2.d0)/(freq*dt)+0.5)+1 ! pico pulso Ricker            
!
if ( free_surf ) then
loop_t_free_surf :do it=1,nt                    ! Loop de evolucao do campo de pressao
!
! Pulso Ricker:
        delt  = (beta*real(it-it0,dp))**2
        fonte = (dx*dx)*exp(-(delt))*(1.0_dp-2.0_dp*delt) ! pulso Ricker
!        
!FD scheme:

!

do ix=DRVLEN,nxx-DRVLEN
   do iz=2,DRVLEN

      gama      = gama_x(ix)+gama_z(iz)
      invpgama  = 1.0_dp/(1.0_dp+gama)
      mgama     = 1.0_dp-gama
      !
      ! calculo do laplaciano:
      laplacian= 2.d0*deriv2(1)*p(iz,ix,2)
      do iconv=1,DRVLEN-1
         laplacian = laplacian +       & 
                     deriv2(iconv+1)*( &
                                       sign(1.d0,real(iz-iconv,dp))*p(max(iz-iconv,iconv),ix,2)+p(iz+iconv,ix,2) + &
                                       p(iz,ix-iconv,2)+p(iz,ix+iconv,2)   &
                                     ) 
      end do
      !
      ! atualizacao do campo de pressao:
      p(iz,ix,1) = invpgama * (                                                      & 
                                2.0_dp*p(iz,ix,2)-mgama*p(iz,ix,1)   +               & 
                                vel(iz,ix)*( laplacian + fonte * sourcemask(iz,ix) ) &
                              )
    end do

    do iz=DRVLEN+1,nzz-DRVLEN

      gama      = gama_x(ix)+gama_z(iz)
      invpgama  = 1.0_dp/(1.0_dp+gama)
      mgama     = 1.0_dp-gama
      !
      ! calculo do laplaciano:
      laplacian= 2.d0*deriv2(1)*p(iz,ix,2)
      do iconv=1,DRVLEN-1
         laplacian = laplacian +       & 
                     deriv2(iconv+1)*( &
                                       p(iz-iconv,ix,2)+p(iz+iconv,ix,2) + &
                                       p(iz,ix-iconv,2)+p(iz,ix+iconv,2)   &
                                     ) 
      end do
      !
      ! atualizacao do campo de pressao:
      p(iz,ix,1) = invpgama * (                                                      & 
                                2.0_dp*p(iz,ix,2)-mgama*p(iz,ix,1)   +               & 
                                vel(iz,ix)*( laplacian + fonte * sourcemask(iz,ix) ) &
                              )
    end do
end do
!
! swap atualiza campo de onda t e t-dt: 
do ix=1,nxx
   do iz=1,nzz
      dswap      = p(iz,ix,2)
      p(iz,ix,2) = p(iz,ix,1)
      p(iz,ix,1) = dswap
   end do
end do
!  
if ( mod(it-1,ndt) == 0.0 ) then
   write(RCID,*) real(p(izr,ixr,2),sp)
end if
!
! armazena instantaneos do campo de onda
if (mod(it-1,nsnap)==0.0) then
    iframe=iframe+1
    if ( mod(iframe,10) == 0 ) write(*,'("writing snapshot # ",i5.5)') iframe
    do ix=ixb,ixe
       itrcmv=itrcmv+1
       vtrace%trcheader%tracl=itrcmv
       vtrace%trcdata(1:nz) = real(p(izb:ize,ix,2),sp)
       ierr=PutSUTrace(SUMVID,itrcmv,vtrace) 
    end do
end if

end do loop_t_free_surf ! fim da evolucao

else
loop_t :do it=1,nt                    ! Loop de evolucao do campo de pressao
!
! Pulso Ricker:
        delt  = (beta*real(it-it0,dp))**2
        fonte = (dx*dx)*exp(-(delt))*(1.0_dp-2.0_dp*delt) ! pulso Ricker
!        
!FD scheme:

!

do ix=DRVLEN,nxx-DRVLEN
   do iz=DRVLEN,nzz-DRVLEN

      gama      = gama_x(ix)+gama_z(iz)
      invpgama  = 1.0_dp/(1.0_dp+gama)
      mgama     = 1.0_dp-gama
      !
      ! calculo do laplaciano:
      laplacian= 2.d0*deriv2(1)*p(iz,ix,2)
      do iconv=1,DRVLEN-1
         laplacian = laplacian +       & 
                     deriv2(iconv+1)*( &
                                       p(iz-iconv,ix,2)+p(iz+iconv,ix,2) + &
                                       p(iz,ix-iconv,2)+p(iz,ix+iconv,2)   &
                                     ) 
      end do
      !
      ! atualizacao do campo de pressao:
      p(iz,ix,1) = invpgama * (                                                      & 
                                2.0_dp*p(iz,ix,2)-mgama*p(iz,ix,1)   +               & 
                                vel(iz,ix)*( laplacian + fonte * sourcemask(iz,ix) ) &
                              )
    end do
end do
!
! swap atualiza campo de onda t e t-dt: 
do ix=1,nxx
   do iz=1,nzz
      dswap      = p(iz,ix,2)
      p(iz,ix,2) = p(iz,ix,1)
      p(iz,ix,1) = dswap
   end do
end do
!  
if ( mod(it-1,ndt) == 0.0 ) then
   write(RCID,*) real(p(izr,ixr,2),sp)
end if
!
! armazena instantaneos do campo de onda
if (mod(it-1,nsnap)==0.0) then
    iframe=iframe+1
    if ( mod(iframe,10) == 0 ) write(*,'("writing snapshot # ",i5.5)') iframe
    do ix=ixb,ixe
       itrcmv=itrcmv+1
       vtrace%trcheader%tracl=itrcmv
       vtrace%trcdata(1:nz) = real(p(izb:ize,ix,2),sp)
       ierr=PutSUTrace(SUMVID,itrcmv,vtrace) 
    end do
end if

end do loop_t ! fim da evolucao
end if
close(RCID)
close(SUMVID)
!
deallocate(p,vel,gama_x,gama_z)
!                     
end program afd_2d
