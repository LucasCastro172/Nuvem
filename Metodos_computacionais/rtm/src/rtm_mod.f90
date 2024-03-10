module rtm_mod
 use mpi
 use number_types,only:fp,sp,dp
 use su_types_mod
 use su_io_mod

implicit none

type,public :: modelagem

!========================================================
!   VARIAVEIS MPI
!========================================================
integer:: adm_myrank
integer:: adm_nproc
integer:: adm_master=0
integer:: adm_source
integer:: adm_dest
integer:: adm_length
integer:: adm_tag = 0
integer:: adm_istatus(MPI_STATUS_SIZE)
integer:: nproc_active

!=======================================================
!  VARIAVEIS CPML
!=======================================================
! Sem born
real(fp),allocatable       :: sigma_z(:),sigma_x(:)
real(fp),allocatable       :: alpha_x(:),alpha_z(:)
real(fp),allocatable 	   :: bx(:),bz(:),ax(:),az(:)
real(fp),allocatable 	   :: psi_x(:,:),psi_z(:,:)
real(fp),allocatable 	   :: eta_x(:,:),eta_z(:,:)
! Com born
real(fp),allocatable       :: sigma_z_b(:),sigma_x_b(:)
real(fp),allocatable       :: alpha_x_b(:),alpha_z_b(:)
real(fp),allocatable 	   :: bx_b(:),bz_b(:),ax_b(:),az_b(:)
real(fp),allocatable 	   :: psi_x_b(:,:),psi_z_b(:,:)
real(fp),allocatable 	   :: eta_x_b(:,:),eta_z_b(:,:)

!=======================================================
!  MODELAGEM ALOCAVEIS
!=======================================================
real(fp),allocatable       :: dvel(:,:)
real(sp),allocatable       :: vel_b(:,:)
real(sp),allocatable       :: velB(:,:)
real(fp),allocatable       :: p(:,:,:)
real(fp),allocatable       :: pb(:,:,:)
real(sp),allocatable       :: imag(:,:)
real(sp),allocatable       :: imagshot(:,:)
real(sp),allocatable       :: imag_filter(:,:)
real(sp),allocatable       :: imag_hess(:,:)
real(sp),allocatable       :: p_hess(:,:)
real(sp),allocatable       :: psrc(:,:)
real(sp),allocatable       :: prcv(:,:)
real(sp),allocatable       :: gradient1d(:)
real(sp),allocatable       :: gradient2d(:,:)
real(fp),allocatable       :: funcao_obj(:)
!=======================================================
!   VARIAVEIS DE ENTRADA
!=======================================================

real(fp)  :: freq,ttotal,dtrec
real(fp)  :: xsmin,xsmax,dxs,zs
real(fp)  :: xg_offmin,xg_offmax,dxg,zg
real(fp)  :: z0,x0,dz,dx
real(fp)  :: funcao_objetivo
integer   :: nz,nx
character(200) :: obs_data

!======================================================
!  VARIAVEIS DA MALHA
!======================================================
integer   :: nxx,nzz
integer   :: izb,ize,ixb,ixe
integer   :: nborda


!=====================================================
!  VARIAVEIS FONTE / ESTABILIDADE / etc...
!=====================================================
real(fp) :: dt, beta
integer  :: ns, ndtrec, ndt,nt, ndt_0,nt_0,dt_0
real(fp) :: pi
integer  :: nframes
integer  :: ndim
real(fp) :: fcost



!====================================================
!  MODELAGEM
!====================================================
real(fp), allocatable :: vel(:,:)
real(fp), allocatable :: model(:)
real(fp)              :: vmax,vmin,vmax_0,vmin_0
type(su_trace)        :: imagtrace, trace
type(su_trace),allocatable        :: shotgather_o(:),shotgather_m(:)
real(fp)              :: scalco, scalel
integer               :: nshots,ntrc
integer               :: isx,isz
integer,allocatable   :: shotidx(:)
integer,allocatable   :: ig1(:),ig2(:)
integer               :: nrec
integer               :: itracl,itracr
integer               :: ilanco





contains
procedure,public :: acquisition_init
procedure,public :: acquisition_seismic2
procedure,public :: evolucao
procedure,public :: evolucao_2
procedure,public :: evolucao_born
procedure,public :: perfis_borda
procedure,public :: alocate_fwi
procedure,public :: input_model

end type


contains	


subroutine acquisition_init(fwi)
class(modelagem)   ::  fwi
integer            ::  ierr
character(200)     ::  parfile,geofile,obsdata,velhdr,task_t
integer,parameter  ::  PARFID = 30
integer,parameter  ::  GEOFID = 31
integer,parameter  ::  OBS    = 32
integer,parameter  ::  SUID   = 33
integer, parameter ::  NAMELEN=200
real(fp),parameter ::  sec2mcs= 1.0e6
integer            ::  ntrcmax,itrc,ishot
real(fp)           ::  dxmax,dxmax_0
real(fp),parameter ::  cfl     = 0.25d0
integer            ::  nsnap
integer,parameter  ::  DRVLEN=8
 real(fp),parameter::  deriv2(DRVLEN)=(/-3.02359410d0,    &
                                        1.75000000d0,     &
                                        -0.291666667d0,   &
                                        0.0648148148d0,   &
                                        -0.0132575758d0,  &
                                        0.00212121212d0,  &
                                        -0.000226625227d0,&
                                        0.0000118928690d0/)
real(fp),parameter :: pi = 3.14159265358979323846264338328_fp
integer :: ix,iz

!* Start MPI *!
   CALL MPI_Init(ierr)
   if (ierr /= 0) then
      stop 'Error in MPI_Init'
   end if

   !* Find process rank *!
   CALL MPI_Comm_rank(MPI_COMM_WORLD, fwi%adm_myrank, ierr);
   if (ierr /= 0) then
      stop 'Error in MPI_Comm_rank'
   end if
    
   !* Find out number of processes *!
   CALL MPI_Comm_size(MPI_COMM_WORLD, fwi%adm_nproc, ierr);
   if (ierr /= 0) then
      stop 'Error in MPI_Comm_size'
   end if    



if ( fwi%adm_myrank == fwi%adm_master ) then
  write(*,'("Enter parameter filename")')
  read(*,'(a)') parfile
  write(*,'("Enter geometry filename")')
  read(*,'(a)') geofile
end if

! Brodcasting input files names
fwi%adm_source = fwi%adm_master
CALL MPI_Bcast(parfile,NAMELEN,MPI_CHARACTER,fwi%adm_source,MPI_COMM_WORLD,ierr)
if ( fwi%adm_myrank == fwi%adm_master ) then
   if(ierr /= 0) then
      stop 'Error broadcasting'
   end if
end if

fwi%adm_source = fwi%adm_master
CALL MPI_Bcast(geofile,NAMELEN,MPI_CHARACTER,fwi%adm_source,MPI_COMM_WORLD,ierr)
if ( fwi%adm_myrank == fwi%adm_master ) then
   if(ierr /= 0) then
      stop 'Error broadcasting'
   end if
end if



open(PARFID,file=parfile,status='old',iostat=ierr, action='read')
read(PARFID,'(A)') velhdr
close(PARFID)



open(GEOFID,file=geofile,status='old',iostat=ierr)
read(GEOFID,*) fwi%freq
read(GEOFID,*) fwi%ttotal
read(GEOFID,*) fwi%dtrec
read(GEOFID,*) fwi%xsmin,fwi%xsmax,fwi%dxs
read(GEOFID,*) fwi%zs
read(GEOFID,*) fwi%xg_offmin,fwi%xg_offmax,fwi%dxg
read(GEOFID,*) fwi%zg  
close(GEOFID)


if ( fwi%adm_myrank == fwi%adm_master ) then
  write(*,'(72("-"))') 
  write(*,'("Vel. model  : ",a)') trim(velhdr)
  write(*,'("AQUISITION GEOMETRY :")')
  write(*,'("Source peak frequency             : ",f10.3," Hz")')fwi%freq
  write(*,'("Sampling interval (in second)     : ",f10.3)')fwi%ttotal
  write(*,'("Simulation time                   : ",f10.3)')fwi%dtrec
  write(*,'("Shot points (xsmin,xsmax,dxs)     : ",f10.3,",",f10.3,",",f10.3)')fwi%xsmin,fwi%xsmax,fwi%dxs
  write(*,'("sources depth                     : ",f10.3)')fwi%zs
  write(*,'("Rec. offsets (xoffmin,xoffmax,dx) : ",f10.3,",",f10.3,",",f10.3)')fwi%xg_offmin,fwi%xg_offmax,fwi%dxg
  write(*,'("receivers depth                   : ",f10.3)')fwi%zg
end if
 
call input_model(fwi,velhdr)


!=================================================================
! verifica condicoes de  estabilidade e dispersao numerica               ANALISADO
!=================================================================
dxmax = 0.25d0*fwi%vmin / (2.5d0*fwi%freq)
if (fwi%dx>dxmax) then
print*,'A frequência está muito alta para o intervalo do grid'
print*,'dxmax=', dxmax
stop
end if
fwi%dt    = cfl *sqrt(2.0/sum(abs(deriv2(:)))) * fwi%dx   / fwi%vmax
!dt    = cfl * dx   / vmax
fwi%ndt   = floor(fwi%dtrec/fwi%dt+0.5) + 1      
!ndt  = ceiling(dtrec/dt)
fwi%ndtrec   = ceiling(fwi%dtrec/fwi%dt)
! intervalo de tempo de registro multiplo inteiro
! do intervalo de modelagem  dt
!
if ( fwi%ndt > 1) then
   fwi%dt = fwi%dtrec/real(fwi%ndt)
else
   fwi%dt = fwi%dtrec
end if
if ( fwi%adm_myrank == fwi%adm_master ) then
print*, 'Tempo total', fwi%ttotal
end if
fwi%nt =ceiling(fwi%ttotal/fwi%dt)          
if( fwi%dx > dxmax ) then
   stop 'freq muito alta'
end if



allocate(fwi%velB(fwi%nzz,fwi%nxx))
allocate(fwi%dvel(fwi%nzz,fwi%nxx))
!==================================================================
!       calculo da refletividade
!==================================================================
!=========== Diferença entre as velocidades ======================
do ix=1,fwi%nxx
	do iz=1,fwi%nzz
!fwi%vel(iz,ix)= 2.d0*(fwi%dt**2)*((fwi%vel(iz,ix)-fwi%vel_b(iz,ix))/fwi%vel_b(iz,ix))*fwi%vel_b(iz,ix)**2 ! Diferença entre as velocidade real e de referencia divido por velocidade de referencia (normalizado)
!fwi%dvel(iz,ix) = (fwi%vel(iz,ix) - fwi%vel(iz,ix))/(fwi%vel(iz+1,ix+1) + fwi%vel(iz,ix))  
fwi%velB(iz,ix) = 2.d0*abs((fwi%vel(iz,ix)-fwi%vel_b(iz,ix))/fwi%vel_b(iz,ix))
	end do
end do




!====================================================================
! NORMALIZAÇÃO DAS VELOCIDADES
!====================================================================


!do ix=1,fwi%nxx
!   do iz=1,fwi%nzz
!      fwi%vel(iz,ix)=(fwi%vel(iz,ix)*fwi%dt_0)**2
!   end do
!end do


! velocidade de propagacao smooth normalizada: 
do ix=1,fwi%nxx
   do iz=1,fwi%nzz
      fwi%vel_b(iz,ix)=(fwi%vel_b(iz,ix)*fwi%dt)**2
   end do
end do




allocate(fwi%imagtrace%trcdata(fwi%nz),stat = ierr) ! traco imagem

! header do traco imagem
  fwi%imagtrace%trcheader%ns=fwi%nz
  fwi%imagtrace%trcheader%dt=int(sec2mcs*fwi%dtrec,2)
  fwi%imagtrace%trcheader%f1=fwi%z0
  fwi%imagtrace%trcheader%f2=fwi%x0
  fwi%imagtrace%trcheader%d1=fwi%dz
  fwi%imagtrace%trcheader%d2=fwi%dx


!=========================================
!Fonte 
 fwi%pi = pi
 fwi%beta= fwi%pi*fwi%freq*fwi%dt
!=========================================
 fwi%ns = ceiling(fwi%ttotal/fwi%dtrec)   
 fwi%nshots=floor((fwi%xsmax-fwi%xsmin)/fwi%dxs)+1
 fwi%nrec=abs(floor((fwi%xg_offmax-fwi%xg_offmin)/fwi%dxg))+1

print*, fwi%nrec
!===================================================================================
! 			ALOCAÇÃO DAS VARIÁVEIS TIPO SU
!===================================================================================

allocate(fwi%shotgather_m(fwi%nrec),stat=ierr)
if ( ierr /= 0 ) then
     stop 'falha de allocacao shotgather(:)'
end if
!
do itrc=1,fwi%nrec
   allocate(fwi%shotgather_m(itrc)%trcdata(fwi%ns),STAT=ierr)
   if ( ierr /= 0 ) then
       stop 'falha de allocacao shotgather(:)%trcdata(:)'
   end if
end do

allocate(fwi%ig1(fwi%nrec),fwi%ig2(fwi%nrec),stat=ierr)
if ( ierr /= 0 ) then
     stop 'falha de allocacao igx ou igz(:,:)'
end if

!===================================================================

! Prints 
if (fwi%adm_myrank == fwi%adm_master) then
  write ( *,'(a)')' '
  write(*,'("Number of shots      : ",I5)')fwi%nshots
  write(*,'("Number of receivers  : ",I5)')fwi%nrec
  write(*,'("Number of samples    : ",I5)')fwi%ns
  write ( *,'(a)') ' '
end if


!=======================================================================
! 		ATRIBUINDO GEOMETRIA DE AQUISIÇÃO
!=======================================================================
 fwi%itracl = 0
 fwi%itracr = 0
 fwi%shotgather_m(1:fwi%nrec)%trcheader%ns=int(fwi%ns,2)
 fwi%shotgather_m(1:fwi%nrec)%trcheader%dt=int(sec2mcs*fwi%dtrec,2)
 fwi%scalco=1.0_fp
 fwi%scalel=1.0_fp
 fwi%shotgather_m(1:fwi%nrec)%trcheader%scalco=int(fwi%scalco,2)
 fwi%shotgather_m(1:fwi%nrec)%trcheader%scalel=int(fwi%scalel,2)
 fwi%shotgather_m(1:fwi%nrec)%trcheader%sdepth=int(fwi%scalel*fwi%zs,4)
 fwi%shotgather_m(1:fwi%nrec)%trcheader%selev=int(-fwi%scalel*fwi%zs,4)
 fwi%shotgather_m(1:fwi%nrec)%trcheader%gelev=int(-fwi%scalel*fwi%zg,4)

!----------------------------------------------
! SUBROTINA PARA ALOCAR CAMPOS 
!----------------------------------------------
 call alocate_fwi(fwi)
!----------------------------------------------
! SUBROTINA PARA ATRIBUIR PERFIS DE ABSRORÇÃO
!----------------------------------------------
 call perfis_borda(fwi)

! fwi%ndim=fwi%nz*fwi%nx

! inervalo entre instantaneos do campo de onda
nsnap  = max(1,fwi%nt/fwi%nframes)
if (fwi%adm_myrank == fwi%adm_master) then
print *, '## of time steps  = ', fwi%nt
print *, 'time interval dt = ',fwi%dt
end if


end subroutine acquisition_init



subroutine acquisition_seismic2(fwi)
class(modelagem) :: fwi
integer          :: ierr
integer,parameter:: SUID = 101
integer,parameter:: SUIDP= 102
integer,parameter:: MVBWID=103
integer,parameter:: ADJFID=104
integer,parameter:: SUIDIMG=105
real(fp),external:: interp_trace
character(200)   :: sudata = " "
character(200)   :: fmovie = " "
character(200)   :: adjoint = " "
character(200)   :: suimag = " "
character(200)   :: ps_hess = " "
character(200)   :: task_t
integer          :: itrcrec,ishot,irec,ircv,itrc,reclen,itrcmv,iframe,itrec,it0,it,iz,ix
real(fp)         :: xs,xg,t
real(fp)         :: delt,fonte,source,fator, izm,izp,ixm,ixp,laplacian
integer          :: iconv
real(sp)         :: imag_hess_max
real(fp),allocatable:: window_x(:),window_z(:)
integer,parameter    :: WDWLEN=5
real(fp),parameter   :: blackmann(WDWLEN)=(/0.893012701892219d0, &
                                             0.630000000000000d0, &
                                             0.340000000000000d0, &
                                             0.130000000000000d0, &
                                             0.026987298107781d0 /)

integer,parameter  :: DRVLEN=8
real(fp),parameter :: deriv2(DRVLEN)=(/-3.02359410d0,    &
                                        1.75000000d0,     &
                                        -0.291666667d0,   &
                                        0.0648148148d0,   &
                                        -0.0132575758d0,  &
                                        0.00212121212d0,  &
                                        -0.000226625227d0,&
                                        0.0000118928690d0/)



fwi%ilanco=0


sudata='dado_observado.su'
  ierr =  OpenSU(SUIDP,trim(sudata),'w',fwi%ns)
  if ( fwi%adm_myrank == fwi%adm_master ) then
    if ( ierr /= 0 ) then
        write(*,'("Error opening SU file: ",A)') trim(sudata)
        stop 
    end if
  end if


 irec= 0
 fwi%itracr   = 0
 fwi%itracl   = 0
 itrcrec      = 0
! task_t= "vetormatriz"
! call transformacao(fwi,task_t)
! call estender_model(fwi)

! velocidade de propagacao normalizada: 

!fwi%velB(:,:)=2*(fwi%vel_b(:,:))*(fwi%vel(:,:)**2)*fwi%dt**2 
!do ix=1,fwi%nxx
!   do iz=1,fwi%nzz
!      fwi%vel(iz,ix)=(fwi%vel(iz,ix)*fwi%dt)**2
!   end do
!end do
! Allocacao de memoria para RTM   
  
  
  allocate(window_x(fwi%nx),stat=ierr)   ! taper na janela de imagem
  if (ierr/=0) then                      
    print*,'Failed to allocate window_x(nx)';stop
  end if

  allocate(window_z(fwi%nz),stat=ierr)   ! taper na janela de imagem
  if (ierr/=0) then                      
     print*,'Failed to allocate window_z(nz)';stop
  end if



  window_z(1:fwi%nz) = 1.0
  do iz=1,WDWLEN
     window_z(1+WDWLEN-iz)  = blackmann(iz)
     window_z(fwi%nz-WDWLEN+iz) = blackmann(iz)
  end do
  window_x(1:fwi%nx) = 1.0
  do ix=1,WDWLEN
     window_x(1+WDWLEN-ix)  = blackmann(ix)
     window_x(fwi%nx-WDWLEN+ix) = blackmann(ix)
  end do



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !  
! ! Inicializacao da imagem RTM  
! 
 fwi%imagshot(1:fwi%nz,1:fwi%nx) = 0.0
  if (fwi%adm_myrank == fwi%adm_master) then
      fwi%imag(1:fwi%nz,1:fwi%nx) = 0.0
      fwi%imag_filter(1:fwi%nz,1:fwi%nx) = 0.0     
  end if
  fwi%prcv(1:fwi%nz,1:fwi%nx)= 0.0
  fwi%psrc(1:fwi%nz,1:fwi%nx)= 0.0


allocate(fwi%psi_x_b(fwi%nzz,fwi%nxx), &
	 fwi%psi_z_b(fwi%nzz,fwi%nxx), &
	 fwi%eta_x_b(fwi%nzz,fwi%nxx), &
	 fwi%eta_z_b(fwi%nzz,fwi%nxx))

!==============================================================================
!                   INICIO LOOP-SHOT 
!==============================================================================
LOOP_SHOTS: do ishot=fwi%adm_myrank+1,fwi%nshots,fwi%adm_nproc

     fwi%itracr = (ishot-1)*fwi%nrec
     fwi%itracl = (ishot-1)*fwi%nrec 
  
     ! Compunting position of the shot 
     xs  = fwi%xsmin + real(ishot-1,fp)*fwi%dxs
     fwi%isz = floor((fwi%zs-fwi%z0)/fwi%dz) + fwi%nborda + 1 
     fwi%isx = floor((xs-fwi%x0)/fwi%dx) + fwi%nborda + 1 
     !print*, isz,isx
  
     do irec=1,fwi%nrec
       !xg  =  xs + fwi%xg_offmin + real(irec-1,fp)*fwi%dxg
        xg  =   fwi%xg_offmin + real(irec-1,fp)*fwi%dxg
        fwi%ig1(irec) = floor((fwi%zg-fwi%z0)/fwi%dz) + fwi%nborda + 1
        fwi%ig2(irec) = floor((xg-fwi%x0)/fwi%dx) + fwi%nborda + 1 
        fwi%shotgather_m(irec)%trcheader%sx=int(fwi%scalco*xs,4)
        fwi%shotgather_m(irec)%trcheader%gx=int(fwi%scalco*xg,4)
        fwi%shotgather_m(irec)%trcheader%offset=int(fwi%scalco*(xg-xs),4)
        fwi%itracl = fwi%itracl + 1
        fwi%itracr = fwi%itracr + 1
        fwi%shotgather_m(irec)%trcheader%tracl=fwi%itracl
        fwi%shotgather_m(irec)%trcheader%tracr=fwi%itracr
        
     end do


reclen = fwi%nz*fwi%nx*SP
fmovie(1:9)='fwi_save_'
write(fmovie(10:14),'(i5.5)') fwi%adm_myrank
fmovie(15:18)='.bin'
open(MVBWID,file=trim(fmovie),  & 
          form='unformatted',   &
          access='direct',      &
          recl=RECLEN,          &
          status='replace',     &
          action='write',       &
          iostat=ierr)
!
if( ierr /= 0) then
    write(*,'("Erro ao abrir arquivo ",(a))') trim(fmovie)
    stop 
end if


write(*,'("Iniciando a propagacao do campo da fonte, tiro ",i5.5," de ",i5.5)') ishot,fwi%nshots

fwi%p(:,:,:)= 0.0_fp
fwi%pb(:,:,:)= 0.0_fp
fwi%psi_x(:,:)= 0.0_fp
fwi%psi_z(:,:)= 0.0_fp
fwi%eta_x(:,:)= 0.0_fp
fwi%eta_z(:,:)= 0.0_fp

fwi%psi_x_b(:,:)= 0.0_fp
fwi%psi_z_b(:,:)= 0.0_fp
fwi%eta_x_b(:,:)= 0.0_fp
fwi%eta_z_b(:,:)= 0.0_fp


!
 itrcmv = 0
 itrec  = 0
 iframe = 0 
 irec   = 0

it0 = floor(sqrt(6.0/fwi%pi)/(fwi%freq*fwi%dt))          ! pico do pulso Ricker 

loop_t :do it=1,fwi%nt
	delt  = (fwi%beta*real(it-it0,fp))**2
        fonte =  exp(-(delt))*(1.0_fp-2.0_fp*delt) ! pulso Ricker
        fwi%p(2,fwi%isz,fwi%isx) =  fwi%p(2,fwi%isz,fwi%isx) + fwi%vel_b(fwi%isz,fwi%isx)*fonte
       
 call evolucao(fwi)
 
 call evolucao_born(fwi)
 !recording traces
 if ( mod(it-1,fwi%ndtrec) == 0.0 ) then
      itrec = itrec + 1
      do irec=1,fwi%nrec
         fwi%shotgather_m(irec)%trcdata(itrec)  = real(fwi%pb(2,fwi%ig1(irec),fwi%ig2(irec)),sp)      
      end do
 end if
  
if (mod(it-1,fwi%ndtrec) == 0.0) then
 iframe=iframe+1
        write(MVBWID,rec=iframe) real(fwi%p(2,fwi%izb:fwi%ize,fwi%ixb:fwi%ixe),sp)
    
end if




end do loop_t
close(MVBWID)
fwi%nframes = iframe


! Writing shot-gather
 itrcrec=(ishot-1)*fwi%nrec 
 do irec=1,fwi%nrec
      itrcrec=itrcrec + 1
      ierr = PutSUTrace(SUIDP,itrcrec,fwi%shotgather_m(irec))
      if ( ierr /= 0 ) then
           write(*,'("Fail to write SU file: ",A)') trim(sudata)
           stop 
       end if
     
 end do







! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!! RETROPROPAGACAO DO CAMPO: INICIO                                   !!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!========================================================================
!              CARREGANDO DADOS
!========================================================================
reclen = fwi%nz*fwi%nx*SP
open(MVBWID,file=trim(fmovie),  & 
          form='unformatted',   &
          access='direct',      &
          recl=RECLEN,          &
          status='replace',     &
          action='read',       &
          iostat=ierr)
!
if( ierr /= 0) then
    write(*,'("Erro ao abrir arquivo ",(a))') trim(fmovie)
    stop 
end if



fwi%p(:,:,:) = 0.d0
iframe = 0 
irec   = 0
 
loop_t_backward :do it = 1,fwi%nt 

  t = real(it-1,fp)*fwi%dt
  
  !---------------------------------------------------- 
  ! source injection:
  do irec=1,fwi%nrec
     
     fonte = interp_trace(fwi%ttotal-t,fwi%ns,0.0,fwi%dtrec,fwi%shotgather_m(irec)%trcdata)
     iz = fwi%ig1(irec)
     ix = fwi%ig2(irec)
 
     
     fwi%p(2,iz,ix)= fwi%p(2,iz,ix) + fwi%vel_b(iz,ix)*fonte
     !p1(iz,ix)= p1(iz,ix) + model%vel(iz,ix)*fonte*(dx*dx)
  enddo
 
 call evolucao(fwi)

!===============================================================================================
! Examinar
!===============================================================================================
 if (mod(it-1,fwi%ndtrec) == 0.0) then

      iframe=iframe+1
                                  
      do ix=1,fwi%nx
         do iz=1,fwi%nz
            !psrc(iz,ix) = real(model%pfield(izb+iz-1,ixb+ix-1,2),sp)
             fwi%psrc(iz,ix) = window_z(iz)*window_x(ix)*real(fwi%p(2,fwi%izb+iz-1,fwi%ixb+ix-1),sp)
             !fwi%psrc(iz,ix) = real(fwi%p(2,fwi%izb+iz-1,fwi%ixb+ix-1),sp)
            !psrc(iz,ix) = real(p2(izb+iz-1,ixb+ix-1),sp)
         end do
      end do
      
      !write(SUIDRTM,rec=iframe) real(psrc(1:nz,1:nx),sp)
      irec = fwi%nframes +1 -iframe   
      read(MVBWID,rec=irec)fwi%prcv(1:fwi%nz,1:fwi%nx)

!------------------------------------------!
!  CONDICAO DE IMAGEM CROSS-CORRELACAO:    !
!------------------------------------------!
     
      do ix=1,fwi%nx
          do iz=1,fwi%nz
          fwi%imagshot(iz,ix) = fwi%imagshot(iz,ix) + fwi%psrc(iz,ix)*fwi%prcv(iz,ix)
        end do
      end do
      
    end if



end do loop_t_backward
close(MVBWID)

end do LOOP_SHOTS
close(SUIDP)
close(SUID)


 CALL MPI_REDUCE(fwi%imagshot,fwi%imag,fwi%nz*fwi%nx,MPI_REAL, &
                 MPI_SUM,fwi%adm_master,MPI_COMM_WORLD,ierr)
 if (ierr /= 0) then
    write(*,'(" MPI_REDUCE failed ")'); stop
 end if  

if ( fwi%adm_myrank == fwi%adm_master ) then


! ! Filtro Laplaciano: remover ruido de baixo numero de onda
do ix=1,fwi%nx
    do iz=1,fwi%nz
         !
        laplacian= deriv2(1)*(fwi%imag(iz,ix)+fwi%imag(iz,ix))
        do iconv=1,DRVLEN-1
            izm=max( 1, iz-iconv)
            izp=min(fwi%nz, iz+iconv)
            laplacian = laplacian +       & 
                        deriv2(iconv+1)*( &
                                          fwi%imag(izm,ix ) + fwi%imag(izp,ix ) &
                                        ) 
         end do
         !
         do iconv=1,DRVLEN-1
            ixm=max( 1, ix-iconv)
            ixp=min(fwi%	nx, ix+iconv)
            laplacian = laplacian +       & 
                        deriv2(iconv+1)*( &
                                          fwi%imag(iz ,ixm) + fwi%imag(iz ,ixp)   &
                                        ) 
         end do
 
         fwi%imag_filter(iz,ix) = real(laplacian)
    end do
 end do
   suimag='rtm_image.su'
   ierr =  OpenSU(SUIDIMG, suimag,'w',fwi%nz)
   IF(ierr /= 0)  PRINT*,"OpenSU falhou para arquivo image!"

   do ix=1,fwi%nx
       fwi%imagtrace%trcheader%tracl = ix
       fwi%imagtrace%trcdata(1:fwi%nz)   = fwi%imag(1:fwi%nz,ix)
       ierr=PutSUTrace(SUIDIMG,ix,fwi%imagtrace) 
   end do
   close(SUIDIMG)  
end if 

end subroutine acquisition_seismic2


!====================================================================================================
!                         EVOLUÇÃO DA ONDA
!====================================================================================================
subroutine evolucao(fwi)
class(modelagem)   :: fwi
integer :: iz,ix
real(fp):: pderiv1_x,pderiv2_x,pderiv1_z,pderiv2_z
integer :: iconv
real(fp):: dswap
real(sp):: deriv1_psi_x,deriv1_psi_z
integer,parameter  :: DRV1LEN = 8
real(fp),parameter :: deriv1(DRV1LEN)=(/0.0d0,            &
                                        0.875d0,           &
                                        -0.291666667d0,    &
                                        0.0972222222d0,    &
                                        -0.0265151515d0,   &
                                        0.00530303030d0,   &
                                        -0.000679875680d0, &
                                        0.0000416250416d0/)

integer,parameter  :: DRVLEN=8
real(fp),parameter :: deriv2(DRVLEN)=(/-3.02359410d0,    &
                                        1.75000000d0,     &
                                        -0.291666667d0,   &
                                        0.0648148148d0,   &
                                        -0.0132575758d0,  &
                                        0.00212121212d0,  &
                                        -0.000226625227d0,&
                                        0.0000118928690d0/)





do ix=DRVLEN,fwi%nxx-DRVLEN
   do iz=DRVLEN,fwi%nzz-DRVLEN
pderiv1_x= deriv1(1)*fwi%p(2,iz,ix)
   do iconv=1,DRV1LEN-1
      pderiv1_x = pderiv1_x +       & 
                  deriv1(iconv+1)*(fwi%p(2,iz,ix+iconv)-fwi%p(2,iz,ix-iconv) ) 
 end do
 pderiv1_x = pderiv1_x/fwi%dx

! Campo auxiliar da primeira equação em x
 fwi%psi_x(iz,ix)= fwi%bx(ix)*fwi%psi_x(iz,ix)+fwi%ax(ix)*pderiv1_x


! Derivada de primeira ordem da cpml (1,iz,ix) Primeiro campo auxiliar
 deriv1_psi_x= deriv1(1)*fwi%psi_x(iz,ix)
do iconv=1, DRV1LEN-1
     deriv1_psi_x= deriv1_psi_x +            &
                         deriv1(iconv+1)*(fwi%psi_x(iz,ix+iconv) - fwi%psi_x(iz,ix-iconv))
end do
deriv1_psi_x=(deriv1_psi_x/fwi%dx)


! Calculo da derivada de segunda ordem no campo de pressão
! Em x
pderiv2_x=deriv2(1)*fwi%p(2,iz,ix)
   do iconv=1,DRVLEN-1
      pderiv2_x = pderiv2_x +       & 
                  deriv2(iconv+1)*(fwi%p(2,iz,ix-iconv)+fwi%p(2,iz,ix+iconv) ) 
end do
pderiv2_x=pderiv2_x/(fwi%dx**2)

! ! Campo auxiliar da segunda equação em x
  fwi%eta_x(iz,ix)= fwi%bx(ix)*fwi%eta_x(iz,ix)+fwi%ax(ix)*(pderiv2_x+deriv1_psi_x)
!----------------------------------------------------------------------------------

! Calculo da derivada de primeira ordem no campo de pressão (presente)
! Em z
pderiv1_z= deriv1(1)*fwi%p(2,iz,ix)
   do iconv=1,DRV1LEN-1
    pderiv1_z = pderiv1_z +       & 
                  deriv1(iconv+1)*( fwi%p(2,iz+iconv,ix)-fwi%p(2,iz-iconv,ix) ) 
   end do
pderiv1_z=pderiv1_z/fwi%dx

! Campo auxiliar da primeira equação em z
 fwi%psi_z(iz,ix)= fwi%bz(iz)*fwi%psi_z(iz,ix)+fwi%az(iz)*pderiv1_z


!----------------------------------------------------------------------------------
! Calculo da derivada de segunda ordem no campo de pressão
! Em z
pderiv2_z=deriv2(1)*fwi%p(2,iz,ix)
   do iconv=1,DRVLEN-1
      pderiv2_z = pderiv2_z +       & 
                  deriv2(iconv+1)*(fwi%p(2,iz-iconv,ix)+fwi%p(2,iz+iconv,ix) ) 
   end do
pderiv2_z=pderiv2_z/(fwi%dx**2)

! Derivada de primeira ordem da cpml (1,iz,ix) Primeiro campo auxiliar
 deriv1_psi_z= deriv1(1)*fwi%psi_z(iz,ix)
do iconv=1,DRV1LEN-1
	deriv1_psi_z= deriv1_psi_z +            &
                         deriv1(iconv+1)*(fwi%psi_z(iz+iconv,ix) - fwi%psi_z(iz-iconv,ix))
end do
deriv1_psi_z=deriv1_psi_z/fwi%dx
! ! Campo auxiliar da segunda equação em z
  fwi%eta_z(iz,ix)= fwi%bz(iz)*fwi%eta_z(iz,ix)+fwi%az(iz)*(pderiv2_z+deriv1_psi_z)

      !
      ! atualizacao do campo de pressao:
      fwi%p(1,iz,ix) = -fwi%p(1,iz,ix) + 2.d0*fwi%p(2,iz,ix) + fwi%vel_b(iz,ix)*(pderiv2_x + &
                    pderiv2_z + deriv1_psi_x + deriv1_psi_z + fwi%eta_x(iz,ix) + &
                                fwi%eta_z(iz,ix))                                        
      !
     fwi%p(3,iz,ix)= (pderiv2_x+pderiv2_z) !+ deriv1_psi_x + deriv1_psi_z + fwi%eta_x(iz,ix) + &
                 !fwi%eta_z(iz,ix))
    end do
end do
!##############################################################################

! swap atualiza campo de onda t e t-dt: 
do ix=1,fwi%nxx
   do iz=1,fwi%nzz
      dswap      = fwi%p(2,iz,ix)
      fwi%p(2,iz,ix) = fwi%p(1,iz,ix)
      fwi%p(1,iz,ix) = dswap
   end do
end do


end subroutine evolucao

!====================================================================================================
!                         EVOLUÇÃO DA ONDA
!====================================================================================================
subroutine evolucao_2(fwi)
class(modelagem)   :: fwi
integer :: iz,ix
real(fp):: pderiv1_x,pderiv2_x,pderiv1_z,pderiv2_z
integer :: iconv
real(fp):: dswap
real(sp):: deriv1_psi_x,deriv1_psi_z
integer,parameter  :: DRV1LEN = 8
real(fp),parameter :: deriv1(DRV1LEN)=(/0.0d0,            &
                                        0.875d0,           &
                                        -0.291666667d0,    &
                                        0.0972222222d0,    &
                                        -0.0265151515d0,   &
                                        0.00530303030d0,   &
                                        -0.000679875680d0, &
                                        0.0000416250416d0/)

integer,parameter  :: DRVLEN=8
real(fp),parameter :: deriv2(DRVLEN)=(/-3.02359410d0,    &
                                        1.75000000d0,     &
                                        -0.291666667d0,   &
                                        0.0648148148d0,   &
                                        -0.0132575758d0,  &
                                        0.00212121212d0,  &
                                        -0.000226625227d0,&
                                        0.0000118928690d0/)
fwi%psi_x(:,:)= 0.0_fp
fwi%psi_z(:,:)= 0.0_fp
fwi%eta_x(:,:)= 0.0_fp
fwi%eta_z(:,:)= 0.0_fp

pderiv2_x =0.0_fp
pderiv2_z =0.0_fp


do ix=DRVLEN,fwi%nxx-DRVLEN
   do iz=DRVLEN,fwi%nzz-DRVLEN
pderiv1_x= deriv1(1)*fwi%p(2,iz,ix)
   do iconv=1,DRV1LEN-1
      pderiv1_x = pderiv1_x +       & 
                  deriv1(iconv+1)*(fwi%p(2,iz,ix+iconv)-fwi%p(2,iz,ix-iconv) ) 
 end do
 pderiv1_x = pderiv1_x/fwi%dx

! Campo auxiliar da primeira equação em x
 fwi%psi_x(iz,ix)= fwi%bx(ix)*fwi%psi_x(iz,ix)+fwi%ax(ix)*pderiv1_x


! Derivada de primeira ordem da cpml (1,iz,ix) Primeiro campo auxiliar
 deriv1_psi_x= deriv1(1)*fwi%psi_x(iz,ix)
do iconv=1, DRV1LEN-1
     deriv1_psi_x= deriv1_psi_x +            &
                         deriv1(iconv+1)*(fwi%psi_x(iz,ix+iconv) - fwi%psi_x(iz,ix-iconv))
end do
deriv1_psi_x=(deriv1_psi_x/fwi%dx)


! Calculo da derivada de segunda ordem no campo de pressão
! Em x
pderiv2_x=deriv2(1)*fwi%p(2,iz,ix)
   do iconv=1,DRVLEN-1
      pderiv2_x = pderiv2_x +       & 
                  deriv2(iconv+1)*(fwi%p(2,iz,ix-iconv)+fwi%p(2,iz,ix+iconv) ) 
end do
pderiv2_x=pderiv2_x/(fwi%dx**2)

! ! Campo auxiliar da segunda equação em x
  fwi%eta_x(iz,ix)= fwi%bx(ix)*fwi%eta_x(iz,ix)+fwi%ax(ix)*(pderiv2_x+deriv1_psi_x)
!----------------------------------------------------------------------------------

! Calculo da derivada de primeira ordem no campo de pressão (presente)
! Em z
pderiv1_z= deriv1(1)*fwi%p(2,iz,ix)
   do iconv=1,DRV1LEN-1
    pderiv1_z = pderiv1_z +       & 
                  deriv1(iconv+1)*( fwi%p(2,iz+iconv,ix)-fwi%p(2,iz-iconv,ix) ) 
   end do
pderiv1_z=pderiv1_z/fwi%dx

! Campo auxiliar da primeira equação em z
 fwi%psi_z(iz,ix)= fwi%bz(iz)*fwi%psi_z(iz,ix)+fwi%az(iz)*pderiv1_z


!----------------------------------------------------------------------------------
! Calculo da derivada de segunda ordem no campo de pressão
! Em z
pderiv2_z=deriv2(1)*fwi%p(2,iz,ix)
   do iconv=1,DRVLEN-1
      pderiv2_z = pderiv2_z +       & 
                  deriv2(iconv+1)*(fwi%p(2,iz-iconv,ix)+fwi%p(2,iz+iconv,ix) ) 
   end do
pderiv2_z=pderiv2_z/(fwi%dx**2)

! Derivada de primeira ordem da cpml (1,iz,ix) Primeiro campo auxiliar
 deriv1_psi_z= deriv1(1)*fwi%psi_z(iz,ix)
do iconv=1,DRV1LEN-1
	deriv1_psi_z= deriv1_psi_z +            &
                         deriv1(iconv+1)*(fwi%psi_z(iz+iconv,ix) - fwi%psi_z(iz-iconv,ix))
end do
deriv1_psi_z=deriv1_psi_z/fwi%dx
! ! Campo auxiliar da segunda equação em z
  fwi%eta_z(iz,ix)= fwi%bz(iz)*fwi%eta_z(iz,ix)+fwi%az(iz)*(pderiv2_z+deriv1_psi_z)

      !
      ! atualizacao do campo de pressao:
      fwi%p(1,iz,ix) = -fwi%p(1,iz,ix) + 2.d0*fwi%p(2,iz,ix) + fwi%vel_b(iz,ix)*(pderiv2_x + &
                    pderiv2_z + deriv1_psi_x + deriv1_psi_z + fwi%eta_x(iz,ix) + &
                                fwi%eta_z(iz,ix))                                        
      !
     !fwi%p(3,iz,ix)= (pderiv2_x+pderiv2_z) !+ deriv1_psi_x + deriv1_psi_z + fwi%eta_x(iz,ix) + &
                 !fwi%eta_z(iz,ix))
    end do
end do
!##############################################################################

! swap atualiza campo de onda t e t-dt: 
do ix=1,fwi%nxx
   do iz=1,fwi%nzz
      dswap      = fwi%p(2,iz,ix)
      fwi%p(2,iz,ix) = fwi%p(1,iz,ix)
      fwi%p(1,iz,ix) = dswap
   end do
end do


end subroutine evolucao_2
!====================================================================================================
!                         EVOLUÇÃO DA ONDA
!====================================================================================================
subroutine evolucao_born(fwi)
class(modelagem)   :: fwi
integer :: iz,ix
real(fp):: pderiv1_x,pderiv2_x,pderiv1_z,pderiv2_z
integer :: iconv
real(fp):: dswap
real(sp):: deriv1_psi_x,deriv1_psi_z
integer,parameter  :: DRV1LEN = 8
real(fp),parameter :: deriv1(DRV1LEN)=(/0.0d0,            &
                                        0.875d0,           &
                                        -0.291666667d0,    &
                                        0.0972222222d0,    &
                                        -0.0265151515d0,   &
                                        0.00530303030d0,   &
                                        -0.000679875680d0, &
                                        0.0000416250416d0/)

integer,parameter  :: DRVLEN=8
real(fp),parameter :: deriv2(DRVLEN)=(/-3.02359410d0,    &
                                        1.75000000d0,     &
                                        -0.291666667d0,   &
                                        0.0648148148d0,   &
                                        -0.0132575758d0,  &
                                        0.00212121212d0,  &
                                        -0.000226625227d0,&
                                        0.0000118928690d0/)


pderiv2_x =0.0_fp
pderiv2_z =0.0_fp


do ix=DRVLEN,fwi%nxx-DRVLEN
   do iz=DRVLEN,fwi%nzz-DRVLEN
pderiv1_x= deriv1(1)*fwi%pb(2,iz,ix)
   do iconv=1,DRV1LEN-1
      pderiv1_x = pderiv1_x +       & 
                  deriv1(iconv+1)*(fwi%pb(2,iz,ix+iconv)-fwi%pb(2,iz,ix-iconv) ) 
 end do
 pderiv1_x = pderiv1_x/fwi%dx

! Campo auxiliar da primeira equação em x
 fwi%psi_x_b(iz,ix)= fwi%bx(ix)*fwi%psi_x_b(iz,ix)+fwi%ax(ix)*pderiv1_x


! Derivada de primeira ordem da cpml (1,iz,ix) Primeiro campo auxiliar
 deriv1_psi_x= deriv1(1)*fwi%psi_x_b(iz,ix)
do iconv=1, DRV1LEN-1
     deriv1_psi_x= deriv1_psi_x +            &
                         deriv1(iconv+1)*(fwi%psi_x_b(iz,ix+iconv) - fwi%psi_x_b(iz,ix-iconv))
end do
deriv1_psi_x=(deriv1_psi_x/fwi%dx)


! Calculo da derivada de segunda ordem no campo de pressão
! Em x
pderiv2_x=deriv2(1)*fwi%pb(2,iz,ix)
   do iconv=1,DRVLEN-1
      pderiv2_x = pderiv2_x +       & 
                  deriv2(iconv+1)*(fwi%pb(2,iz,ix-iconv)+fwi%pb(2,iz,ix+iconv) ) 
end do
pderiv2_x=pderiv2_x/(fwi%dx**2)

! ! Campo auxiliar da segunda equação em x
  fwi%eta_x_b(iz,ix)= fwi%bx(ix)*fwi%eta_x_b(iz,ix)+fwi%ax(ix)*(pderiv2_x+deriv1_psi_x)
!----------------------------------------------------------------------------------

! Calculo da derivada de primeira ordem no campo de pressão (presente)
! Em z
pderiv1_z= deriv1(1)*fwi%pb(2,iz,ix)
   do iconv=1,DRV1LEN-1
    pderiv1_z = pderiv1_z +       & 
                  deriv1(iconv+1)*( fwi%pb(2,iz+iconv,ix)-fwi%pb(2,iz-iconv,ix) ) 
   end do
pderiv1_z=pderiv1_z/fwi%dx

! Campo auxiliar da primeira equação em z
 fwi%psi_z_b(iz,ix)= fwi%bz(iz)*fwi%psi_z_b(iz,ix)+fwi%az(iz)*pderiv1_z


!----------------------------------------------------------------------------------
! Calculo da derivada de segunda ordem no campo de pressão
! Em z
pderiv2_z=deriv2(1)*fwi%pb(2,iz,ix)
   do iconv=1,DRVLEN-1
      pderiv2_z = pderiv2_z +       & 
                  deriv2(iconv+1)*(fwi%pb(2,iz-iconv,ix)+fwi%pb(2,iz+iconv,ix) ) 
   end do
pderiv2_z=pderiv2_z/(fwi%dx**2)

! Derivada de primeira ordem da cpml (1,iz,ix) Primeiro campo auxiliar
 deriv1_psi_z= deriv1(1)*fwi%psi_z_b(iz,ix)
do iconv=1,DRV1LEN-1
	deriv1_psi_z= deriv1_psi_z +            &
                         deriv1(iconv+1)*(fwi%psi_z_b(iz+iconv,ix) - fwi%psi_z_b(iz-iconv,ix))
end do
deriv1_psi_z=deriv1_psi_z/fwi%dx
! ! Campo auxiliar da segunda equação em z
  fwi%eta_z_b(iz,ix)= fwi%bz(iz)*fwi%eta_z_b(iz,ix)+fwi%az(iz)*(pderiv2_z+deriv1_psi_z)

      !
      ! atualizacao do campo de pressao:
      fwi%pb(1,iz,ix) = -fwi%pb(1,iz,ix) + 2.d0*fwi%pb(2,iz,ix) + fwi%vel_b(iz,ix)*(pderiv2_x + &
                    pderiv2_z + deriv1_psi_x + deriv1_psi_z + fwi%eta_x_b(iz,ix) + &
                                fwi%eta_z_b(iz,ix)) + &
                    fwi%velB(iz,ix)*fwi%p(3,iz,ix)*fwi%vel_b(iz,ix)                                            
      !
     !fwi%p(3,iz,ix)= (pderiv2_x+pderiv2_z) !+ deriv1_psi_x + deriv1_psi_z + fwi%eta_x(iz,ix) + &
                 !fwi%eta_z(iz,ix))
    end do
end do
!##############################################################################

! swap atualiza campo de onda t e t-dt: 
do ix=1,fwi%nxx
   do iz=1,fwi%nzz
      dswap      = fwi%pb(2,iz,ix)
      fwi%pb(2,iz,ix) = fwi%pb(1,iz,ix)
      fwi%pb(1,iz,ix) = dswap
   end do
end do



end subroutine evolucao_born









subroutine perfis_borda(fwi)
implicit none
class(modelagem)  :: fwi
real(fp)        :: RC,RC_cpml
real(sp)        :: sigma,sigma_b, sigma_0,alpha0,sigma_0_b,alpha0_b
integer         :: ix,iz
!=========================================================================================
!Bordas                              	ANALISADA
!=========================================================================================
!                                CPML
! mascara de absorcao: PREENCHIMENTO DO VALOR DE GAMA NAS BORDAS PARA A ABSORÇÃO
!                      DO CAMPO DE ONDA
!==========================================================================================
!Bordas

allocate(fwi%sigma_x(fwi%nxx),        &
         fwi%sigma_z(fwi%nzz),        &
         fwi%alpha_x(fwi%nxx),        &
         fwi%alpha_z(fwi%nzz)) 

allocate(fwi%sigma_x_b(fwi%nxx),        &
         fwi%sigma_z_b(fwi%nzz),        &
         fwi%alpha_x_b(fwi%nxx),        &
         fwi%alpha_z_b(fwi%nzz)) 

! phi
fwi%sigma_x(:) = 0.0_fp
fwi%sigma_z(:) = 0.0_fp


! phi
!fwi%sigma_x_b(:) = 0.0_fp
!fwi%sigma_z_b(:) = 0.0_fp

! Alpha  
fwi%alpha_x(:)= fwi%pi*fwi%freq ! alpha_x(xtotal)
fwi%alpha_z(:)= fwi%pi*fwi%freq ! alpha_z(ztotal) 
!fwi%alpha_x_b(:)= fwi%pi*fwi%freq ! alpha_x(xtotal)
!fwi%alpha_z_b(:)= fwi%pi*fwi%freq ! alpha_z(ztotal)
RC=10d-7
RC_cpml= log(RC)
!sigma_0= -(3.0d0)*vmax*RC_cpml/(2*real(Nborda,dp))
sigma_0= -(1.5)*fwi%vmax*RC_cpml/(fwi%dx*real(fwi%nborda,fp))
!sigma_0_b= -(1.5)*fwi%vmax_0*RC_cpml/(fwi%dx*real(fwi%nborda,fp))

do ix=1,fwi%nborda
      sigma = sigma_0*((real((fwi%nborda-ix+1),fp)/fwi%nborda))**2
      !sigma_b = sigma_0_b*((real((fwi%nborda-ix+1),fp)/fwi%nborda))**2
      fwi%sigma_x(ix) = sigma
      fwi%sigma_x(fwi%nxx-ix+1) = sigma
    
      !fwi%sigma_x_b(ix) = sigma_b
      !fwi%sigma_x_b(fwi%nxx-ix+1) = sigma_b
      !print*,'Sigma=', sigma
end do
do iz=1,fwi%nborda
      sigma = sigma_0*((real(fwi%nborda-iz+1,fp)/fwi%nborda))**2
      fwi%sigma_z(iz) = sigma
      fwi%sigma_z(fwi%nzz-iz+1) = sigma
end do

do ix=1,fwi%nborda
   alpha0= fwi%pi*fwi%freq*real((1.d0-(real(fwi%nborda-ix+1,fp)/fwi%nborda)))
   !alpha0_b= fwi%pi*fwi%freq*real((1.d0-(real(fwi%nborda-ix+1,fp)/fwi%nborda)))
   fwi%alpha_x(ix)= alpha0
   fwi%alpha_x(fwi%nxx-ix+1)= alpha0
   !fwi%alpha_x_b(ix)= alpha0_b
   !fwi%alpha_x_b(fwi%nxx-ix+1)= alpha0_b
   !print*,'alpha=', alpha0
end do 

do iz=1,fwi%nborda
 alpha0= fwi%pi*fwi%freq*real((1.d0-(real(fwi%nborda-iz+1,fp)/fwi%nborda)))
 !alpha0_b= fwi%pi*fwi%freq*real((1.d0-(real(fwi%nborda-iz+1,fp)/fwi%nborda)))
 fwi%alpha_z(iz)= alpha0
 fwi%alpha_z(fwi%nzz-iz+1)= alpha0
 !fwi%alpha_z_b(iz)= alpha0_b
 !fwi%alpha_z_b(fwi%nzz-iz+1)= alpha0_b
end do

fwi%bx(:)=0.0
fwi%bz(:)=0.0
fwi%ax(:)=0.0
fwi%az(:)=0.0


!fwi%bx_b(:)=0.0
!fwi%bz_b(:)=0.0
!fwi%ax_b(:)=0.0
!fwi%az_b(:)=0.0

do ix=1,fwi%nxx
fwi%bx(ix)= exp(-(fwi%alpha_x(ix)+fwi%sigma_x(ix))*fwi%dt)
fwi%ax(ix)= (fwi%sigma_x(ix)/(fwi%alpha_x(ix)+fwi%sigma_x(ix)) ) * (fwi%bx(ix)-1.0d0)
end do

!do ix=1,fwi%nxx
!fwi%bx_b(ix)= exp(-(fwi%alpha_x_b(ix)+fwi%sigma_x_b(ix))*fwi%dt_0)
!fwi%ax_b(ix)= (fwi%sigma_x_b(ix)/(fwi%alpha_x_b(ix)+fwi%sigma_x_b(ix)) ) * (fwi%bx_b(ix)-1.0d0)
!end do
   
do iz=1,fwi%nzz
fwi%bz(iz)= exp(-(fwi%alpha_z(iz)+fwi%sigma_z(iz))*fwi%dt)
fwi%az(iz)= (fwi%sigma_z(iz)/(fwi%alpha_z(iz)+fwi%sigma_z(iz)) )* (fwi%bz(iz)-1.0d0)
end do

!do iz=1,fwi%nzz
!fwi%bz_b(iz)= exp(-(fwi%alpha_z_b(iz)+fwi%sigma_z_b(iz))*fwi%dt_0)
!fwi%az_b(iz)= (fwi%sigma_z_b(iz)/(fwi%alpha_z_b(iz)+fwi%sigma_z_b(iz)))*(fwi%bz_b(iz)-1.0d0)
!end do

end subroutine perfis_borda



!=========================================================================================================
! 				ALOCACAO DE MEMORIA 
!=========================================================================================================
subroutine alocate_fwi(fwi)
implicit none
class(modelagem)   :: fwi
integer          :: ierr

allocate (  fwi%p(3,fwi%nzz,fwi%nxx),           &
	    fwi%pb(3,fwi%nzz,fwi%nxx),           &
            fwi%imag(fwi%nz,fwi%nx),           &  
	    fwi%imagshot(fwi%nz,fwi%nx),&   ! alocacao de memoria na malha incluindo bordas de absorcao
            fwi%imag_filter(fwi%nz,fwi%nx),           &  
            fwi%p_hess(fwi%nz,fwi%nx),&   ! alocacao de memoria na malha incluindo bordas de absorcao
            fwi%imag_hess(fwi%nz,fwi%nx),           & 
            fwi%prcv(fwi%nz,fwi%nx),           &  
            fwi%psrc(fwi%nz,fwi%nx),           &  
            fwi%bx(fwi%nxx),fwi%bz(fwi%nzz),fwi%ax(fwi%nxx), &
            fwi%az(fwi%nzz),                 &
            fwi%psi_x(fwi%nzz,fwi%nxx),fwi%psi_z(fwi%nzz,fwi%nxx), &
            fwi%eta_x(fwi%nzz,fwi%nxx),fwi%eta_z(fwi%nzz,fwi%nxx), &!
            fwi%gradient2d(fwi%nz,fwi%nx),fwi%gradient1d(fwi%nz*fwi%nx), &
            stat=ierr)
if (ierr/=0) then                      ! informa falha na alocacao de memoria
    stop 'falha de alocaçao da malha'
end if

end subroutine alocate_fwi

!============================================================================================
subroutine input_model(fwi,velhdr)
implicit none
class(modelagem)   :: fwi
integer,parameter  :: nborda = 30
integer, parameter :: VELFID = 10
integer, parameter :: PARFI  = 11
integer            :: ierr,irec
integer            :: ix,iz
real(sp)           :: swap
character(200),intent(in)    :: velhdr
character(200)     :: veldat,veldat_smooth


!========================================================================
! 	ABRINDO ARQUIVO DO MODELO DE VELOCIDADE
!========================================================================
open(VELFID, file=trim(velhdr),iostat=ierr) 
if ( ierr /= 0 ) then
     write(*,'("Erro ao abrir arquivo: ",A)') trim(velhdr)
     stop 
end if

read(VELFID,*) fwi%nz,fwi%nx,fwi%z0,fwi%x0,fwi%dz,fwi%dx
read(VELFID,'(A)') veldat
read(VELFID,'(A)') veldat_smooth
close(VELFID)

fwi%nborda = nborda

   fwi%nxx=fwi%nx+2*fwi%nborda
   fwi%ixb=fwi%nborda+1
   fwi%ixe=fwi%nborda+fwi%nx
   fwi%nzz=fwi%nz+2*fwi%nborda
   fwi%izb=fwi%nborda+1
   fwi%ize=fwi%nborda+fwi%nz

allocate(fwi%vel(fwi%nzz,fwi%nxx))


open(PARFI,file=trim(veldat), &
           form='unformatted',&
           access='direct',   &
           recl= sp, iostat=ierr)
  if ( ierr /= 0 ) then
       write(*,'("Failed to open file: ",A)') trim(veldat)
       stop 
  end if

irec = 0
  ! leitura do modelo de velocidade
  do ix=fwi%ixb,fwi%ixe
     do iz=fwi%izb,fwi%ize
        irec = irec + 1
        read(PARFI,rec=irec) swap
        fwi%vel(iz,ix) = real(swap,fp)
     end do
  end do

fwi%vmax=fwi%vel(fwi%izb,fwi%ixb)                      !  
fwi%vmin=fwi%vel(fwi%izb,fwi%ixb)                      !
do ix=fwi%ixb,fwi%ixe                          ! 
   do iz=fwi%izb,fwi%ize                       ! 
      fwi%vmin=min(fwi%vmin,fwi%vel(iz,ix))        !   
      fwi%vmax=max(fwi%vmax,fwi%vel(iz,ix))        !  
   end do
end do

!fwi%vmax= 5848.d0
!fwi%vmin= 1423.d0
! preenche o modelo de velocidades nas bordas 
!
do ix=1,fwi%nborda
   do iz=fwi%izb,fwi%ize
      fwi%vel(iz,ix)          = fwi%vel(iz,fwi%ixb)
      fwi%vel(iz,fwi%ixe+ix) = fwi%vel(iz,fwi%ixe)
   end do
end do

do ix=1,fwi%nxx
   do iz=1,fwi%izb-1
      fwi%vel(iz,ix)=fwi%vel(fwi%izb,ix)
   end do
   do iz=fwi%ize+1,fwi%nzz
      fwi%vel(iz,ix)=fwi%vel(fwi%ize,ix)
   end do
end do
CLOSE(PARFI)
!=============================================================
!          VELOCITY MODEL (SMOOTH) 
!=============================================================
open(PARFI,file=trim(veldat_smooth), &
           form='unformatted',&
           access='direct',   &
           recl= sp, iostat=ierr)
  if ( ierr /= 0 ) then
       write(*,'("Failed to open file: ",A)') trim(veldat_smooth)
       stop 
  end if
allocate(fwi%vel_b(fwi%nzz,fwi%nxx))
irec = 0
  ! leitura do modelo de velocidade
  do ix=fwi%ixb,fwi%ixe
     do iz=fwi%izb,fwi%ize
        irec = irec + 1
        read(PARFI,rec=irec) swap
        fwi%vel_b(iz,ix) = real(swap,fp)
     end do
  end do
! preenche o modelo de velocidades nas bordas 

do ix=1,fwi%nborda
   do iz=fwi%izb,fwi%ize
      fwi%vel_b(iz,ix)          = fwi%vel_b(iz,fwi%ixb)
      fwi%vel_b(iz,fwi%ixe+ix) = fwi%vel_b(iz,fwi%ixe)
   end do
end do

do ix=1,fwi%nxx
   do iz=1,fwi%izb-1
      fwi%vel_b(iz,ix)=fwi%vel_b(fwi%izb,ix)
   end do
   do iz=fwi%ize+1,fwi%nzz
      fwi%vel_b(iz,ix)=fwi%vel_b(fwi%ize,ix)
   end do
end do

fwi%vmax_0=fwi%vel_b(fwi%izb,fwi%ixb)                      !  
fwi%vmin_0=fwi%vel_b(fwi%izb,fwi%ixb)                      !
do ix=fwi%ixb,fwi%ixe                          ! 
   do iz=fwi%izb,fwi%ize                       ! 
      fwi%vmin_0=min(fwi%vmin_0,fwi%vel_b(iz,ix))        !   
      fwi%vmax_0=max(fwi%vmax_0,fwi%vel_b(iz,ix))        !  
   end do
end do

CLOSE(PARFI)



end subroutine input_model
end module rtm_mod

