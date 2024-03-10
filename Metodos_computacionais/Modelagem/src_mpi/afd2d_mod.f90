module afd2d_mod 
!#################################################################### 
! 
! Prototipo para modelagem acustica 2D
! 
!
!####################################################################
!
use mpi
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
integer,parameter::fp=dp
!
integer,parameter :: NAMELEN=512
integer,parameter :: MDLFID=12
integer,parameter :: RCID=15
integer,parameter :: SUDATFID=20
integer,parameter :: DATFID=11
integer,parameter :: SUMVID=16
integer,parameter :: MVFID=19
!
!Parametros
!
real(dp),parameter :: pi      = 3.14159265358979323846264338328_dp
real(dp),parameter :: twopi   = 2.0_dp*pi
real(dp),parameter :: cfl     = 0.25d0
real,parameter     :: sec2mcs = 1.0e6
real,parameter     :: sec2ms  = 1.0e3
real(dp),parameter :: dtrec   = 0.004d0  ! ms
real(dp),parameter :: dtsave  = 0.008d0  ! ms

integer,parameter  :: DRV1LEN = 1
real(dp),parameter :: deriv1(DRV1LEN)=(/1.0d0/)
!
! operadores de diferencas finitas:
integer,parameter  :: DRVLEN=7
real(dp),parameter :: deriv2(DRVLEN)=(/-3.171097979234166, &
                                        1.883703563645937, &
                                       -0.392013173326410, &
                                        0.126419594523541, &
                                       -0.043436978346425, &
                                        0.013633723622292, &
                                       -0.002757740501852 /)
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

! blackmann window 
integer,parameter    :: WDWLEN=5
real(dp),parameter   :: blackmann(WDWLEN)=(/0.893012701892219d0, &
                                            0.630000000000000d0, &
                                            0.340000000000000d0, &
                                            0.130000000000000d0, &
                                            0.026987298107781d0/)

! point source injection window length
integer, parameter :: WDWSRC=6
!
!
! Parametros: 
! =================
! nframes: numero de frames gravados em disco
! nborda : numero de pontos utilizados nas bordas de absorcao
!          Fronteira absorvente aplicada nas bordas da malha
!
integer,parameter  :: nframes = 100
integer,parameter  :: nborda  = 100    !  comp. de onda maximo
real(dp),parameter :: R0pml   = 1.0d-3    ! R_0
!
!integer, parameter :: nfields=14
integer, parameter :: nfields=10
!
! For mpi server-client parallel computations:
!
!# MPI data structure
!  integer:: mpimyrank
!  integer:: mpinproc
!  integer:: mpimaster
!  integer:: mpisource
!  integer:: mpidest
!  integer:: mpilength
!  integer:: mpitag 
!  integer:: mpistat(MPI_STATUS_SIZE)
!  integer:: mpiactive

type, public:: afd2d

!# MPI data structure
  integer:: mpimyrank
  integer:: mpinproc
  integer:: mpimaster
  integer:: mpisource
  integer:: mpidest
  integer:: mpilength
  integer:: mpitag 
  integer:: mpistat(MPI_STATUS_SIZE)
  integer:: mpiactive
!
integer  :: n1,n2,nt
real(dp) :: o1,o2,d1,d2,dt,freq
real(dp) :: vmin,vmax
!
real(dp), allocatable:: gamma1(:), gamma2(:), gradwind1(:),  gradwind2(:)! fronteiras de absorcao
real(dp), allocatable:: wsrc1(:), wsrc2(:)   ! fronteiras de absorcao
real(dp), allocatable:: wrcv1(:,:), wrcv2(:,:)   ! fronteiras de absorcao
real(dp), allocatable:: source_pulse(:) ! fronteiras de absorcao
!
real(dp),allocatable:: vel(:,:)           ! modelo de velocidade
real(dp),pointer:: pfield1(:,:)       ! campo de onda
real(dp),pointer:: pfield2(:,:)       ! campo de onda
real(dp),pointer:: p1(:,:)=>null()        ! campo de onda
real(dp),pointer:: p2(:,:)=>null()        ! campo de onda
real(dp),pointer:: fswap(:,:)=>null()     ! auxiliar pointer
!
type(su_trace) :: vtrace     ! traco instantaneo do campo de onda
type(su_trace) :: trace
real(dp) :: dswap
!
! acquisition
character(:),allocatable :: datafile
character(NAMELEN):: outfile=""
real(dp) :: datanorm
integer  :: mode
integer  :: nshots,ntrec,ngmax,ng
integer  :: ishot_beg,ishot_end,step_shots
integer, allocatable  :: shotid(:)
real(dp)              :: s1,s2
real(dp), allocatable :: g1(:),g2(:)
type(su_trace),allocatable :: shotgather(:) 
type(su_trace),allocatable :: survey(:) 
real(dp) :: ttotal,t0rec,dtrec
integer  :: ndtrec,ndtsave
integer  :: i1min,i1max,i2min,i2max
integer  :: i1beg,i1end,i2beg,i2end
!
character(NAMELEN) :: bornfile
!
real(dp), allocatable :: wdata(:)
integer  :: MVFID
character(NAMELEN) :: moviefile
!
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
character(NAMELEN):: velhdr,velbin
!
integer :: workflow=0 
!
!
contains
!
procedure :: input_model
procedure :: input_data
procedure :: set_fd_solver
procedure :: afd_solver
!
end type
!
interface interp_trace
        module procedure interp_trace_sp
        module procedure interp_trace_dp
end interface
!
contains
!
   subroutine mpi_begin(this,stat)
class(afd2d),intent(inout) :: this
integer, intent(out) :: stat
    !* Start MPI *!
    CALL MPI_Init(stat)
    if ( stat /= 0 ) return

    !* Find process rank *!
    CALL MPI_Comm_rank(MPI_COMM_WORLD, this%mpimyrank, stat);
    if ( stat /= 0 ) return

    !* Find out number of processes *!
    CALL MPI_Comm_size(MPI_COMM_WORLD, this%mpinproc, stat);
    !* create goups and topology *!
    this%mpimaster = 0
!
end subroutine


subroutine mpi_end(this,stat)
class(afd2d),intent(inout) :: this
integer, intent(out) :: stat
!     
    CALL MPI_Finalize(stat)
!
end subroutine mpi_end

!
   function sinc(x)
   !
   ! definicao da funcao sinc(x) = sin(pi x ) / ( pi x )
   !
   implicit none
   integer,parameter :: fp=kind(1.d0)
   real(fp),intent(in) :: x
   real(fp) :: sinc
   !
   real(fp),parameter::pi=3.14159265358979323846264338328d0

   if ( abs(x) > 0.d0 ) then
      sinc = sin(pi*x) / (pi*x)
   else
      sinc = 1.d0
   end if
   end function sinc


   function lanczos(x,n)
! interpolation of source time history 
! Lanczos interpolation kernel
! L(x,n) = sinc(x)*sinc(x/n)
   !
   implicit none
   integer,parameter :: fp=kind(1.d0)
   real(fp),intent(in) :: x
   integer, intent(in) :: n
   real(fp) :: lanczos
   real(fp) :: a
   !
   a = real(n,fp)
   if ( abs(x) < a ) then
      lanczos = sinc(x)*sinc(x/a)
   else
      lanczos = 0.d0
   end if
   end function lanczos

   function interp_trace_sp(t,nt,t0,dt,trace)
   !
   ! Lanczos interpolation: windowed sinc 
   !
   implicit none
   integer,parameter   :: fp=kind(1.d0)
   integer,intent(in)  :: nt        ! nt = numero de amostras no arranjo trace(:);
   real(fp),intent(in) :: t,t0,dt   ! t0 = amostra inicial;
                                    ! dt = intervalo de amostragem original;
                                    ! t  = tempo de amostragem  
   real(sp),intent(in) :: trace(nt) ! arranjo de dados 
   real(fp) :: interp_trace_sp         ! valor interpolado trace(t)
   !
   real(fp),parameter :: tol = 1.0d-3
   integer,parameter  :: nlag = 3
   integer :: it,itb,ite
   real(fp) :: aux0,aux1

   aux0 = (t - t0) / dt
   it  = floor(aux0)+1
   !
   if ( abs(aux0-real(it-1,fp)) < tol ) then
     interp_trace_sp = real(trace(it),fp)
   end if
   !
   interp_trace_sp = 0.d0
   if ( t0 <= t .and. t <= t0+real(nt-1,fp)*dt ) then
      itb = max(it-nlag,1)
      ite = min(it+nlag,nt)
      do it=itb,ite
         aux1 = aux0 - real(it-1,fp)
         interp_trace_sp = interp_trace_sp + real(trace(it),fp)*lanczos(aux1,nlag)
      end do
   end if
   end function interp_trace_sp

   function interp_trace_dp(t,nt,t0,dt,trace)
   !
   ! Lanczos interpolation: windowed sinc 
   !
   implicit none
   integer,parameter   :: fp=kind(1.d0)
   integer,intent(in)  :: nt        ! nt = numero de amostras no arranjo trace(:);
   real(fp),intent(in) :: t,t0,dt   ! t0 = amostra inicial;
                                    ! dt = intervalo de amostragem original;
                                    ! t  = tempo de amostragem  
   real(dp),intent(in) :: trace(nt) ! arranjo de dados 
   real(fp) :: interp_trace_dp         ! valor interpolado trace(t)
   !
   real(fp),parameter :: tol = 1.0d-3
   integer,parameter  :: nlag = 3
   integer :: it,itb,ite
   real(fp) :: aux0,aux1

   aux0 = (t - t0) / dt
   it  = floor(aux0)+1
   !
   if ( abs(aux0-real(it-1,fp)) < tol ) then
     interp_trace_dp = real(trace(it),fp)
   end if
   !
   interp_trace_dp = 0.d0
   if ( t0 <= t .and. t <= t0+real(nt-1,fp)*dt ) then
      itb = max(it-nlag,1)
      ite = min(it+nlag,nt)
      do it=itb,ite
         aux1 = aux0 - real(it-1,fp)
         interp_trace_dp = interp_trace_dp + real(trace(it),fp)*lanczos(aux1,nlag)
      end do
   end if
   end function 
!
function bessel_i0(x)
real(fp), intent(in) :: x
real(fp) :: bessel_i0

integer, parameter :: np=7
integer, parameter :: nq=9
real(fp), dimension(np), parameter :: p(1:np)=(/1.0000000d0,  &
                                                3.5156229d0,  &
                                                3.0899424d0,  &
                                                1.2067492d0,  &
                                                0.2659732d0,  &
                                                0.3607680d-1, &
                                                0.4581300d-2/)

real(fp), dimension(nq), parameter :: q(1:nq)=(/0.39894228d0,  &
                                                0.13285920d-1, &
                                                0.22531900d-2, &
                                               -0.15756500d-2, &
                                                0.91628100d-2, &
                                               -0.20577060d-1, &
                                                0.26355370d-1, &
                                               -0.16476330d-1, &
                                                0.39237700d-2/)

real(fp), parameter :: xcut=3.75d0

real(fp) :: y,ax
ax=abs(x)
if ( ax < xcut ) then
y=(ax/xcut)**2
bessel_i0=p(1)+y*(p(2)+y*(p(3)+y*(p(4)+y*(p(5)+y*(p(6)+y*p(7))))))
else
y=xcut/ax
bessel_i0=(exp(ax)/sqrt(ax))*(q(1)+y*(q(2)+y*(q(3)+y*(q(4)+y*(q(5)+y*(q(6)+y*(q(7)+y*(q(8)+y*q(9)))))))))
end if

end function

      function wkeiser(x,lx,b)
! interpolation of spatial source distribution 
! Keiser interpolation kernel
! K(x,h) = sinc(x/h)*keiser(x,h)
   !
   implicit none
   integer,parameter :: fp=kind(1.d0)
   real(fp),intent(in)  :: x
   real(fp), intent(in) :: lx
   real(fp), intent(in) :: b
   real(fp) :: wkeiser
   real(fp) :: a
   real(fp), parameter :: one=1.d0
   !
   a = x/lx
   if ( abs(a) <= one  ) then
      wkeiser = sinc(a)*bessel_i0(b*sqrt(one-a*a))/ bessel_i0(b)
   else
      wkeiser = 0.d0
   end if
   end function wkeiser

!
! band limited spatial delta function using Keiser window
  subroutine pointsource_wkeiser(this,xs1,xs2)
  class(afd2d),intent(inout) :: this
  real(fp),intent(in) :: xs1,xs2
  !
  ! locals:
  real(fp) :: s1,s2,x1,x2,ds1,ds2
  integer  :: i1,i2,j1s,j2s,shift,iconv,iexit
  !
  integer, parameter :: nkeiser4=4
  real(fp),parameter :: bkeiser4=6.31d0
  !
  integer, parameter :: nkeiser8=8
  real(fp),parameter :: bkeiser8=12.52d0
  !
  !
  integer  :: nkeiser
  real(fp) :: hwidth 
  real(fp) :: bkeiser
  !
     hwidth=2.d0*min(this%d1,this%d2)

     s1=(xs1-this%o1)/this%d1+real(nborda,fp)
     s2=(xs2-this%o2)/this%d2+real(nborda,fp)

     j1s=floor(s1)
     j2s=floor(s2)

     ds1=s1-real(j1s,fp)
     ds2=s2-real(j2s,fp)

     this%wsrc1(:)  = 0.d0
     this%wsrc2(:)  = 0.d0
    
     bkeiser=bkeiser4
     nkeiser=floor(hwidth/this%d1)+1

     do iconv=-nkeiser,nkeiser
        shift=max(min(j1s+iconv,this%n1+nborda+nkeiser),nborda-nkeiser)
        x1= (real(iconv,fp)+ds1)*this%d1
        this%wsrc1(shift)=wkeiser(x1,hwidth,bkeiser) !/this%d1
     end do

     nkeiser=floor(hwidth/this%d2)+1
     do iconv=-nkeiser,nkeiser
        shift=max(min(j2s+iconv,this%n2+nborda+nkeiser),nborda-nkeiser+1)
        x2= (real(iconv,fp)+ds2)*this%d2
        this%wsrc2(shift)=wkeiser(x2,hwidth,bkeiser) !/this%d2
     end do
!
  end subroutine 

 !
! band limited spatial delta function using Keiser window
  subroutine pointreceiver_wkeiser(this,ircv,xr1,xr2)
  class(afd2d),intent(inout) :: this
  integer,intent(in) :: ircv
  real(fp),intent(in) :: xr1,xr2
  !
  ! locals:
  real(fp) :: r1,r2,x1,x2,dr1,dr2
  integer  :: i1,i2,j1,j2,shift,iconv,iexit
  !
  integer, parameter :: nkeiser4=4
  real(fp),parameter :: bkeiser4=6.31d0
  !
  integer, parameter :: nkeiser8=8
  real(fp),parameter :: bkeiser8=12.52d0
  !
  !
  integer  :: nkeiser
  real(fp) :: hwidth 
  real(fp) :: bkeiser
  !
     hwidth=2.d0*min(this%d1,this%d2)

     r1=(xr1-this%o1)/this%d1+real(nborda,fp)
     r2=(xr2-this%o2)/this%d2+real(nborda,fp)

     j1=floor(r1)
     j2=floor(r2)

     dr1=r1-real(j1,fp)
     dr2=r2-real(j2,fp)

     this%wrcv1(:,ircv)  = 0.d0
     this%wrcv2(:,ircv)  = 0.d0
    
     bkeiser=bkeiser4
     nkeiser=floor(hwidth/this%d1)+1

     do iconv=-nkeiser,nkeiser
        shift=max(min(j1+iconv,this%n1+nborda),nborda+1)
        x1= (real(iconv,fp)+dr1)*this%d1
        this%wrcv1(shift,ircv)=wkeiser(x1,hwidth,bkeiser) !/this%d1
     end do

     nkeiser=floor(hwidth/this%d2)+1
     do iconv=-nkeiser,nkeiser
        shift=max(min(j2+iconv,this%n2+nborda),nborda+1)
        x2= (real(iconv,fp)+dr2)*this%d2
        this%wrcv2(shift,ircv)=wkeiser(x2,hwidth,bkeiser) !/this%d2
     end do
!
  end subroutine 
 
!
! procedure: input_model(this)
! description
!
! read velocity model
! extends model boundaries
! allocate fields 
! 
! 
!
! procedure: input_data(this)
! description
!
! read data file
! allocate containers to read common shot
! check if acquisition fits into velocity model dimensions
! 
! procedure: set_fd_parameters(this)
! check dispersion and stability
! compute time evolution step
! compute ABC auxiliary parameters
!
! procedure: initial_conditions(this)
! set all fields to zero
!
! procedure: wfield_update(this)
! compute wavefield evolution
!
! procedure: modeling(this,adj)
! if adj=0 performs forward simulation
! if adj=1 performs reverse time simulation

! procedure: born_modeling(this,adj)
! if adj=0 performs forward born simulation
! if adj=1 performs adjoint born simulation

! procedure: cleanup(this)
! deallocate memory
! finalize MPI
!

subroutine input_model(this,velhdr,stat)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
character(*),intent(in)    :: velhdr
integer, intent(out)       :: stat

integer::MDLFID
character(500) :: velbin
integer  :: i1,i2,idir
real(sp) :: swap
logical :: any_failures=.false.

open(MDLFID, file=trim(velhdr),iostat=stat) 
if ( stat /= 0 ) then
     stat=1
     any_failures=.true.
     call assert([.not.any_failures],[error_message('Error opening velocity header file '//trim(velhdr))])
     return
end if
!
read(MDLFID,*) this%n1,this%n2,this%o1,this%o2,this%d1,this%d2
read(MDLFID,'(a)') velbin
close(MDLFID)

! 
allocate ( this%vel(this%n1+2*nborda,this%n2+2*nborda),                 &   ! alocacao de memoria na malha incluindo bordas de absorcao
           this%pfield1(this%n1+2*nborda,this%n2+2*nborda),             &   
           this%pfield2(this%n1+2*nborda,this%n2+2*nborda),             &   ! 
           this%gamma1(this%n1+2*nborda),this%gamma2(this%n2+2*nborda), &   !
           this%gradwind1(this%n1+2*nborda), &
           this%gradwind2(this%n2+2*nborda), &
           this%wsrc1(this%n1+2*nborda),this%wsrc2(this%n2+2*nborda),stat=stat)
if ( stat /= 0 ) then 
     stat=2
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('Allocation failed for AFD2D type')])
     return
end if
!
this%p1=>this%pfield1(:,:)
this%p2=>this%pfield2(:,:)

! acessa o arquivo de velocidades
idir=max(index(velhdr,'\',.true.), index(velhdr,'/',.true.))
velbin=velhdr(1:idir)//trim(velbin) ! adiciona caminho do diretorio
open(MDLFID,file=trim(velbin), form='unformatted',   &
        access='stream',status='old',action='read',iostat=stat)
if ( stat /= 0 ) then
     stat=3
     any_failures=.true.
     call assert([.not.any_failures], &
          [error_message('Error opening velocity binary file '//trim(velbin))])
     return
end if
!
! leitura do modelo de velocidade
  do i2=nborda+1,this%n2+nborda
     do i1=nborda+1,this%n1+nborda
        read(MDLFID) swap
        this%vel(i1,i2) = real(swap,dp) 
     end do
  end do
close(MDLFID) !  
!
! model limits:
this%i1min=1+nborda
this%i1max=this%n1+nborda
this%i2min=1+nborda
this%i2max=this%n2+nborda
!
! modeling grid limits
this%i1beg=max(this%i1min-WDWSRC-nborda,1)
this%i1end=min(this%i1max+WDWSRC+nborda,this%n1+2*nborda)
this%i2beg=max(this%i2min-WDWSRC-nborda,1)
this%i2end=min(this%i2max+WDWSRC+nborda,this%n2+2*nborda)
!
! calcula os vmax e vmin 
this%vmax=this%vel(nborda+1,nborda+1)                 !  
this%vmin=this%vel(nborda+1,nborda+1)                 !
do i2=nborda+1,this%n2+nborda
   do i1=nborda+1,this%n1+nborda
      this%vmin=min(this%vmin,this%vel(i1,i2))        !   
      this%vmax=max(this%vmax,this%vel(i1,i2))        !  
   end do
end do
!
! preenche o modelo de velocidades nas bordas 
!
do i2=1,nborda
   do i1=nborda+1,nborda+this%n1
      this%vel(i1,i2)                = this%vel(i1,nborda+1)
      this%vel(i1,i2+nborda+this%n2) = this%vel(i1,nborda+this%n2)
   end do
end do
!
do i2=1,this%n2+2*nborda
   do i1=1,nborda
      this%vel(i1,i2)=this%vel(nborda+1,i2)
   end do
   do i1=nborda+this%n1+1,this%n1+2*nborda
      this%vel(i1,i2)=this%vel(nborda+this%n1,i2)
   end do
end do
end subroutine input_model

subroutine input_regular_data(this,ns,o1s,o2s,ds, &
                                   nr,o1r,o2r,dr, &
                                   nt,ot,dt,scalco,scalel,datafile)
!
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
!
integer, intent(in) :: ns  ! number of sources
real(sp), intent(in) :: o1s ! sources depth
real(sp), intent(in) :: o2s ! first source horizontal coordinate
real(sp), intent(in) :: ds  ! distance between sources 
integer, intent(in) :: nr  ! number of receivers per shot
real(sp), intent(in) :: o1r ! receivers depth
real(sp), intent(in) :: o2r ! minimum offset
real(sp), intent(in) :: dr  ! distance between receivers
integer, intent(in) :: nt  ! number of samples per trace
real(sp), intent(in) :: ot  ! initial recording time
real(sp), intent(in) :: dt  ! time sampling interval
!
real(sp), intent(inout) :: scalco
real(sp), intent(inout) :: scalel
!
character(*),intent(inout) :: datafile
!
integer :: ntrc,is,ir
!
integer :: iexit,stat_alloc
!
type(su_trace) :: trace
!
real(sp) :: xs,xr
integer :: nsamples,itrc
!
logical :: any_failures=.false.

ntrc = ns*nr

this%datafile=trim(datafile)
this%nshots=ns
this%ngmax=nr
allocate(this%wdata(this%nshots),stat=stat_alloc)
if ( stat_alloc /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('Allocation error for wdata(:) at procedure input_regular_data')])
end if
allocate(this%shotid(this%nshots+1),stat=stat_alloc)
if ( stat_alloc /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('Allocation error for shotid(:) at procedure input_regular_data ')])
end if

this%wdata(:)=1.d0

trace%trcheader%ns=nt
trace%trcheader%dt=int(dt*sec2ms,2)

this%ntrec  = nt
this%t0rec  = ot
this%dtrec  = dt! intervalo de amostragem dos dados
this%ttotal = ot+real(nt-1)*this%dtrec

allocate(trace%trcdata(nt),stat=stat_alloc)
trace%trcdata(1:nt)=0.0
!
trace%trcheader%scalco=int(scalco)
trace%trcheader%scalel=int(scalel)
if ( scalco < 0.0  ) scalco = 1.0/abs(scalco)
if ( scalel < 0.0 ) scalel = 1.0/abs(scalel)
!
iexit =  OpenSU(SUDATFID,trim(this%datafile),'w',nsamples)
if ( iexit /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('Error opening seismic data file '//trim(this%datafile))])
end if
!
ntrc=0
do is=1,ns
   this%shotid(is)=ntrc+1
   xs = o2s+real(is-1)*ds
   trace%trcheader%fldr   = is
   trace%trcheader%sx     = int(scalco*xs)
   trace%trcheader%sdepth = int(scalel*o1s)
   do ir=1,nr
      itrc=itrc+1
      xr = xs+o2r+real(ir-1)*dr
      trace%trcheader%tracl  = itrc
      trace%trcheader%tracr  = itrc
      trace%trcheader%tracf  = ir
      trace%trcheader%gx     = int(scalco*xr)
      trace%trcheader%gelev  =-int(scalel*o1r)
      iexit=PutSUTrace(SUDATFID,itrc,trace)
   end do
end do
this%shotid(ns+1)=ntrc+1
close(SUDATFID)
end subroutine


subroutine input_data(this,sufile,stat)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
character(*),intent(in)    :: sufile
integer, intent(out)       :: stat
!
type(su_trace) :: trace
real(dp) :: scalco,scalel
integer  :: ntrc,is
integer  :: itrc,isrc
logical :: any_failures=.false.
real    :: dnorm
!

this%datafile=trim(sufile)


! ler arquivo common-shot:
stat =  OpenSU(SUDATFID,trim(this%datafile),'r',this%ntrec)
if ( stat /= 0 ) then
     stat = 1
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('Error opening seismic data file '//trim(this%datafile))])
     return 
end if
!
!
allocate(trace%trcdata(this%ntrec), stat=stat) 
if ( stat /= 0 ) then
     stat = 2
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('Allocation error for trace at the procedure input_data')] )
     return
end if
!
! conta numero de tracos no gather
stat = GetSUTrace(SUDATFID,1,trace)
if ( stat /= 0 ) then
     stat = 3
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('error reading trace at the procedure input_data')] )
     return
end if
!
this%dtrec  = real(trace%trcheader%dt)/sec2mcs ! intervalo de amostragem dos dados
this%ttotal = real(trace%trcheader%ns-1)*this%dtrec
!
! fator de escala para coordenadas scalco
scalco = 1.0
if ( trace%trcheader%scalco < 0 ) then
     scalco = 1.0 / abs(real(trace%trcheader%scalco)) !SEG-Y standard
else if ( trace%trcheader%scalco > 0 ) then
     scalco = real( trace%trcheader%scalco) ! SEG-Y standard
end if
!
! fator de escala para elevacoes e profundidade de fontes e receptores scalel
scalel = 1.0
if ( trace%trcheader%scalel < 0 ) then
     scalel = 1.0 / abs(real(trace%trcheader%scalel)) !SEG-Y standard
else if ( trace%trcheader%scalel > 0 ) then
     scalel = real( trace%trcheader%scalel) ! SEG-Y standard
end if
!
this%nshots = 1
ntrc   = 1 
isrc     = trace%trcheader%sx
loop_data_headers:DO 
   stat = GetSUTrace(SUDATFID,ntrc,trace)
!
   IF( stat == 0 ) THEN
      ntrc = ntrc + 1
      if( trace%trcheader%sx /= isrc ) then
          this%nshots = this%nshots + 1
          isrc    = trace%trcheader%sx
       end if
   ELSE
       EXIT loop_data_headers 
    ENDIF 
!
ENDDO loop_data_headers


allocate(this%wdata(this%nshots),stat=stat)
if ( stat /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('allocation failed for wdata(:) at the procedure input_data')] )
end if
!
allocate(this%shotid(this%nshots+1),stat=stat)
if ( stat /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('allocation failed for shotid(:) at the procedure input_data')] )
end if
!
stat        = GetSUTrace(SUDATFID,1,trace)
isrc        = trace%trcheader%sx
this%nshots = 1
this%shotid(1)  = 1
ntrc        = 1 
this%wdata(:)=0.d0
DO 
   stat = GetSUTrace(SUDATFID,ntrc,trace)
!
   IF( stat == 0 ) THEN
      ntrc = ntrc + 1
      if( trace%trcheader%sx /= isrc ) then
          this%nshots = this%nshots + 1
          isrc    = trace%trcheader%sx
          this%shotid(this%nshots) = ntrc
      else
          this%wdata(this%nshots)=this%wdata(this%nshots) + &
                  real(dot_product(trace%trcdata,trace%trcdata),dp)
      end if
   ELSE
       EXIT 
    ENDIF 
ENDDO
this%shotid(this%nshots+1) = ntrc
!
allocate(this%survey(ntrc-1),STAT=stat)
if ( stat /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('allocation failed for survey(:) at  procedure input_data')] )
     return
end if
!
dnorm=0
do itrc=1,ntrc-1
   allocate(this%survey(itrc)%trcdata(this%ntrec),STAT=stat)
if ( stat /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('allocation failed for survey(:)%trcdata(:) at  procedure input_data')] )
     return
end if
   stat = GetSUTrace(SUDATFID,itrc,this%survey(itrc))
if ( stat /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('error reading data at procedure input_data')] )
     return
end if
   dnorm=dnorm+ &
   dot_product(this%survey(itrc)%trcdata,this%survey(itrc)%trcdata)
end do
CLOSE(SUDATFID)
!
print *,'myrank ',this%mpimyrank,dnorm
!
do is=1,this%nshots
!if ( this%wdata(is) > 0.d0 ) then
!this%wdata(is)=1.d0/sqrt(this%wdata(is))
!end if
this%wdata(is)=1.d0
end do
this%ngmax=0
do isrc=1,this%nshots
   this%ngmax = max(this%ngmax,this%shotid(isrc+1)-this%shotid(isrc))
end do
!
allocate(this%shotgather(this%ngmax),stat=stat)
if ( stat /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('allocation failed for shotgather(:) at  procedure input_data')] )
     return
end if
!
do itrc=1,this%ngmax
   allocate(this%shotgather(itrc)%trcdata(this%ntrec),STAT=stat)
   if ( stat /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('allocation failed for shotgather(:)%trcdata(:) at procedure input_data')] )
     return
   end if
end do
!
allocate(this%g1(this%ngmax),this%g2(this%ngmax),stat=stat)
if ( stat /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('allocation failed for g1(:), g2(:) at procedure input_data')] )
end if
!
   allocate(this%wrcv1(this%n1+2*nborda,this%ngmax),stat=stat)
    if( stat /= 0 ) then
       !print *, 'shot index ',shotindex,iexit
       any_failures=.true.
       call assert([.not.any_failures], &
       [error_message('Error allocating wrcv1(:,:)!')])
    endif

   allocate(this%wrcv2(this%n2+2*nborda,this%ngmax),stat=stat)
    if( stat /= 0 ) then
       !print *, 'shot index ',shotindex,iexit
       any_failures=.true.
       call assert([.not.any_failures], &
       [error_message('Error allocating wrcv2(:,:)!')])
    endif


end subroutine

subroutine time_update(this)
class(afd2d),intent(inout) :: this

integer  :: i1,i2,ic
real(dp) :: laplacian,pderiv1,pderiv2,pmlderiv1,pmlderiv2
real(dp) :: gama,invgama,mgama

this%fswap=>null()

!******************************************************************************************************************
! atualizacao do campo de onda:
!
do i2=this%i2beg+DRVLEN,this%i2end-DRVLEN
do i1=this%i1beg+DRVLEN,this%i1end-DRVLEN

   gama      = this%gamma1(i1)+this%gamma2(i2)
   !
   ! calculo do laplaciano:
   laplacian= 2.d0*deriv2(1)*this%p2(i1,i2)
   do ic=1,DRVLEN-1
      laplacian =  laplacian +  & 
                   deriv2(ic+1)*( this%p2(i1-ic,i2)+this%p2(i1+ic,i2) + &
                                  this%p2(i1,i2-ic)+this%p2(i1,i2+ic) ) 
   end do
   !
   ! FD update:
   this%p1(i1,i2) = (2.d0-gama)*this%p2(i1,i2)-(1.d0-gama)*this%p1(i1,i2) + & 
                          this%vel(i1,i2)*laplacian     
    end do
end do
!
! permutacao de arranjos no time:
this%fswap=>this%p2
this%p2=>this%p1
this%p1=>this%fswap

end subroutine


subroutine set_source_pulse(this)
class(afd2d),intent(inout) :: this
!
real(dp) :: beta,aux0,aux
integer  :: it,stat
!
allocate(this%source_pulse(this%ntrec),stat=stat)
!
beta=pi*this%freq*this%dtrec
aux0=2.d0/(this%freq*this%dtrec)
do it=1,this%ntrec
   aux = beta*(real(it-1,dp)-aux0)
   this%source_pulse(it) = exp(-aux*aux)*(1.d0-2.d0*aux*aux)
end do
!
!print *, 'pulse ',mpimyrank,dot_product(this%source_pulse,this%source_pulse) 
!
end subroutine

subroutine source_injection(this,adj,t)
class(afd2d),intent(inout) :: this
integer, intent(in) :: adj
real(dp), intent(in) :: t
!
real(dp) :: pulse, source
integer :: i1,i2,ig
!

if ( adj == 0 ) then

   pulse=interp_trace(t,this%ntrec,0.d0,this%dtrec,this%source_pulse)
!
   source=0.d0
   do i2=this%i2min,this%i2max
   do i1=this%i1min-WDWSRC,this%i1max+WDWSRC
      this%p2(i1,i2) = this%p2(i1,i2) + this%vel(i1,i2)*pulse*this%wsrc1(i1)*this%wsrc2(i2)
      source=source+this%p2(i1,i2)*this%p2(i1,i2)
   end do
   end do
   !print *,'myrank = ',mpimyrank, source

else

   do ig=1,this%ng
      pulse=interp_trace(t,this%ntrec,0.d0,this%dtrec,this%shotgather(ig)%trcdata)
!
      do i2=this%i2min,this%i2max
      do i1=this%i1min-WDWSRC,this%i1max+WDWSRC
      ! corrigir 
         this%p2(i1,i2) = this%p2(i1,i2) + &
         this%vel(i1,i2)*pulse*this%wrcv1(i1,ig)*this%wrcv2(i2,ig)
      end do
      end do
    end do
end if

end subroutine

subroutine set_shots(this,shot_beg,shot_end,shot_step)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
integer, intent(in) :: shot_beg,shot_end,shot_step
logical :: any_failures=.false.

this%ishot_beg=min(max(shot_beg,1),this%nshots)
this%ishot_end=max(min(shot_end,this%nshots),1)
this%step_shots=max(shot_step,1)

end subroutine



subroutine close_cshot_file()
    close(SUDATFID)
end subroutine 




subroutine get_cshot_record(this,shotindex,datafile,iexit)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
integer, intent(in) :: shotindex
character(:),allocatable, intent(inout) :: datafile
integer, intent(out)::iexit

! local
integer :: itrc,ig,ns
real(dp) :: scalco,scalel,s1,s2,shotspread,x1min,x1max,x2min,x2max
logical :: any_failures=.false.
logical :: fidopen


   this%ng = this%shotid(shotindex+1)-this%shotid(shotindex)

    this%wrcv1(:,:)=0.d0
    this%wrcv2(:,:)=0.d0

   ig = 0
   shotspread=0.d0
   loop_cshot_gather:do itrc=this%shotid(shotindex), this%shotid(shotindex+1)-1
      ig = ig + 1
      this%shotgather(ig)%trcheader=this%survey(itrc)%trcheader
      this%shotgather(ig)%trcdata(:)=this%survey(itrc)%trcdata(:)
      if ( ig == 1 ) then
           scalco=real(this%shotgather(ig)%trcheader%scalco,dp)
           if ( scalco < 0.d0 ) then
           scalco=1.d0/abs(scalco)
           else if ( scalco < 1.d0 ) then
              scalco=1.d0
           end if
           scalel=real(this%shotgather(ig)%trcheader%scalel,dp)
           if ( scalel < 0.d0 ) then
                scalel=1.d0/abs(scalel)
           else if ( scalel < 1.d0 ) then
                scalel=1.d0
           end if
           s1 =  real(this%shotgather(ig)%trcheader%sdepth,dp)*scalel
           s2 =  real(this%shotgather(ig)%trcheader%sx,dp)*scalco
      end if
!
!  armazena posicao do traco ma malha do modelo
!  converte unidades usando as palavras scalco e scalel do header do traco
!
      this%g1(ig) = -real(this%shotgather(ig)%trcheader%gelev,dp)*scalel
      this%g2(ig) =  real(this%shotgather(ig)%trcheader%gx,dp)*scalco
      shotspread  =  max(shotspread,abs(this%g2(ig)-s2))

      call pointreceiver_wkeiser(this,ig,this%g1(ig),this%g2(ig))
      ! apply data weight
      this%shotgather(ig)%trcdata(:)=this%wdata(shotindex)*this%shotgather(ig)%trcdata(:)

      this%datanorm=this%datanorm + &
                  real(dot_product(this%shotgather(ig)%trcdata,this%shotgather(ig)%trcdata),dp)
   end do   loop_cshot_gather

! definir janela lateral de modelagem
   x2min  = max(min(minval(this%g2(1:this%ng)),s2)-shotspread,this%o2)
   x2max  = min(max(maxval(this%g2(1:this%ng)),s2)+shotspread, &
                this%o2+real(this%n2-1)*this%d2)
!
! definir janela vertical para computar gradiente
   x1min=max(maxval(this%g1(1:this%ng)),s1)+2.d0*minval(this%vel(nborda+1,nborda+1:nborda+this%n2))/this%freq
   x1max=this%o1+real(this%n1-1,dp)*this%d1
!
   call pointsource_wkeiser(this,s1,s2)

   end subroutine


subroutine put_cshot_record(this,shotindex,datafile,iexit)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
integer, intent(in) :: shotindex
character(:),allocatable, intent(inout) :: datafile
integer, intent(out)::iexit

! local
integer :: FID
integer :: itrc,ig,ns
real(dp) :: scalco,scalel,s1,s2,shotspread,x2min,x2max
logical :: any_failures=.false.

   this%ng = this%shotid(shotindex+1)-this%shotid(shotindex)
   ig = 0
   shotspread=0.d0
   loop_cshot_gather:do itrc=this%shotid(shotindex), this%shotid(shotindex+1)-1
      ig = ig + 1
      !
      ! apply data weight
      !this%shotgather(ig)%trcdata(:)=this%wdata(shotindex)*this%shotgather(ig)%trcdata(:)
      !this%shotgather(ig)%trcdata(:)=this%shotgather(ig)%trcdata(:)
      this%datanorm=this%datanorm + &
                  real(dot_product(this%shotgather(ig)%trcdata,this%shotgather(ig)%trcdata),dp)
      this%survey(itrc)%trcdata(:)=this%shotgather(ig)%trcdata(:)

   end do   loop_cshot_gather

   end subroutine


subroutine open_cshot_file(this,datafile,ioaction,iexit)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
character(:),allocatable, intent(inout) :: datafile
character(1),intent(in) :: ioaction
integer, intent(out)::iexit

logical :: any_failures=.false.

select case (ioaction)
case('r')
   iexit =  OpenSU(SUDATFID,trim(datafile),'r',this%ntrec)
   if( iexit /= 0 ) then
       !print *, 'shot index ',shotindex,iexit
       iexit = 1
       any_failures=.true.
       call assert([.not.any_failures], &
       [error_message('Error opening common-shot file to read '//trim(datafile))])
       return
   endif
case('a')
   iexit =  OpenSU(SUDATFID,trim(datafile),'a',this%ntrec)
   if( iexit /= 0 ) then
       !print *, 'shot index ',shotindex,iexit
       iexit = 1
       any_failures=.true.
       call assert([.not.any_failures], &
       [error_message('Error opening common-shot file to read '//trim(datafile))])
       return
   endif
end select

end subroutine


subroutine read_cshot_record(this,shotindex,datafile,iexit)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
integer, intent(in) :: shotindex
character(:),allocatable, intent(inout) :: datafile
integer, intent(out)::iexit

! local
integer :: itrc,ig,ns
real(dp) :: scalco,scalel,s1,s2,shotspread,x1min,x1max,x2min,x2max
logical :: any_failures=.false.
logical :: fidopen


   this%ng = this%shotid(shotindex+1)-this%shotid(shotindex)

    this%wrcv1(:,:)=0.d0
    this%wrcv2(:,:)=0.d0


   call open_cshot_file(this,datafile,'r',iexit)
   ig = 0
   shotspread=0.d0
   loop_cshot_gather:do itrc=this%shotid(shotindex), this%shotid(shotindex+1)-1
      ig = ig + 1
      iexit = GetSUTrace(SUDATFID,itrc,this%shotgather(ig))
      if( iexit /= 0 ) then
          iexit = 2
       any_failures=.true.
       call assert([.not.any_failures], &
       [error_message('Error reading common-shot trace '//trim(datafile))])
          return
      endif

      if ( ig == 1 ) then
           scalco=real(this%shotgather(ig)%trcheader%scalco,dp)
           if ( scalco < 0.d0 ) then
           scalco=1.d0/abs(scalco)
           else if ( scalco < 1.d0 ) then
              scalco=1.d0
           end if
           scalel=real(this%shotgather(ig)%trcheader%scalel,dp)
           if ( scalel < 0.d0 ) then
                scalel=1.d0/abs(scalel)
           else if ( scalel < 1.d0 ) then
                scalel=1.d0
           end if
           s1 =  real(this%shotgather(ig)%trcheader%sdepth,dp)*scalel
           s2 =  real(this%shotgather(ig)%trcheader%sx,dp)*scalco
      end if
!
!  armazena posicao do traco ma malha do modelo
!  converte unidades usando as palavras scalco e scalel do header do traco
!
      this%g1(ig) = -real(this%shotgather(ig)%trcheader%gelev,dp)*scalel
      this%g2(ig) =  real(this%shotgather(ig)%trcheader%gx,dp)*scalco
      shotspread  =  max(shotspread,abs(this%g2(ig)-s2))

      call pointreceiver_wkeiser(this,ig,this%g1(ig),this%g2(ig))
      ! apply data weight
      this%shotgather(ig)%trcdata(:)=this%wdata(shotindex)*this%shotgather(ig)%trcdata(:)

      this%datanorm=this%datanorm + &
                  real(dot_product(this%shotgather(ig)%trcdata,this%shotgather(ig)%trcdata),dp)
   end do   loop_cshot_gather
   
   call close_cshot_file()

! definir janela lateral de modelagem
   x2min  = max(min(minval(this%g2(1:this%ng)),s2)-shotspread,this%o2)
   x2max  = min(max(maxval(this%g2(1:this%ng)),s2)+shotspread, &
                this%o2+real(this%n2-1)*this%d2)
!
! definir janela vertical para computar gradiente
   x1min=max(maxval(this%g1(1:this%ng)),s1)+2.d0*minval(this%vel(nborda+1,nborda+1:nborda+this%n2))/this%freq
   x1max=this%o1+real(this%n1-1,dp)*this%d1
!
   call pointsource_wkeiser(this,s1,s2)

   end subroutine


subroutine write_cshot_record(this,shotindex,datafile,iexit)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
integer, intent(in) :: shotindex
character(:),allocatable, intent(inout) :: datafile
integer, intent(out)::iexit

! local
integer :: FID
integer :: itrc,ig,ns
real(dp) :: scalco,scalel,s1,s2,shotspread,x2min,x2max
logical :: any_failures=.false.

   this%ng = this%shotid(shotindex+1)-this%shotid(shotindex)
   call open_cshot_file(this,datafile,'a',iexit)
   ig = 0
   shotspread=0.d0
   loop_cshot_gather:do itrc=this%shotid(shotindex), this%shotid(shotindex+1)-1
      ig = ig + 1
      !
      ! apply data weight
      !this%shotgather(ig)%trcdata(:)=this%wdata(shotindex)*this%shotgather(ig)%trcdata(:)
      !this%shotgather(ig)%trcdata(:)=this%shotgather(ig)%trcdata(:)
      this%datanorm=this%datanorm + &
                  real(dot_product(this%shotgather(ig)%trcdata,this%shotgather(ig)%trcdata),dp)
      iexit = PutSUTrace(SUDATFID,itrc,this%shotgather(ig));
      if( iexit /= 0 ) then
          iexit = 2
       any_failures=.true.
       call assert([.not.any_failures], &
       [error_message('Error writing common-shot trace '//trim(datafile))])
          return
      endif

   end do   loop_cshot_gather
   
   call close_cshot_file()

   end subroutine


   subroutine open_moviefile(this,shotindex,ioaction)
   use assertion_mod, only:assert,error_message
   class(afd2d),intent(inout) :: this
   integer, intent(in) :: shotindex
   character(*), intent(in) :: ioaction

   logical :: any_failures=.false.
   integer :: reclen,ioexit

   reclen = (this%i1max-this%i1min+1)*(this%i2max-this%i2min+1)*SP

   this%MVFID=MVFID
   this%moviefile(1:9)='afdmovie_'
   write(this%moviefile(10:14),'(i5.5)') shotindex
   this%moviefile(15:18)='.bin'
   select case(ioaction)

   case('write')
   open(this%MVFID,file=trim(this%moviefile),  & 
          form='unformatted',     &
          access='direct',        &
          recl=RECLEN,            &
          status='replace',       &
          action=ioaction,        &
          iostat=ioexit)
   if( ioexit /= 0 ) then
       ioexit = 1
       any_failures=.true.
       call assert([.not.any_failures], &
       [error_message('Error opening file '//trim(this%moviefile)//' to write!')])
       return
   endif

  case('read')
   open(this%MVFID,file=trim(this%moviefile),  & 
          form='unformatted',     &
          access='direct',        &
          recl=RECLEN,            &
          status='old',           &
          action=ioaction,        &
          iostat=ioexit)
   if( ioexit /= 0 ) then
       ioexit = 1
       any_failures=.true.
       call assert([.not.any_failures], &
       [error_message('Error opening file '//trim(this%moviefile)//' to read!')])
       return
   endif
  end select
!
end subroutine

subroutine record_field_at_receivers(this,itrec)
class(afd2d),intent(inout) :: this
integer, intent(in)        :: itrec
!
!
integer :: i1,i2,ig,ig1
real(dp) :: recfield

! model = 0 , record the modeled field
! model = 1 , record the residuals: observed field - modeled field

select case(this%mode)
case(0)
   do ig=1,this%ng
      !
      recfield = 0.d0
      do i2=this%i2min,this%i2max
      do i1=this%i1min-WDWSRC,this%i1max+WDWSRC
      ! corrigir 
         recfield = recfield + &
         this%wrcv1(i1,ig)*this%wrcv2(i2,ig)*this%p2(i1,i2)
      end do
      end do

      this%shotgather(ig)%trcdata(itrec) = recfield

    end do
case(1)
   do ig=1,this%ng
      !
      recfield = 0.d0
      do i2=this%i2min,this%i2max
      do i1=this%i1min-WDWSRC,this%i1max+WDWSRC
      ! corrigir 
         recfield = recfield + &
         this%wrcv1(i1,ig)*this%wrcv2(i2,ig)*this%p2(i1,i2)
      end do
      end do

      this%shotgather(ig)%trcdata(itrec) = this%shotgather(ig)%trcdata(itrec)-recfield

    end do
end select
end subroutine


subroutine afd_solver(this,datafile,iexit)
class(afd2d),intent(inout) :: this
character(:),allocatable, intent(inout) :: datafile
integer, intent(out)  :: iexit
!
! local
integer  :: ishot,ig,itrc,ns,it,itrec,i1,i2,frameindex
real(dp) :: scalco,scalel,s1,s2,shotspread,datnorm
real(dp) :: t
real(sp) :: dnorm
integer :: ishot_beg
integer :: ishot_end

! forward modeling
!
loop_shots:do ishot=this%step_shots*this%mpimyrank+this%ishot_beg, &
                    this%ishot_end,this%step_shots*this%mpinproc
!
print *,'myrank :', this%mpimyrank, ishot
!
! initial conditions: 
this%p1(:,:) = 0.d0
this%p2(:,:) = 0.d0
!
itrec=0
loop_tforward: do it=1,this%nt                 ! evolucao do campo de pressao
!
t=real(it-1,dp)*this%dt
!
call source_injection(this,0,t)
!
call time_update(this)
!
if ( mod(it-1,this%ndtrec) == 0.0 ) then
     itrec=itrec+1
     call record_field_at_receivers(this,itrec)
end if
!
end do loop_tforward
!call write_cshot_record(this,ishot,this%datafile,iexit)
call put_cshot_record(this,ishot,datafile,iexit)
end do loop_shots

end subroutine



subroutine set_fd_solver(this,freq,iexit)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
real(dp), intent(in) :: freq
integer, intent(out) :: iexit
logical :: any_failures=.false.
!
!
real(dp) :: mu,alpha,gama,dxmax
integer  :: i1,i2,ndt,ntsave
! verifica condicoes de  estabilidade e dispersao numerica
!
this%freq=freq
dxmax       = this%vmin / (6.d0*this%freq)
this%dt     = cfl *sqrt(2.d0/sum(abs(deriv2(:))))*this%d1/this%vmax
this%ndtrec = floor(this%dtrec/this%dt) + 1      
!
! intervalo de tempo de registro multiplo inteiro
! do intervalo de modelagem  dt
!
if ( this%ndtrec > 1) then
   this%dt = this%dtrec/real(this%ndtrec)
else
   this%dt = this%dtrec
end if
!
! intervalo para salvar campos de onda
this%ndtsave = floor(dtsave/this%dt) + 1      

this%nt=ceiling(this%ttotal/this%dt)  

!ntsave=this%nt/this%ndtsave

! make modeling time a multiple of frame saving time interval
!this%nt=(ntsave+1)*this%ndtsave

if( this%d1 > dxmax ) then
    iexit=1
       any_failures=.true.
       call assert([.not.any_failures], &
       [error_message('Error stability conditions failed !')])
   return
else
   iexit=0
end if
!
! velocidade de propagacao normalizada: 
!
mu = this%dt/this%d1
do i2=1,this%n2+2*nborda
   do i1=1,this%n1+2*nborda
      this%vel(i1,i2)=(this%vel(i1,i2)*mu)**2
   end do
end do
!
! mascara de absorcao:
alpha = -0.75d0*log(R0pml)*this%dt*this%vmax/(real(Nborda,dp)*this%d1)
this%gamma1(:) = 0.d0
this%gamma2(:) = 0.d0
do i2=1,nborda
   gama = abs(real(i2,dp)/real(nborda,dp))
   gama = alpha*(gama-sin(twopi*gama)/twopi)
   this%gamma2(nborda+1-i2) = gama
   this%gamma2(this%n2+nborda+i2) = gama
end do

do i1=1,nborda-WDWSRC
   gama = abs(real(i1,dp)/real(nborda,dp))
   gama = alpha*(gama-sin(twopi*gama)/twopi)
   this%gamma1(nborda+1-WDWSRC-i1) = gama
   this%gamma1(this%n1+nborda+WDWSRC+i1) = gama
end do
call set_source_pulse(this)

end subroutine


!          
subroutine write_survey(this,iexit)
use assertion_mod, only:assert,error_message
class(afd2d),intent(inout) :: this
integer, intent(out) :: iexit
!
integer :: ishot,itrc
logical  :: any_failures=.false.
!
iexit =  OpenSU(SUDATFID,trim(this%datafile),'a',this%ntrec)
if ( iexit /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('Error opening seismic data file '//trim(this%datafile))])
end if

loop_shots:do ishot=this%step_shots*this%mpimyrank+this%ishot_beg, &
                    this%ishot_end,this%step_shots*this%mpinproc
!
   if ( ishot <= this%ishot_end ) then
        do itrc=this%shotid(ishot),this%shotid(ishot+1)-1
           iexit=PutSUTrace(SUDATFID,itrc,this%survey(itrc))
!
if ( iexit /= 0 ) then
     any_failures=.true.
     call assert([.not.any_failures], &
     [error_message('Error writing modeled seismic data file '//trim(this%datafile))])
end if

        end do
   else
        exit loop_shots
   end if
end do loop_shots

close(SUDATFID)
 

end subroutine
!
subroutine afd_solver_driver()
!
type(afd2d)::afd
!
integer,parameter :: NAMELEN=1024
integer,parameter :: PARFID=72
!
integer::ns,nr,nt
real(sp) :: o1s,o2s,ds
real(sp) :: o1r,o2r,dr
real(sp) :: ot,dt
real(sp) :: scalco=100.0
real(sp) :: scalel=100.0
!
character(NAMELEN) :: velhdr,sudata
!
real(dp) :: freq
!
integer :: iexit,stat
!
character(NAMELEN) :: parfile
!
real(dp),allocatable :: grad(:,:)
real(dp) :: grad_norm,data_norm
!
call mpi_begin(afd,iexit)
!
!entrada do modelo de velcidades
if ( afd%mpimyrank == afd%mpimaster ) then
     read(*,'(a)') parfile
end  if
!
call mpi_bcast(parfile,NAMELEN,MPI_CHARACTER,afd%mpimaster,MPI_COMM_WORLD,iexit)
!
open(PARFID,file=trim(parfile),status='old',action='read',iostat=iexit)
!
!write(*,'("Entre vel. model header file: ")',advance='no')
read(PARFID,'(A)') velhdr
!
!write(*,'("Entre freq. pulso fonte: ")',advance="no")
read(PARFID,*) freq
!
read(PARFID,*) ns,o1s,o2s,ds
read(PARFID,*) nr,o1r,o2r,dr
read(PARFID,*) nt,ot,dt
!
!write(*,'("Entre SU data file: ")',advance='no')
read(PARFID,'(A)') sudata
!
close(PARFID)
!
afd%mode=0
!
if ( afd%mpimyrank == afd%mpimaster ) then
write(*,'("INPUT DATA :")')
write(*,'("Vel. model : ",a)') trim(velhdr)
write(*,'("Data file: ",a)') trim(sudata)
write(*,'("Freq. pulso fonte: ",f7.2," Hz")') freq
end if
!
call input_model(afd,velhdr,iexit)
!
call input_regular_data(afd,ns,o1s,o2s,ds, &
                            nr,o1r,o2r,dr, &
                            nt,ot,dt,scalco,scalel,sudata)
!
call set_fd_solver(afd,freq,iexit)
print *,'set_fd_solver done'
!
!
call  write_survey(afd,iexit)
!
call mpi_end(afd,iexit)
!
end subroutine


end module afd2d_mod 
