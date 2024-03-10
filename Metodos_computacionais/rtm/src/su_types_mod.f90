MODULE SU_Types_mod
  ! 
  ! Data types for SU Headers FORMAT: SU Trace_Header 
  ! 
  ! OBS : Extended Headers formats are not implemented yet, 
  !       but can be added using specific data types as bellow 
  ! 
  ! References: 
  ! 
  ! 
  ! Seismic Unix Manual 
  ! 
  ! 
  ! Author : Jesse Costa  
  ! e-mail : jesse@ufpa.br 
  ! 
  !
  INTEGER,PRIVATE,PARAMETER :: SP=KIND(1.0)
  INTEGER,PRIVATE,PARAMETER :: I4B  = SELECTED_INT_KIND(9)
  INTEGER,PRIVATE,PARAMETER :: I2B  = SELECTED_INT_KIND(4)
  INTEGER,PRIVATE,PARAMETER :: I1B  = SELECTED_INT_KIND(2)

  REAL(SP)  ,PRIVATE,PARAMETER :: NR32 = 0.0 
  INTEGER(4),PRIVATE,PARAMETER :: NA32 = 0 
  INTEGER(2),PRIVATE,PARAMETER :: NA16 = 0 

  TYPE SU_TraceHeader 
     SEQUENCE 
     INTEGER(4):: tracl                                  = NA32     !  1-4    
     INTEGER(4):: tracr                                  = NA32     !  5-8  
     INTEGER(4):: fldr                                   = NA32     !  9-12 
     INTEGER(4):: tracf                                  = NA32     ! 13-16 
     INTEGER(4):: ep                                     = NA32     ! 17-20 
     INTEGER(4):: cdp                                    = NA32     ! 21-24 
     INTEGER(4):: cdpt                                   = NA32     ! 25-28 
     INTEGER(2):: trid                                   = NA16     ! 29-30 
     INTEGER(2):: nvs                                    = NA16     ! 31-32 
     INTEGER(2):: nhs                                    = NA16     ! 33-34 
     INTEGER(2):: duse                                   = NA16     ! 35-36 
     INTEGER(4):: offset                                 = NA32     ! 37-40 
     INTEGER(4):: gelev                                  = NA32     ! 41-44 
     INTEGER(4):: selev                                  = NA32     ! 45-48 
     INTEGER(4):: sdepth                                 = NA32     ! 49-52 
     INTEGER(4):: gdel                                   = NA32     ! 53-56 
     INTEGER(4):: sdel                                   = NA32     ! 57-60 
     INTEGER(4):: swdep                                  = NA32     ! 61-64 
     INTEGER(4):: gwdep                                  = NA32     ! 65-68 
     INTEGER(2):: scalel                                 = NA16     ! 69-70 
     INTEGER(2):: scalco                                 = NA16     ! 71-72 
     INTEGER(4):: sx                                     = NA32     ! 73-76 
     INTEGER(4):: sy                                     = NA32     ! 77-80 
     INTEGER(4):: gx                                     = NA32     ! 81-84 
     INTEGER(4):: gy                                     = NA32     ! 85-88 
     INTEGER(2):: counit                                 = NA16     ! 89-90 
     INTEGER(2):: wevel                                  = NA16     ! 91-92 
     INTEGER(2):: swevel                                 = NA16     ! 93-94 
     INTEGER(2):: sut                                    = NA16     ! 95-96 
     INTEGER(2):: gut                                    = NA16     ! 97-98 
     INTEGER(2):: sstat                                  = NA16     ! 99-100 
     INTEGER(2):: gstat                                  = NA16     !101-102 
     INTEGER(2):: tstat                                  = NA16     !103-104 
     INTEGER(2):: laga                                   = na16     !105-106 
     INTEGER(2):: lagb                                   = na16     !107-108 
     INTEGER(2):: delrt                                  = na16     !109-110 
     INTEGER(2):: muts                                   = na16     !111-112 
     INTEGER(2):: mute                                   = na16     !113-114 
     INTEGER(2):: ns                                     = na16     !115-116 
     INTEGER(2):: dt                                     = na16     !117-118 
     INTEGER(2):: gain                                   = na16     !119-120 
     INTEGER(2):: igc                                    = na16     !121-122 
     INTEGER(2):: igi                                    = na16     !123-124 
     INTEGER(2):: corr                                   = na16     !125-126 
     INTEGER(2):: sfs                                    = na16     !127-128 
     INTEGER(2):: sfe                                    = na16     !129-130 
     INTEGER(2):: slen                                   = na16     !131-132 
     INTEGER(2):: styp                                   = na16     !133-134 
     INTEGER(2):: stas                                   = na16     !135-136 
     INTEGER(2):: stae                                   = na16     !137-138 
     INTEGER(2):: tatyp                                  = na16     !139-140 
     INTEGER(2):: afilf                                  = na16     !141-142 
     INTEGER(2):: afils                                  = na16     !143-144 
     INTEGER(2):: nofilf                                 = na16     !145-146 
     INTEGER(2):: nofils                                 = na16     !147-148 
     INTEGER(2):: lcf                                    = na16     !149-150 
     INTEGER(2):: hcf                                    = na16     !151-152 
     INTEGER(2):: lcs                                    = na16     !153-154 
     INTEGER(2):: hcs                                    = na16     !155-156 
     INTEGER(2):: year                                   = na16     !157-158 
     INTEGER(2):: day                                    = na16     !159-160 
     INTEGER(2):: hour                                   = na16     !161-162 
     INTEGER(2):: minute                                 = na16     !163-164 
     INTEGER(2):: sec                                    = na16     !165-166 
     INTEGER(2):: timbas                                 = na16     !167-168 
     INTEGER(2):: trwf                                   = na16     !169-170 
     INTEGER(2):: grnors                                 = na16     !171-172 
     INTEGER(2):: grnofr                                 = na16     !173-174 
     INTEGER(2):: grnlof                                 = na16     !175-176 
     INTEGER(2):: gaps                                   = na16     !177-178 
     INTEGER(2):: otrav                                  = na16     !179-180 
     REAL(4)   :: d1                                     = nr32     !181-184 
     REAL(4)   :: f1                                     = nr32     !185-188 
     REAL(4)   :: d2                                     = nr32     !189-192 
     REAL(4)   :: f2                                     = nr32     !193-196 
     REAL(4)   :: ungpow                                 = nr32     !197-200
     REAL(4)   :: unscale                                = nr32     !201-204 
     INTEGER(4):: ntr                                    = na32     !205-208 
     INTEGER(2):: mark                                   = na16     !209-210 
     INTEGER(2):: shortpad                               = na16     !211-212 
     INTEGER(2),DIMENSION(14) :: unassigned              = na16     !213-240 unassigned 

  END TYPE SU_TraceHeader

  TYPE su_trace
    TYPE(su_traceheader)     :: trcheader
    REAL(4),ALLOCATABLE :: trcdata(:)
END TYPE su_trace

END MODULE su_types_mod
