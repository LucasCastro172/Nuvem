MODULE su_io_mod 
! 
! SU format I/O using random acess
! to read and write traces
! 
! References: 
! 
! References: 
! Seismic Un*x Manual (http://www.cwp.mines.edu/cwpcodes/) 
! 
! 
! Author : Jesse Costa  (UFPA)
! e-mail : jesse@ufpa.br 
!
! WARNINGS:
! This module assumes byte record length, as does gfortran.
! When using ifort compiler make sure the use 
! "-assume byterecl" in the compiler options. 
!
! Assumes regular data, same number of samples for each trace 
! 
! Does not work if SU is installed with XDR option on
! assumes data is stored in machine native data type
!
USE SU_Types_mod 
IMPLICIT NONE 

PRIVATE
 
CHARACTER(1) :: BYTE
INTEGER :: BYTELEN,RECLEN,NS

interface GetSUtrace
   module procedure getSUtraceSplit
   module procedure getSUtraceType
end interface

interface PutSUtrace
   module procedure PutSUtraceSplit
   module procedure PutSUtraceType
end interface

PUBLIC::                   &
         OpenSU          , &
         GetSUtrace      , &
         PutSUtrace                  

CONTAINS 
  

FUNCTION OpenSU(SUid,SUFname,ioaction,nsamples) RESULT(ierr)
!============================================================
!
! Open an SU format file
!
!============================================================
INTEGER,      INTENT(IN) :: SUid
CHARACTER*(*),INTENT(IN)    :: SUFname
CHARACTER(1), INTENT(IN)    :: ioaction
INTEGER,      INTENT(INOUT) :: nsamples

INTEGER :: IERR
TYPE(SU_TraceHeader)::TraceHeader
character(512)::iomesg

INQUIRE (IOLENGTH = BYTELEN) BYTE

SELECT CASE(ioaction)

CASE('r')

     RECLEN = 240*BYTELEN
     OPEN(                                  &
           UNIT       = SUid              , & 
           FILE       = SUFname           , &
           IOSTAT     = IERR              , & 
           STATUS     ='old'              , & 
           FORM       ='unformatted'      , & 
           ACCESS     ='direct'           , & 
           RECL       = RECLEN            , &
           IOMSG      = iomesg            , &
           ACTION     ='read'               & 
          );
      if ( ierr /= 0 ) then 
              print *, trim(iomesg)
              return
      end if
 
     ierr=GetSUFirstTraceHeader(suid,TraceHeader)
     if (ierr /= 0) return

     ns       = TraceHeader%ns
     nsamples = ns

     CLOSE(SUid)

     RECLEN = (240 + 4*ns)*BYTELEN
          OPEN(                             &
           UNIT       = SUid              , & 
           FILE       = SUFname           , &
           IOSTAT     = IERR              , & 
           STATUS     ='old'              , & 
           FORM       ='unformatted'      , & 
           ACCESS     ='direct'           , & 
           RECL       = RECLEN            , &
           IOMSG      = iomesg            , &
           ACTION     ='read'               & 
          ); 

      if ( ierr /= 0 ) then 
              print *, trim(iomesg)
              return
      end if
CASE('a')

     RECLEN = 240*BYTELEN
     OPEN(                                  &
           UNIT       = SUid              , & 
           FILE       = SUFname           , &
           IOSTAT     = IERR              , & 
           STATUS     ='old'              , & 
           FORM       ='unformatted'      , & 
           ACCESS     ='direct'           , & 
           RECL       = RECLEN            , &
           IOMSG      = iomesg            , &
           ACTION     ='read'               & 
          );
      if ( ierr /= 0 ) then 
              print *, trim(iomesg)
              return
      end if
 
     ierr=GetSUFirstTraceHeader(suid,TraceHeader)
     if (ierr /= 0) return

     ns       = TraceHeader%ns
     nsamples = ns

     CLOSE(SUid)

     RECLEN = (240 + 4*ns)*BYTELEN
          OPEN(                             &
           UNIT       = SUid              , & 
           FILE       = SUFname           , &
           IOSTAT     = IERR              , & 
           STATUS     ='old'              , & 
           FORM       ='unformatted'      , & 
           ACCESS     ='direct'           , & 
           RECL       = RECLEN            , &
           IOMSG      = iomesg            , &
           ACTION     ='write'              & 
          ); 

      if ( ierr /= 0 ) then 
              print *, trim(iomesg)
              return
      end if


CASE('w')

     ns     = nsamples
     RECLEN = (240 + 4*ns)*BYTELEN
     OPEN(                                 &
           UNIT       = SUid             , & 
           FILE       = SUFname          , &
           IOSTAT     = IERR             , & 
           STATUS     = 'replace'        , & 
           FORM       = 'unformatted'    , & 
           ACCESS     = 'direct    '     , &
           RECL       = RECLEN           , &
           IOMSG      = iomesg            , &
           ACTION     = 'write'            & 
          );                     
        
      if ( ierr /= 0 ) then 
              print *, trim(iomesg)
              return
      end if
END SELECT
 
END FUNCTION OpenSU


FUNCTION GetSUFirstTraceHeader(suid,TraceHeader)  RESULT(ierr)
!========================================================== 
!
!  READ  TraceHeader 
!
!========================================================== 
INTEGER,INTENT(IN) :: suid 
TYPE(SU_TraceHeader),INTENT(INOUT)       :: TraceHeader 
 
INTEGER :: ierr 
! 
!  READ SU TraceHeader  	 	 	
! fdmig_teste.par
INQUIRE(UNIT=suid,IOSTAT=IERR)
IF(IERR == 0) THEN
READ(suid,REC=1,ERR=1000,IOSTAT=ierr) TraceHeader
END IF
1000 RETURN
 
END FUNCTION GetSUFirstTraceHeader 


FUNCTION GetSUTraceSplit(suid,trcrec,TraceHeader,sudata) RESULT(ierr)
!========================================================== 
!
!  READ  SU Trace 
!
!========================================================== 
INTEGER,INTENT(IN) :: suid
INTEGER,INTENT(IN) :: trcrec
TYPE(SU_TraceHeader),INTENT(OUT)  :: TraceHeader
REAL(4),INTENT(OUT),DIMENSION(ns) :: sudata

INTEGER::IERR

INQUIRE(UNIT=suid,IOSTAT=IERR)

IF(IERR == 0) THEN
   READ(suid,rec=trcrec,ERR=1000,IOSTAT=ierr) TraceHeader,sudata(1:ns)
END IF

1000 RETURN

END FUNCTION GetSUTraceSplit

FUNCTION GetSUTraceType(suid,trcrec,sutrace) RESULT(ierr)
!========================================================== 
!
!  READ  SU Trace 
!
!========================================================== 
INTEGER,INTENT(IN) :: suid
INTEGER,INTENT(IN) :: trcrec
TYPE(SU_Trace),INTENT(INOUT)  :: sutrace

INTEGER::IERR

INQUIRE(UNIT=suid,IOSTAT=IERR)

IF(IERR == 0) THEN
   READ(suid,rec=trcrec,ERR=1000,IOSTAT=ierr) sutrace%trcheader,sutrace%trcdata(1:ns)
END IF

1000 RETURN

END FUNCTION GetSUTraceType

FUNCTION PutSUTraceSplit(suid,trcrec,TraceHeader,sudata) RESULT(ierr)
INTEGER,INTENT(IN) :: suid
INTEGER,INTENT(IN) ::trcrec
TYPE(SU_TraceHeader),INTENT(IN) :: TraceHeader
REAL(4),INTENT(IN),DIMENSION(ns) :: sudata

INTEGER::IERR

INQUIRE(UNIT=suid,IOSTAT=IERR)

IF(IERR == 0) THEN
   WRITE(suid,REC=trcrec,ERR=1000,IOSTAT=ierr) TraceHeader,sudata
END IF

1000 RETURN

END FUNCTION PutSUTraceSplit

FUNCTION PutSUTraceType(suid,trcrec,sutrace) RESULT(ierr)
INTEGER,INTENT(IN) :: suid
INTEGER,INTENT(IN) ::trcrec
TYPE(SU_Trace),INTENT(IN) :: sutrace

INTEGER::IERR

INQUIRE(UNIT=suid,IOSTAT=IERR)

IF(IERR == 0) THEN
   WRITE(suid,REC=trcrec,ERR=1000,IOSTAT=ierr) sutrace%trcheader,sutrace%trcdata(1:ns)
END IF

1000 RETURN

END FUNCTION PutSUTraceType
 
END MODULE su_io_mod
